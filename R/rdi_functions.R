# geneCounts = list(v=c(IGHV=65,TRAV=43,TRBV=52),
#                   j=c(IGHJ=6, TRAJ=65,TRBJ=14))



############# MAIN FUNCTIONS ###############

#' Calculate RDI dissimilarity matrix 
#'
#' @name rdi
#' @description Wrapper function for calculating RDIs
#' 
#' @param genes         vector (or matrix) containing the gene calls for each sequence.
#'                      If genes is a matrix, counts will be calculated for each 
#'                      column of 'genes', and the resulting count matrices will be
#'                      concatenated.
#' @param seqAnnot      matrix containing repertoire annotations. Must be same length
#'                      as 'genes'.
#' @param params        list; contains parameters to pass to child functions.
#'                      Should contain \code{countParams} and \code{distParams} lists,
#'                      which contain parameters for \link{calcVDJcounts} and
#'                      \link{calcRDI}, respectively. See Details.
#' @param ...           other parameters to pass to \link{calcVDJcounts} and 
#'                      \link{calcRDI}.
#'                      
#' @details 
#' This function is a wrapper for the two core functions of RDI, 
#' \link{calcVDJcounts} and \link{calcRDI}. To control the function of 
#' both \code{calcVDJcounts} and \code{calcRDI}, additional parameters can be specified
#' either directly in the RDI function call, or parameters for the individual functions 
#' can be wrapped up into lists of parameters and passed into the \code{params} parameter.
#' \code{params} should be a list containing at least one of two parameter lists:
#' \code{countParams} and \code{distParams}, which  
#' will be passed to \code{calcVDJcounts} and \code{calcRDI}, respectively. An example
#' analysis is included below. 
#' 
#' 
#' @return 
#' Dissimilarity structure, as calculated by dist. In addition to the standard
#' attributes returned by dist, two additional attributes are defined as follows:
#' \tabular{rl}{
#'   \emph{nseq}  \tab integer, the number of sequences used after subsampling  
#'                               the repertoires\cr
#'   \emph{ngenes}  \tab integers, the number of genes in each column of "genes" 
#'                                 that were included in at least one repertoire.
#' }
#' 
#' @examples
#' 
#' #create genes
#' genes = sample(letters, 10000, replace=TRUE)
#' 
#' #create sequence annotations
#' seqAnnot = data.frame(donor = sample(1:4, 10000, replace=TRUE),
#'                       visit = sample(c("V1","V2","V3"), 10000, replace=TRUE),
#'                       cellType = sample(c("B","T"), 10000, replace=TRUE)
#'                      )
#'                      
#' #parameters
#' params = list(
#'   countParams = list(
#'     select = list(
#'       visit = c("V1","V3"),
#'       cellType = "B"
#'     ),
#'     combine = list(
#'       visit = "V[13]"
#'     ),
#'     simplifyNames = FALSE
#'   ),
#'   distParams = list(
#'     constScale=FALSE
#'   )
#' )
#'
#' ##calculate RDI
#' d = rdi(genes, seqAnnot, params)
#' 
#' ##plot using hierarchical clustering
#' plot(hclust(d))
#' 
#' @export
rdi <- function(genes, seqAnnot, params=NULL, ...){

  ##get VDJ counts for each leaf
  vdjCounts = do.call(calcVDJcounts, 
                      c(list(genes=genes, seqAnnot=seqAnnot),
                        params$countParams,list(...)))
  
  ##calculate RDI matrix
  distMat = do.call(calcRDI,
                    c(list(vdjCounts=vdjCounts), 
                      params$distParams,list(...)))
  
  return(distMat)
}

#' Calculate repertoire counts
#' 
#' @name calcVDJcounts
#' @description Create count matrices for a set of repertoires
#'
#' @param   genes           vector (or matrix) containing the gene calls for each sequence.
#'                          If genes is a matrix, counts will be calculated for each 
#'                          column of 'genes', and the resulting count matrices will be
#'                          concatenated. See Details.
#' @param   seqAnnot        matrix containing repertoire annotations.
#' @param   select          a list containing definitions of which repertoires to use. 
#'                          See Details.
#' @param   combine         a list defining repertoires that should be grouped together.
#'                          See Details. 
#' @param   vdjDrop         a list specifying specific genes to exclude from analysis.
#'                          See Details.
#' @param   splitBy         the columns in seqAnnot to use for splitting repertoires. 
#'                          Default is to use all columns in seqAnnot.
#' @param   simplifyNames   logical; if true, any columns of seqAnnot where all selected
#'                          (and collapsed) sequences share the same value will not be 
#'                          used to make the names of sequenceCodes.
#' @param   splitCommas     logical; if true, seqAnnot is assumed to contain.
#'                          comma-separated lists of sequence annotations, which will be 
#'                          split up before generating sequence codes. Note: setting this 
#'                          to TRUE will make processing much slower.
#' @param   ...             additional parameters; these are ignored by this function.
#'
#' @details 
#' In most cases, \code{genes} will be a single vector or one-column matrix. However, 
#' there are some cases where a row of \code{seqAnnot} corresponds to two (or more) genes
#' (e.g. the V and J gene segments of a single immune sequence). Rather than make multiple
#' rows for each gene, the \code{calcVDJcounts} function provides the option to provide 
#' a multi-column matrix for \code{genes}. The counts for each column will be tallied 
#' separately, and are then concatenated.
#' 
#' To ensure equal variance across all repertoires, the default RDI metric uses 
#' subsampling to ensure that all repertoires have the same number of sequences. The 
#' default RDI metric subsamples all repertoires to the size of the smallest repertoire, 
#' which may result in a loss of power for comparisons between larger repertoires.
#' In order to increase power for various tests, it is often useful to only calculate 
#' the repertoire counts for a subset of the repertoires in seqAnnot. This can be done by
#' using the \code{select} and \code{combine} parameters to specify which 
#' repertoires to include in the analysis.
#' 
#' Both parameters are lists containing entries
#' with the same name as one of the columns of seqAnnot. For \code{select}, each entry is 
#' a vector defining which values to include (e.g., to include only Visit 1 and 3, you 
#' might specify \code{select=list(visit=c("V1","V3"))}, where the \code{'visit'} column 
#' in seqAnnot contains the values \code{"V1"},\code{"V2"}, and \code{"V3"}). In this 
#' case, any rows of \code{genes} and \code{seqAnnot} that come from a repertoire not 
#' specified in \code{select} will be discarded. By default, if a \code{select} code is 
#' not specified for a column in \code{seqAnnot}, all values from that column will be 
#' included.
#' 
#' The \code{combine} parameter works in a similar fashion, but instead of a vector
#' describing which parameters to include, you can specify a vector of regular 
#' expressions, and any values of the \code{seqAnnot} column that match the regular
#' expression will be combined into a single repertoire (e.g. to combine visits 1 and 3 
#' into a single repertoire, you might specify \code{combine=list(visit="V[13]")}). 
#' 
#' The \code{vdjDrop} parameter is also useful for limiting sequences. Like 
#' \code{select} and \code{combine}, this is a named list, with entries corresponding to 
#' the columns of \code{genes}. Each entry of \code{vdjDrop} is a vector of gene segment
#' names to remove from the analysis. All sequences containing those genes are removed 
#' from the analysis before subsampling.
#'
#' Once unwanted rows have been removed, the columns of \code{seqAnnot} are concatenated 
#' to generate "repertoire" labels for each row. The repertoire labels are then used
#' to split the rows of \code{genes}, and gene prevalence is tallied within a repertoire.
#' By default, columns of \code{seqAnnot} that are constant after subsetting will not be
#' included in the label. However, this can be controlled by the \code{simplifyNames} 
#' parameter. If \code{simplifyNames} is FALSE, all columns of  \code{seqAnnot} are 
#' included when generating labels.
#'
#' @return
#' A matrix where each row represents a gene, and each column represents a 
#' repertoire.
#' 
#' 
#' @examples 
#' 
#' #create genes
#' genes = sample(letters, 10000, replace=TRUE)
#' 
#' #create sequence annotations
#' seqAnnot = data.frame(donor = sample(1:4, 10000, replace=TRUE),
#'                       visit = sample(c("V1","V2","V3"), 10000, replace=TRUE),
#'                       cellType = sample(c("B","T"), 10000, replace=TRUE)
#'                      )
#'                      
#' ##generate repertoire counts for all repertoires
#' cts = calcVDJcounts(genes,seqAnnot) 
#' 
#' ##Only include visit 1
#' cts = calcVDJcounts(genes,seqAnnot, select=list(visit="V1"))
#' 
#' 
#' ## Just T cell repertoires, combining visit 1 and 3 together, and dropping visit 2
#' cts = calcVDJcounts(genes,seqAnnot, 
#'                     select=list(cellType="T", visit=c("V1","V3")), 
#'                     combine=list(visit="V[13]")) 
#' 
#' @export
calcVDJcounts = function(genes, seqAnnot, select=NULL, combine=NULL, vdjDrop=NULL,
                         splitBy=NULL, simplifyNames=TRUE, splitCommas=FALSE, ...){
  ### Check parameter formatting
  ##coerce to data frames
  genes = as.data.frame(genes, stringsAsFactors = F)
  if(!is.data.frame(genes)){stop("cannot coerce genes to class data.frame")}
  seqAnnot = as.data.frame(seqAnnot, stringsAsFactors = F)
  if(!is.data.frame(seqAnnot)){stop("cannot coerce seqAnnot to class data.frame")}
  
  ##check that genes matches seqAnnot
  if(nrow(genes)!=nrow(seqAnnot)){stop("size of genes array does not match seqAnnot")}
  
  ##generate sequence codes
  codeData = createSequenceCodes(genes, seqAnnot, select,combine, vdjDrop,
                                 splitBy,simplifyNames, splitCommas)
  genes.sub = codeData$genes
  sequenceCode = codeData$sequenceCode
  
  
  ##generate table
  vdjCounts = lapply(1:ncol(genes.sub),
                     function(i) { table(genes.sub[, i], sequenceCode) })
  vdjCounts = do.call(rbind, vdjCounts)
  
  attr(vdjCounts, "ngenes") = apply(genes.sub, 2, 
                            function(x) { x = unique(x); sum(!is.na(x)) })

  ##edit labels -- remove extra spaces, etc
  labs = colnames(vdjCounts)

  labs = gsub(" +"," ", labs)    ##get rid of extra spaces
  labs = gsub("(^ | $)","",labs) ##get rid of leading and trailing spaces
  colnames(vdjCounts) = labs
  
  return(vdjCounts) 
}

# Create sequence codes for use in calcVDJcounts
#
# @param   genes          vector (or matrix) containing the gene calls for each sequence.
#                         If genes is a matrix, counts will be calculated for each 
#                         column of 'genes', and the resulting count matrices will be
#                         concatenated.
# @param   seqAnnot       matrix containing repertoire annotations for each sequence.
#                         Must be the same lenght as genes. 
# @param   select         a list containing definitions of which repertoires to use. 
#                         See Details.
# @param   combine        a list defining repertoires that should be grouped together.
# @param   vdjDrop        a list specifying specific genes to exclude from analysis.
# @param   splitBy        the columns in seqAnnot to use for splitting repertoires
# @param   simplifyNames  logical; if true, any columns of seqAnnot where the selected
#                         (and collapsed) sequences share the same value will not be 
#                         used to make the names of sequenceCodes.
# @param   splitCommas    logical; if true, seqAnnot is assumed to contain
#                         comma-separated lists of sequence annotations, which will be 
#                         split up before generating sequence codes.
# 
# @return
# A list containing two objects: genes and sequenceCodes. Genes is a subset of
# the input genes matrix containing only rows that were selected by "select".
# sequenceCode contains a vector of repertoire IDs for each sequence in genes.
createSequenceCodes = function(genes, seqAnnot, select=NULL, combine=NULL,vdjDrop=NULL,  
                               splitBy=NULL, simplifyNames=TRUE,splitCommas=FALSE){
  ##validate inputs
  if(is.null(splitBy)){splitBy=colnames(seqAnnot)}
  
  ########
  ##generate identifying sample codes for each read
  
  ##create codes for each distinguishing column
  sequenceCodeCols = lapply(splitBy, function(n){
    code = select[[n]]
    x = as.character(seqAnnot[,n])
    
    if(splitCommas){
      #split codes into the comma-separated values
      x = strsplit(x,",")
      i.orig = rep(1:length(x), sapply(x,length))
      x = unlist(x)
      names(x) = i.orig
    }else{
      names(x) = 1:length(x)
    }
    
    ##drop everything not specified by "code"  (this takes a while)
    if(!is.null(code)){
      keep = x %in% code
      x=x[keep]
    }
    
    ##combine together codes from this column
    if(!is.null(combine[[n]])){
      ##get the combineNames
      comb.name = sapply(combine[[n]],regexToName)
      ##for each thing to combine, set it to 
      for(i in 1:length(combine[[n]])){
        simplify = grepl(combine[[n]][i],x)
        x[simplify] = comb.name[i]
      }
    }
    
    if(splitCommas){
      ##drop seqs that are replicates (created by the combine step)
      isDuplicated = duplicated(data.frame(x,names(x)))
      x = x[!isDuplicated]
      x = split(x,names(x))
    }
    
    x
  })
  inAll = Reduce(intersect,lapply(sequenceCodeCols, names))
  
  ##check if we get any sequences (if not, return an error)
  if(length(inAll)==0){
    i = names(which(sapply(sequenceCodeCols,length)==0))
    if(length(i>0)){
      stop("No sequences found matching the following input codes: ",unlist(select[i]))
    }
    stop("No sequences found matching the input codes (codes are mutually exclusive)");
  }
  
  ##get rid of all sequences that were masked or removed from one of the columns
  sequenceCodeCols = lapply(sequenceCodeCols,function(x){x[inAll]})
  
  ##mask genes specified by vdjDrop
  if(!is.null(vdjDrop)){
    for(dr in vdjDrop)
      genes = apply(genes,2,function(x){x[grep(dr,x)] = NA;x})
  }
  
  ##for labeling purposes, remove any parts of sequenceCodeCols that are constant
  if(simplifyNames){
    codeParts = lapply(sequenceCodeCols, function(x){unique(unlist(x))})
    isUnique = sapply(codeParts,length)==1
    sequenceCodeCols = sequenceCodeCols[!isUnique]
  }
  
  if(splitCommas){
    if(length(sequenceCodeCols)>1){
      ##create sample codes by pasting together sequenceCodeCols 
      #  and expanding out cells with more than one label
      lens = lapply(sequenceCodeCols, function(x){sapply(x,length)})
      sequenceCodeCols.x = lapply(1:length(lens), function(i){
        toX = Reduce("*",lens[-i])
        rep(sequenceCodeCols[[i]], toX)
      })
      ids = unlist(lapply(sequenceCodeCols.x[[1]],names))
      sequenceCode=do.call(paste, lapply(sequenceCodeCols.x,unlist))
    
    }else{ ##if sequenceCodeCols is length 1, just unlist 
      ids = unlist(lapply(sequenceCodeCols[[1]],names))
      sequenceCode = unlist(sequenceCodeCols[[1]])
    }
    
    ##make sure "genes" matches up
    genes.m = genes[as.numeric(ids),,drop=F]
  }else{
    ##paste together sequenceCodeCols
    sequenceCode=do.call(paste, sequenceCodeCols)
    genes.m = genes[as.numeric(inAll),,drop=F]
  }  
  
  list(genes=genes.m, sequenceCode=sequenceCode)
}


# Transform and scale VDJ counts
# 
# @param vdjCounts   a matrix of vdj counts, as created by calcVDJCounts
# @param constScale  logical; if \code{TRUE}, vdjCounts will be scaled such that the sum of
#                    each column will be equal to 500 (an arbitrary constant). Otherwise,
#                    the columns will be scaled to the average count of all the columns.
# @param xform       The function to use for transforming the VDJ counts. Default is the
#                    arcsinh transformation.
# 
# @return The transformed and normalized VDJ counts
transformVDJCounts = function(vdjCounts, constScale=TRUE, xform=asinh_xform){
  #normalize by total counts
  normVDJCounts = t(vdjCounts)/colSums(vdjCounts)
  
  ct.mu = mean(colSums(vdjCounts))
  ##if not making a hexplot, put onto same scale
  if(constScale){
    ct.mu = 500
  }
  
  apply(normVDJCounts*ct.mu,2,xform)
}


#' Calculate repertoire distances
#' 
#' @name calcRDI
#' @description Calculate repertoire distances from a matrix of vdjCounts
#' 
#' @param vdjCounts     a matrix of repertoire counts, as created by calcVDJCounts
#' @param distMethod    one of c("euclidean","cor") determining how to calculate the 
#'                      distance from the matrix of vdj counts. See Details.
#' @param subsample     logical; if true, all repertoires are subsampled to be equal 
#'                      size. 
#' @param nIter        value defining how many times the subsampling should be 
#'                      repeated. Only used if subsample is TRUE.
#' @param constScale    logical; if \code{TRUE}, vdjCounts will be scaled such that the sum of
#'                      each column will be equal to 500 counts (an arbitrary constant). 
#'                      Otherwise, the columns will be scaled to the average count of all 
#'                      the columns.
#' @param units         One of "lfc" or "pct". This determines the method used for 
#'                      transforming the repertoire counts. See Details.
#' @param   ...         additional parameters; these are ignored by this function.    
#'                    
#' @details  
#'  There are two options for distance methods, "euclidean" and "cor". Euclidean refers to
#'  standard euclidean distance, and is the standard for the RDI measure described in 
#'  (Bolen et al. Bioinformatics 2016). In contrast, cor refers to a correlation-based 
#'  distance metric, where the distance is defined as \code{(1-correlation)} between each 
#'  column of vdjCounts. 
#'  
#'  The \code{units} parameter is used to determine the transformation function for the
#'  repertoire counts. If \code{units='lfc'} (default), then the arcsinh transformation 
#'  is applied to the count matrix, resulting in a distance metric which 
#'  will scale with the average log fold change of each gene. In contrast,
#'  \code{units='pct'} will result in no transformation of the count matrix, and distances
#'  will be proportional to the average percent change of each gene, instead. Note that
#'  "units" is a bit of a misnomer, as the distance metric doesn't actually represent the
#'  true log-fold or percent change in the repertoires. In order to actually estimate 
#'  these parameters, refer to the \link{rdiModel} and \link{convertRDI}
#'  functions.
#'  
#' @return
#' A dissimilarity structure containing distances between repertoires, averaged
#' across each subsampe run.
#' In addition to the standard attributes in a dist object, 
#' three additional attributes are defined as follows:
#' \tabular{rl}{
#'   \emph{ngenes}  \tab integers, the number of genes in each column of "genes" 
#'                                 that were included in at least one repertoire.\cr
#'   \emph{nseq}  \tab integer, the number of sequences used after subsampling  
#'                               the repertoires. If \code{subsample=FALSE}, this is not
#'                               defined.\cr
#'   \emph{units} \tab string, either "lfc" or "pct", depending on the "units" in the 
#'                     original call
#' }
#' 
#' @examples 
#' 
#' #create genes
#' genes = sample(letters, 10000, replace=TRUE)
#' 
#' #create sequence annotations
#' seqAnnot = data.frame(donor = sample(1:4, 10000, replace=TRUE),
#'                       cellType = sample(c("B","T"), 10000, replace=TRUE)
#'                      )
#' ##generate repertoire counts
#' cts = calcVDJcounts(genes,seqAnnot) 
#' 
#' ##calculate RDI 
#' d = calcRDI(cts)
#' 
#' ##calculate RDI in percent space
#' d_pct = calcRDI(cts,units="pct")
#' 
#' ##convert RDI to actual 'lfc' estimates and compare
#' dtrue = convertRDI(d)$pred
#' plot(d, dtrue)
#' 
#' @export
calcRDI = function(vdjCounts, distMethod=c("euclidean", "cor"), 
                   subsample=TRUE, nIter=100, constScale=TRUE, 
                   units=c("lfc","pct"), ...){
  ##match arguments
  distMethod = match.arg(distMethod)
  
  units = match.arg(units)
  if(units == "lfc"){xform=asinh_xform}
  else {xform=base::identity}
  
  ###############
  ## frequency code
  #if(distMethod=="euclidean" || distMethod=="cor"){
    
    ##determine distance function
    dfunc = stats::dist
    if(distMethod=="cor"){dfunc = function(x){as.dist(cor(t(x)))}}
    
    if(subsample){
    
      ##get number of sequences in each repertoire
      allSizes = colSums(vdjCounts)
      ##find the smallest sized repertoire
      size = min(allSizes)
      
      distMat = matrix(0, ncol(vdjCounts),ncol(vdjCounts))
      colnames(distMat) = colnames(vdjCounts)
      distMat = as.dist(distMat)
      ##sum distMats
      for(i in 1:nIter){
        
        ##subsample vdjCounts
        vdjCounts.pt = sapply(1:length(allSizes), function(j){
          x = rep(0, nrow(vdjCounts))
          for(i in 1:nrow(vdjCounts)){
            x[i] = rhyper(1, 
                          vdjCounts[i,j], 
                          allSizes[j]-sum(vdjCounts[1:i,j]), 
                          size-sum(x))
          }
          x
        })
        
        ##transform counts
        normVDJCounts.trans = transformVDJCounts(vdjCounts.pt,
                                                 constScale=constScale,xform=xform)
        
        ####calc distance -- average subsample runs together 
        tryCatch({
            distMat = distMat + dfunc(normVDJCounts.trans)
          },
          ##I'm silencing the "standard deviation is zero" warning, because we don't care.
          warning=function(w){if(!grepl("standard deviation is zero",w)){warning(w)}}
        )
      }
      ##and divide by number of subsamples
      distMat = distMat/nIter
      attr(distMat,"nseq") = size
    }else{
      ##transform counts
      normVDJCounts.trans = transformVDJCounts(vdjCounts,
                                               constScale=constScale,xform=xform)
      distMat = dist(normVDJCounts.trans)
      
    }
  
  #}
  
  ##########
  ## Jaccard distance code
  ##compare types with more than 1 sequence for each patient
#   if(distMethod=="jaccard"){
#     foundRecomb = lapply(1:ncol(vdjCounts), function(i){
#       x = rownames(vdjCounts)[vdjCounts[,i]>1]
#       #if(length(x)>100){return(x[order(vdjCounts[vdjCounts[,i]>1,i],decreasing=T)][1:100])}
#       return(x)
#     })
#     names(foundRecomb) = colnames(vdjCounts)
#     
#     jaccardMat = matrix(0, length(foundRecomb),length(foundRecomb))
#     for(i in 1:length(foundRecomb)){
#       for(j in 1:length(foundRecomb)){
#         a = foundRecomb[[i]]; b = foundRecomb[[j]]
#         jac = length(intersect(a,b))/length(union(a,b))
#         if(length(a)==0 & length(b)==0){jac=0}
#         jaccardMat[i,j] = jac
#       }
#     }
#     distMat = as.dist(1-jaccardMat)
#   }
  
  ##carry over the ngenes attribute if it exists
  if(is.null(attr(vdjCounts,"ngenes"))){
    attr(distMat,"ngenes") = nrow(vdjCounts)
  }else{
    attr(distMat,"ngenes") = attr(vdjCounts,"ngenes")
  }
  
  ##fix the "call" attribute
  attr(distMat, "call") = match.call()
  ##add units attribute
  attr(distMat, "units") = units
  
  distMat
}


#' Convert RDI measures
#' @name convertRDI
#' @description Method to convert RDI values to fold/percent change
#' 
#' @param d       Distance matrix (as produced by \link{calcRDI}), or a vector of 
#'                distances.
#' @param models  Set of RDI models, as produced by \link{rdiModel}. If \code{NULL},
#'                RDI models will be calculated based on the attributes in the distance
#'                matrix.
#' @param calcSD  logical; if \code{TRUE}, standard deviations for each estimate will be returned.
#' 
#' @details  
#' The convertRDI function works by first generating a model for the RDI values at a given
#' repertoire size and feature count using the \link{rdiModel} function (see that 
#' method's help file for more details). The RDI models predict the average 
#' log-fold/percent change across a range of RDI values, and allows us to convert RDI to
#' a more stable and interpretable metric.  
#' 
#' In addition to the average log-fold or percent change value, \link{rdiModel} 
#' also generates models for the standard deviation at each RDI value. This is useful for
#' understanding the confidence intervals around the fold change estimate. 
#'                
#' @return
#' A list containing either one or two features:
#' 
#' \tabular{rl}{
#'   \emph{pred}  \tab The converted predictions; same length as \code{d}. \cr
#'   \emph{sd}  \tab If \code{calcSD==T}, a set of standard deviation estimates for each
#'                prediction. 
#' }
#' 
#' @examples 
#' 
#' #create genes
#' genes = sample(letters, 10000, replace=TRUE)
#' #create sequence annotations
#' seqAnnot = data.frame(donor = sample(1:4, 10000, replace=TRUE))
#' #calculate RDI
#' d = rdi(genes, seqAnnot)
#' 
#' ##convert RDI to actual 'lfc' estimates and compare
#' dtrue = convertRDI(d)$pred
#' plot(d, dtrue)
#' 
#' ##look at SD ranges around lfc estimates
#' dtrue = convertRDI(d, calcSD=TRUE)
#' 
#' ##plot using ggplot2
#' library(ggplot2)
#' x = as.numeric(d)
#' y = as.numeric(dtrue$pred)
#' sd = as.numeric(dtrue$sd)
#' qplot(x,y)+geom_errorbar(aes(x=x, ymin=y-sd, ymax=y+sd))
#' 
#' @export
convertRDI <- function(d, models=NULL, calcSD=FALSE){
  if(is.null(models)){
    if(!("dist" %in% class(d)) || 
       !all(c("nseq","ngenes","units") %in% names(attributes(d))))
    {
      stop("Required attributes 'nseq', 'ngenes', or 'units' not found in distances. Cannot generate RDI models." )
    }
    models = rdiModel(attr(d,"nseq"),attr(d,"ngenes"),units=attr(d,"units"))
  }
  
  mean = predict(models$rev.fit, d)$y
  mean[mean<0] = 0
  if(calcSD){
    sd = predict(models$rev.fit$varFit, d)$y
    sd[sd<0] = 0
    return(list(pred=mean,sd=sd))
  }
  return(list(pred=mean))
  
}


#' RDI Models
#' @description Generate models equating RDI values to true differences in underlying 
#' prevalence values
#' 
#' @param n           the repertoire size.
#' @param ngenes      numeric vector indicating the number of genes in each chain.
#'                    If \code{baseVects} is not provided, this parameter will be used to 
#'                    generate a base prevalence vector for each of the genes.
#' @param baseVects   A vector or list of vectors representing the total prevalence of 
#'                    each gene (for each chain) in the dataset. 
#'                    Differential datasets will be created from alterations of
#'                    this vector. If not provided, a base vector will be randomly .
#'                    generated at each subsample step containing the number of genes
#'                    specified by ngenes.
#' @param nIter      The number of iterations (i.e. number of datasets to generate).
#' @param nSample     The number of samples to generate for each subsample. Each sample 
#'                    will have a different true fold change, but the same starting vector
#' @param units       String; either "lfc" or "pct", depending on what transform was used
#'                    in the original RDI calculations. See Details.
#' @param constScale  logical; if \code{TRUE}, vdjCounts will be scaled such that the sum of
#'                    each column will be equal to 500 (an arbitrary constant). Otherwise,
#'                    the columns will be scaled to the average count of all the columns.
#' 
#' @details 
#' This method uses simulated sequencing datasets to estimate the RDI values for datasets
#' with a known true deviation. 
#' 
#' Briefly, a baseline probability vector (either randomly generated or supplied by the
#' \code{baseVects} parameter) is randomly perturbed, and the difference between the 
#' baseline vector and the perturbed vector is calculated. Then, \code{nSample} sequencing
#' datasets of size n are randomly drawn from both the baseline vector and the perturbed 
#' vector, and the RDI distance between all datasets calculated. This process is repeated 
#' \code{nIter} times, resulting in a dataset of RDI values and matched true differences.
#' A set of spline models is then fit to the data: one from RDI to true difference, and 
#' another from true difference to RDI value, allowing for bi-directional conversions. 
#' 
#' If a baseline probability vector is not provided, one will be generated from an 
#' empirical model of gene segment prevalence. However, for best performance, this is not 
#' recommended. Estimates of true fold change is very sensitive to the distribution of 
#' features in your count dataset, and it is important that your baseline vector match
#' your overall dataset as accurately as possible. The best baseline vector is almost 
#' always the average feature prevalence across all repertoires in a dataset, although
#' manually generated baseline vectors may also work well.
#' 
#' The units used for the RDI model should always match the units used to generate your
#' RDI values. For more details on units, refer to the details of \link{calcRDI}.
#'
#' @return
#' A list containing three objects:
#' \tabular{rl}{
#'   \emph{fit}      \tab an object of class \code{"smooth.spline"}, based on a spline 
#'                   model with the true difference (lfc or pct) as the independent (x)
#'                   and RDI as the dependent (y). Used for converting from true difference
#'                   to RDI.  \cr
#'   \emph{rev.fit}  \tab an object of class \code{"smooth.spline"}. 
#'                   The opposite of fit. Used for converting from RDI to true difference.
#'                    \cr
#'   \emph{units}    \tab one of \code{c("lfc","pct")}, representing the units of the true
#'                   difference values.
#' }
#'
#' @seealso  \link{rdiAxis}, \link{rdiLadder}, \link{plotRDIladder}
#'
#' @examples 
#' #create genes
#' genes = sample(letters, 10000, replace=TRUE)
#' #create sequence annotations
#' seqAnnot = data.frame(donor = sample(1:4, 10000, replace=TRUE))
#' #calculate RDI
#' d = rdi(genes, seqAnnot)
#' 
#' ##create a "baseVect" with the same probability as our features
#' ##since we sampled uniformly, the base vector has equal probability
#' baseVect = rep(1/length(letters),length(letters))
#' 
#' ##generate an RDI model
#' m = rdiModel(attr(d, "nseq"), baseVects=baseVect)
#' 
#' ##plot the spline model
#' plot(m$fit, xlab="log fold change",ylab="RDI",type='l')
#' 
#' ##convert RDI to log fold change
#' mean = predict(m$rev.fit, d)$y
#' mean[mean<0] = 0
#'
#' @export
rdiModel = function(n, ngenes=NULL, baseVects=NULL, nIter=50, nSample=20,
                    units=c("lfc","pct"), constScale=TRUE){
  
  fcMeans = bootstrapRDI(n, ngenes, baseVects, nIter, nSample, 
                         units, constScale)
  
  ##get the true dists and calculate a spline model
  td_all = as.numeric(names(fcMeans))
  fit = smooth.spline(td_all, fcMeans, spar=1)
  rev.fit = smooth.spline(fcMeans,td_all,spar=1)
  
  ##calculate sds from residuals
  resid = (fit$yin - fit$y) / (1-fit$lev)  ##jackknife residuals
  rev.resid = (rev.fit$yin - rev.fit$y) / (1-rev.fit$lev)
  ##use splines to find point estimates of the variance across fcMeans
  ##note: I have no idea if this is reasonable, but this _is_ the exact formula for 
  ## the variance around a mean. In this case, the residuals are (mu-X), and spline
  ## finds the average of those points across the range of fcMeans
  var.fit = smooth.spline(fit$x, resid^2, spar=1)
  var.revfit = smooth.spline(rev.fit$x, rev.resid^2, spar=1)
  
  ##all estimates below the variance of the baseline are set to baseline
  var.fit$y[var.fit$y<var(fcMeans[td_all==0])] = var(fcMeans[td_all==0])
  var.fit$fit$coef[var.fit$fit$coef<var(fcMeans[td_all==0])] = var(fcMeans[td_all==0])
  
  ##add these sd vals to the fit
  fit$varFit = var.fit
  rev.fit$varFit = var.revfit
  
  out = list(fit=fit, rev.fit=rev.fit, units=units)
  class(out) = c(class(out),"rdiModel")
  return(out)
}


#' RDI ladder
#' 
#' @name rdiLadder
#' @description Function for creating the RDI ladder for a specific number of sequences
#' 
#' @param n           the repertoire size; alternatively, an rdiModel object as created by
#'                    \link{rdiModel}.
#' @param ngenes      numeric vector indicating the number of genes in each chain.
#'                    If baseVect is not provided, this parameter will be used to 
#'                    generate a base prevalence vector for each of the genes.
#' @param baseVects   A vector or list of vectors representing the total prevalence of 
#'                    each gene (for each chain) in the dataset. See \link{rdiModel} for details.
#' @param diffPoints  numeric vector; each value specifices either a log2 fold change or 
#'                    percent deviation value (depending on the 'units') at which 
#'                    the RDI ladder will be calculated.
#' @param units       String; either "lfc" or "pct", depending on what transform was used
#'                    in the original RDI calculations. See Details.
#' @param ...         Additional parameters to be passed to \link{rdiModel}
#' 
#' 
#' @details  
#' Because RDI values vary according to the number of genes and size of the repertoires,
#' they are not useful as numbers by themselves. Instead, it is useful to compare them 
#' with estimates of the true difference between the two repertoires. 
#' This function uses the models generated by \link{rdiModel} to generate estimated
#' RDI values corresponding to a set of pre-defined true distance (log-fold change or 
#' percent) values. This function is primarily meant to be used in conjunction with
#' \link{plotRDIladder} in order to add a useful reference point for RDI values. 
#' 
#' The units used for the RDI model should always match the units used to generate your
#' RDI values. For more details on units, refer to the details of \link{calcRDI}
#' 
#' @return 
#' A list of the same length as diffPoints, with each entry in the list containing
#' the mean RDI value and standard deviation corresponding to a given true difference value.
#' 
#' @seealso  \link{plotRDIladder}, \link{rdiModel}, \link{rdiAxis}
#' 
#' @export
rdiLadder = function(n, ngenes=NULL, baseVects=NULL, diffPoints=NULL, 
                     units=c("lfc","pct"), ...){
  ##check the units
  units = match.arg(units)
  
  ##get models
  if("rdiModel" %in% class(n)){
    models = n
  }else{
    models = rdiModel(n,ngenes,baseVects,units=units,...)
  }  
  
  ##define ladder points
  if(is.null(diffPoints)){
    if(units=="lfc"){
      diffPoints = log2(c(1,1.2,1.5,2,3,4))
      names(diffPoints) = paste0(c(1,1.2,1.5,2,3,4),"x")
    }else{
      diffPoints = pretty(range(models$fit$x))
      names(diffPoints) = paste0(diffPoints,"%")
    }
  }
  
  ##use models to calculate means and sds for each point
  ladder = lapply(diffPoints, function(fc){
    c(mean=predict(models$fit,fc)$y, sd=sqrt(predict(models$fit$varFit,fc)$y))
  })
  names(ladder) = diffPoints
  ladder
}


# The main function that bootstraps artifical data to create the RDI ladder/model fit
bootstrapRDI <- function(n,nGenes=NULL, baseVects=NULL, nIter=50, nSample=20,
                         units=c("lfc","pct"), constScale=TRUE){
  ##make sure one of nGenes or baseVect is specified
  if(is.null(nGenes) && is.null(baseVects)){
    stop("Must provide either nGenes or baseVects")
  }
  if(!is.null(baseVects)){
    if(!is.list(baseVects)){
      baseVects = list(baseVects)
    }
    baseVects = lapply(baseVects, function(x){x/sum(x)})
  }
  
  ##check units param
  units = match.arg(units)
  
  ##set the transform
  if(units=="lfc"){xform=asinh_xform
  }else{           xform=identity}
  
  ##set the seed fold/pct change for each subsample run.
  if(units=="lfc"){
    seed = 2^(runif(nIter, 0, 3))
  }else{
    seed = 2^(runif(nIter, 0, 10))
  }
  fn.sd = log2(seed)/sqrt(2/pi)
    
  ##iterate each bootstrap
  dists = lapply(1:nIter, function(bs){
    ##define the baseline vector
    if(is.null(baseVects)){
      ##generate a random base vector for this iteration
      baseVects = lapply(nGenes, rGene)
    }
    if(is.null(names(baseVects))){names(baseVects) = 1:length(baseVects)}
    
    ##get the baseline counts
    baseCts = NULL
    for(ns in 1:nSample){
      baseCts.pt = lapply(baseVects, function(truePcts){
        x = rep(0, length(truePcts))
        for(i in 1:length(x)){
          x[i] = rbinom(1, n-sum(x), truePcts[i]/sum(truePcts[i:length(truePcts)]))
        }
        x
      })
      baseCts.pt = as.matrix(unlist(baseCts.pt))
      baseCts = merge.table(baseCts, baseCts.pt)
    }
    ##normalize and transform base counts
    baseCts.x = transformVDJCounts(baseCts,constScale = constScale, xform=xform)
    
    ##get the fold-changed counts
    diffCts = NULL
    trueDiffs = NULL
    ##nSample number of samples for each subsample run
    for(ns in 1:nSample){
      diffInfo = lapply(baseVects,function(baseVect){
        truePcts = baseVect
        #create edited vector
        truePcts = 2^(log2(truePcts) + rnorm(length(truePcts), 0, fn.sd[bs]))
        truePcts = truePcts/sum(truePcts)
        
        ##define the difference score for this run
        ##either the average absolute lfc, or the linear distance between the two vectors
        if(units == "lfc"){ diff = mean(abs(log2(truePcts)-log2(baseVect))) }
        else{               diff = mean(abs(baseVect-truePcts))*100 }
        #else{              diff = dist(rbind(baseVect,truePcts)) }
        
        x = rep(0, length(truePcts))
        for(i in 1:length(x)){
          x[i] = rbinom(1, n-sum(x), truePcts[i]/sum(truePcts[i:length(truePcts)]))
        }
        list(cts=x,diff=diff)
      })
      diffCts.pt = as.matrix(unlist(lapply(diffInfo,"[[","cts")))
      diffCts = merge.table(diffCts, diffCts.pt)
      
      ##if there is more than one base vect, combine the trueDiff info together. 
      ##For "lfc"s, this is the weighted average of the log diffs
      ##For "pct"s, this is the weighted average of the squared diffs
      tds = sapply(diffInfo, "[[","diff")
      if(units=="lfc"){ 
        td = sum(log2(tds) * sapply(baseVects,length) ) / length(unlist(baseVects))
        td = 2^td
      }else{
        td = sum(  tds^2   * sapply(baseVects,length) ) / length(unlist(baseVects))     
        td = sqrt(td)
      }
      
      trueDiffs = c(trueDiffs, td)
    }
    ##normalize and transform counts
    diffCts.x = transformVDJCounts(diffCts,constScale = constScale, xform=xform)
    
    ##get the distance between the baseCounts
    baseDist = dist(baseCts.x)
    
    ##for the rest of the points, calculate it using matrix multiplication
    diffDist = pdist(baseCts.x,diffCts.x)@dist
    ##diffDists is a vector of all comparisons between baseCts and diffCts
    ##so we repeat the trueDiffs multiple times
    tdVect = rep(trueDiffs, times=nSample)
    list(baseDist=baseDist, diffDist=diffDist,trueDiff=tdVect)
  })
  
  ##calculate means for no change
  baseMeans = sapply(dists, function(d){mean(d[[1]])})
  
  ##calculate means for each fold change value
  dd = do.call(c, lapply(dists,"[[",2))
  td = do.call(c, lapply(dists,"[[",3))
  fcMeans = sapply(split(dd,td), function(x){mean(x)})
  
  ##combine with base mean
  fcMeans = c(baseMeans, fcMeans)
  names(fcMeans)[1:length(baseMeans)] = 0
  
  fcMeans
}

############## Helper Functions  ############################

asinh_xform = function(x, a=1, b=1, c=0){
  asinh(a+b*x) + c
}

##uses geneDist in sysdata.rda
rGene = function (n) {
  unifPts = runif(n)
  aVal = approx(geneDist$x, geneDist$y, unifPts, rule = 2)
  pct = 2^aVal$y
  pct/sum(pct)
}

regexToName = function(x) {
  ##collapse bracket stuff
  s = str_match(x,"(.*)\\[(.*)\\](.*)")
  if(is.na(s[1])){return(x)}
  prefix = regexToName(s[2])
  combine = s[3]
  suffix = regexToName(s[4])
  if(nchar(combine)>1){
    for(i in nchar(combine):2){
      combine = paste0(substr(combine,1,i-1),"/",substr(combine,i,nchar(combine)))
    }
  }
  return(paste0(prefix,combine,suffix))
}

merge.table = function(x, y, by=c("row.names", "col.names"), missing=0){
  if(is.null(x)){return(y)}
  if(is.null(y)){return(x)}
  
  by = match.arg(by)
  if(by=="row.names"){
    by.x = row.names(x)
    by.y = row.names(y)
    
    x.miss = by.y[!(by.y %in% by.x)]
    if(length(x.miss)>0){
      x = rbind(x, matrix(missing, length(x.miss), ncol(x), 
                          dimnames=list(x.miss, colnames(x))))
    }
    y.miss = by.x[!(by.x %in% by.y)]
    if(length(y.miss)>0){
      y = rbind(y, matrix(missing, length(y.miss), ncol(y), 
                          dimnames=list(y.miss, colnames(y))))
    }
    return(cbind(x, y[match(rownames(x),rownames(y)),,drop=F]))
    
  }else{
    by.x = colnames(x)
    by.y = colnames(y)
    
    x.miss = by.y[!(by.y %in% by.x)]
    if(length(x.miss)>0){
      x = cbind(x, matrix(missing, ncol(x), length(x.miss),
                          dimnames=list(rownames(x),x.miss)))
    }
    y.miss = by.x[!(by.x %in% by.y)]
    if(length(y.miss)>0){
      y = cbind(y, matrix(missing, ncol(y), length(y.miss), 
                          dimnames=list(rownames(y),y.miss)))
    }
    return(rbind(x, y[,match(colnames(x),colnames(y)),drop=F]))
  }
  
}

strwidth.rot = function(s, cex=NULL){
  if(is.null(cex)){cex = par("cex")}
  xusr <- par("usr") 
  xh <- strwidth(s, cex = par("cex")) 
  yh <- strheight(s, cex = par("cex")) * 5/3 
  
  yh <- xh/(xusr[2]-xusr[1])* par("pin")[1] 
  yh <- yh/ par("pin")[2] * (xusr[4]-xusr[3]) 
  yh
}

strheight.rot = function(s, cex=NULL){
  if(is.null(cex)){cex = par("cex")}
  xusr <- par("usr") 
  xh <- strwidth(s, cex = cex) 
  yh <- strheight(s, cex = cex) * 5/3 
  
  xh <- yh/(xusr[4]-xusr[3])*par("pin")[2] 
  xh <- xh/ par("pin")[1] * (xusr[2]-xusr[1])
  
  xh
}