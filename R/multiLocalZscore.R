#' Multiple Local Z-Score test
#'
#' Perform a multiple permutation test and local z-score calculation between a region set and 
#' list of regions set (or GrangedList)
#' 
#'
#' @usage listOvelapPermtest( Alist, Blist, sampling=FALSE, fraction=0.15,
#' ranFun="randomizeRegions", ...)
#'
#' @param A Region Set of any accepted formats by  \code{\link{regioneR}} package
#' (\code{\link{GenomicRanges}}, \code{\link{data.frame}} etc...)
#' @param Blist list of Region Set of any accepted formats by \code{\link{regioneR}} package
#' (\code{\link{GenomicRanges}}, \code{\link{data.frame}} etc...)
#' @param sampling Boolean, if is true the function will use only a sample of
#' each element of Alist to perform the test
#' @param fraction Numeric, if sampling==TRUE is the fraction of the region sets
#' used to perform the test
#' @param ranFun Function, choose the randomization strategy used for the test (default)
#' for details see \code{\link{regioneR}} 
#' @param universe Region Set of any accepted formats by  \code{\link{regioneR}}, using only when resamplinRegions function is
#' selected (default = NULL) 
#' @param ...  further arguments to be passed to other methods.
#' @param adj_pv_method Charachter, the method used for the calculation of the adjusted p-value, 
#' to choose between the options of \code{\link{p.adjust}}. (default = "BH")
#' @param  max_pv Numeric, the z-scores associate a p-values higher of this parameter will be transform in 0. (default =0.05)
#' @param  verbose Boolean, if verbose
#' @details  the permutation test core used in this function permit to change
#' \code{"randomize.function"} \code{\link{randomizeRegions}},
#' \code{\link{circularRandomizeRegions}}, \code{\link{resampleRegions}} or a
#' custom function), but use only an \code{"evaluation.function"}
#' \code{\link{numOverlaps}}
#' @return ...
#'
#' @seealso  \code{\link{permTest}}
#'
#' @examples  ...
#'
#' @export listOverlapPermTest



multiLocalZscore <- function(A, Blist=NULL,ranFun=randomizeRegions,sampling=FALSE,
                            fraction=0.15, universe=NULL, window=2000, step=100,
                            mc.cores=2,adj_pv_method="BH",max_pv=0.05,genome="hg19",...){
  
  paramList<-list(A=deparse(substitute(A)),
                  Blist=deparse(substitute(Blist)),
                  sampling=deparse(substitute(sampling)),
                  fraction=deparse(substitute(fraction)),
                  #min_sampling=deparse(substitute(fraction)),
                  ranFun=deparse(substitute(ranFun)),
                  universe=deparse(substitute(universe)),
                  window=window,
                  step=step,
                  adj_pv_method=adj_pv_method,
                  max_pv=deparse(substitute(max_pv)),
                  mc.cores=mc.cores)
  
  A<-toGRanges(A)
  
  if(sampling==TRUE){
    A<-A[sample(length(A),round(length(A)*fraction))] 
  }
  
  if(paramList$ranFun=="resampleRegions" & is.null(universe)){
    if (is.null(universe)){
      print("resampleRegions function need that universe parameters in not NULL universe will created using all the regions present in Blist")
      universe<-createUniverse(Blist) # check well this option
    }
  }
  
  function.list <- createFunctionsList(FUN = numOverlaps, param.name = "B", values = Blist)
  
  pt <- permTest(A = granges(A),
                 evaluate.function = function.list,
                 randomize.function = ranFun,genome=genome,
                 count.once = TRUE,
                 universe=universe)
  
  
  lZs <- list()
  for(i in 1:length(pt)){
    lZs[[i]] <- localZScore(A = A,  pt = pt[[i]], count.once=TRUE,window = window, step = step)
    names(lZs)[i]<-names(pt[i])
  }
  
  
  Nreg<-length(A)
  p_values<-do.call(c,pt)
  pval<-p_values[grep(".pval",names(p_values))]
  means_pemuted<-lapply(p_values[grep(".permuted",names(p_values))],mean)
  sd_pemuted<-lapply(p_values[grep(".permuted",names(p_values))],sd)
  z_score<-p_values[grep(".zscore",names(p_values))]
  observed<-p_values[grep(".observed",names(p_values))]
  localZs<-do.call(c,lZs)
  shiftedZs<-localZs[grep("shifted.z.scores",names(localZs))]
  shifts<-localZs[grep(".shifts",names(localZs))][[1]]
  
  
  if(is.null(names(Blist))){names(Blist)<-1:length(Blist)}
  
  tab<-data.frame(name = names(Blist),
                  p_value = unlist(pval),
                  z_score =unlist(z_score),
                  mean_perm_test=unlist(means_pemuted),
                  sd_perm_test=unlist(sd_pemuted), 
                  n_overlaps=unlist(observed))
  tab$norm_zscore<-tab$z_score/sqrt(Nreg)
  maxzscores<-(Nreg-tab$mean_perm_test)/tab$sd_perm_test
  tab$ranged_zscore<-tab$z_score/maxzscores
  tab$adj.p_value<-round(p.adjust(tab$p_value,method=adj_pv_method),digits = 4)
  
  names(shiftedZs)<-names(Blist)
  
  multiLZ.obj<-list(param=paramList,
                  nreg=Nreg,
                  resumeTab=tab,
                  max_zscores=maxzscores,
                  shifts=shifts,
                  shifed_ZSs=shiftedZs
  )
  
  class(multiLZ.obj)<-"multiLocalZscore"
  return (multiLZ.obj)
}




