#' Multiple Permutation test
#'
#' Perform a multiple permutation test between each elements of 2 list of Region
#' set the resut will be store in a multiOvelapPermTest class
#'
#' @usage multiPermTest(Alist, Blist, sampling=FALSE, fraction=0.15, min_sampling=1000,
#'                            ranFun="randomizeRegions",universe=NULL, adj_pv_method="BH", 
#'                            max_pv=0.05,verbose=FALSE,mc.cores=2 ...)
#' 
#' 
#' @param Alist GRangesList or list of Region Set of any accepted formats by  \code{\link{regioneR}} package
#' (\code{\link{GenomicRanges}}, \code{\link{data.frame}} etc...)
#' @param Blist GRangesList or list of Region Set of any accepted formats by  \code{\link{regioneR}} package
#' (\code{\link{GenomicRanges}}, \code{\link{data.frame}} etc...)
#' @param sampling Boolean, if is true the function will use only a sample of
#' each element of Alist to perform the test (default = FALSE)
#' @param fraction Numeric, if sampling==TRUE is the fraction of the region sets
#' used to perform the test (default = 0.15)
#' @param min_sampling numeric, minimum number of regions accepted after the sampling, if the number of the sampled 
#' regions is less than min_sampling will be used min_sampling value as number of regions
#' @param ranFun c("randomizeRegions","circularRandomizeRegions",
#' "resampleRegions"), choose the randomization strategy used for the test see  \code{\link{regioneR}}
#' @param universe (default = NULL) using only when resamplinRegions function is
#' selected
#' @param adj_pv_method Charachter, the method used for the calculation of the adjusted p-value, 
#' to choose between the options of \code{\link{p.adjust}}. (default = "BH")
#' @param  max_pv Numeric, the z-scores associate a p-values higher of this parameter will be transform in 0. (default =0.05)
#' @param  verbose Boolean, if verbose
#' @param ...  further arguments to be passed to other methods. 
#' 
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
#' @export multiOverlapPermTest



multiPermTest<-function(Alist,Blist=NULL, 
                               sampling=FALSE,fraction=0.15, min_sampling=1000,
                               ranFun=randomizeRegions, evalFun=numOverlaps, universe=NULL,
                               adj_pv_method="BH", max_pv=0.05, pv="adj.pv", 
                               verbose=FALSE,subEx=NA,mc.cores=2,...){
  
  paramList<-list(Alist=deparse(substitute(Alist)),
                  Blist=deparse(substitute(Blist)),
                  sampling=deparse(substitute(sampling)),
                  fraction=deparse(substitute(fraction)),
                  min_sampling=deparse(substitute(fraction)),
                  ranFun=deparse(substitute(ranFun)),
                  universe=deparse(substitute(universe)),
                  adj_pv_method=adj_pv_method,
                  max_pv=deparse(substitute(max_pv)),
                  mc.cores=mc.cores)
  
  if(is.null(Blist)){Blist<-Alist}
  list.tabs<-list()
  list.pt<-list()
  if (sampling==TRUE){  
    Alist<-subList(Alist,min_sampling=min_sampling,fraction=fraction)  
  }
  
  if(paramList$ranFun=="resampleRegions" & is.null(universe)){
    if (is.null(universe)){
      warning("resampleRegions function need that universe parameters in not NULL universe will created using all the regions present in Blist")
      universe<-createUniverse(Alist) # check well this option
    }
  }
  
  if(is.null(names(Alist))){names(Alist)<-paste0("Alist-",1:length(Alist))}
  if(is.null(names(Blist))){names(Blist)<-paste0("Blist-",1:length(Blist))}
  
  for ( i in 1:length(Alist)){
    print(names(Alist[i]))
    A<-Alist[[i]]
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=evalFun, param.name="B", values=Blist)
    ptm <- proc.time()
    
    pt <- permTest( A=A, evaluate.function=func.list, 
                    randomize.function=randomizeRegions,count.once=TRUE,mc.cores=mc.cores,...)
    
  
  
  time<-proc.time() - ptm
  time<-time[3]/60
  if (verbose==TRUE){print(paste0(" run in ",time,"  minute"))}
  tab<-data.frame()
  
  for (j in 1:length(pt)){
    if(pt[[j]]$zscore==0 | is.na(pt[[j]]$zscore) | is.nan((pt[[j]]$zscore))){ # modifica per non andare in errore
      zscore.norm<-0
      zscore.std<-0
    }else{
      zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
    }
    
    vec<-data.frame(order.id=j,
                    name=new.names[j],
                    n_regions=length(Blist[[j]]),
                    z_score=pt[[j]]$zscore,
                    p_value=pt[[j]]$pval,
                    n_overlaps=pt[[j]]$observed,
                    mean_perm_test=mean(pt[[j]]$permuted),
                    sd__perm_test=sd(pt[[j]]$permuted))
    tab<-rbind(vec,tab)
  }
  tab$norm_zscore<-tab$z_score/sqrt(length(A))
  tab$ranged_zscore<-rangedVector(tab$norm_zscore)
  tab$adj.p_value<-round(p.adjust(tab$p_value,method=adj_pv_method),digits = 4)
  tab<-funRemove(tab,max_pv=max_pv,pv=pv,subEx=subEx)
  if (verbose==TRUE){print(tab)}  # remember to activate only if verbose.... 
  list.tabs[[i]]<-tab
  names(list.tabs)[i]<-names(Alist)[i]
}
list.tabs<-c(list(parameters=as.list(paramList)),list.tabs)
class(list.tabs) <- "multiPermTest"


return(list.tabs)
}



