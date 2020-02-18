#' Multiple Permutation test
#'
#' Perform a multiple permutation test between each element of 2 list of Region
#' set
#'
#' @usage listOvelapPermtest( Alist, Blist, sampling=FALSE, fraction=0.15,
#' ranFun="randomizeRegions", ...)
#'
#' @param Alist list of Region Set of any accepted formats by  \code{\link{regioneR}} package
#' (\code{\link{GenomicRanges}}, \code{\link{data.frame}} etc...)
#' @param Blist list of Region Set of any accepted formats by \code{\link{regioneR}} package
#' (\code{\link{GenomicRanges}}, \code{\link{data.frame}} etc...)
#' @param sampling boolean if is true the function will use only a sample of
#' each element of Alist to perform the test
#' @param fraction numeric, if sampling==TRUE is the fraction of the region sets
#' used to perform the test
#' @param ranFun c("randomizeRegions","circularRandomizeRegions",
#' "resampleRegions") choose the randomization strategy used for the test
#' @param universe (default = NULL) using only when resamplinRegions function is
#' selected
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
#' @export listOverlapPermTest



listOverlapPermTest2<-function(Alist,Blist,sampling=FALSE,fraction=0.15,min_sampling=1000,
                              ranFun="randomizeRegions",universe=NULL,
                              verbose=FALSE,mc.cores=2,...){
  list.tabs<-list()
  list.pt<-list()
  if (sampling==TRUE){   # I am making only the subsampling of Alist check if it is clever or stupid
    Alist<-subList(Alist,min_sampling=min_sampling,fraction=fraction)  
  }
  
  for ( i in 1:length(Alist)){
    print(names(Alist[i]))
    A<-Alist[[i]]
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=numOverlaps, param.name="B", values=Blist)
    ptm <- proc.time()
    
    if(ranFun=="randomizeRegions"){
      pt <- permTest( A=A, evaluate.function=func.list,  #aggiungere TryCatch
                      randomize.function=randomizeRegions,count.once=TRUE,mc.cores=mc.cores)
    }
    if(ranFun=="circularRandomizeRegions"){
      pt <- permTest(A=A,evaluate.function=func.list,
                     randomize.function=circularRandomizeRegions,count.once=TRUE,mc.cores=mc.cores)
    }
    if(ranFun=="resampleRegions"){
      if (is.null(universe)){
        print("resampleRegions function need that universe parameters in not NULL universe will created using all the regions present in Alist")
        uniList<-data.frame()
        for(u in 1:length(Alist)){
          df<-toDataframe(Alist[[u]])[,1:3]
          uniList<-rbind(uniList,df)
        }
        universe<-uniList
      }
      pt <- permTest(A=A,evaluate.function=func.list,
                     randomize.function=resampleRegions,universe=universe,
                     count.once=TRUE,mc.cores=mc.cores)

    }

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
    if (verbose==TRUE){print(tab)}  # remenber to activate only if verbose.... 
    list.tabs[[i]]<-tab
    names(list.tabs)[i]<-names(Alist)[i]
  }
  
  funRemove<-function(x,max_pv=0.05){     # funzione per azzerare gli zScore che non passano in test 
    x$z_score[x$p_value>max_pv]<-0        # da portare fuori
    x$norm_zscore[x$p_value>max_pv]<-0
    x$ranged_zscore[x$p_value>max_pv]<-0
    return(x)
  }
  lapply(list.tabs,funRemove)
  return(list.tabs)
}



