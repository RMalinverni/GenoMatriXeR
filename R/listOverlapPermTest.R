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

listOverlapPermTest<-function(Alist,Blist,sampling=FALSE,fraction=0.15,
                              ranFun="randomizeRegions",universe=NULL,mc.cores=2,...){
  list.tabs<-list()
  list.pt<-list()
  for ( i in 1:length(Alist)){
    print(names(Alist[i]))
    A<-Alist[[i]]
    if(sampling==TRUE){
      A<-A[sample(length(A),round(length(A)*fraction))] #controllare se è vero
    }
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=numOverlaps, param.name="B", values=Blist)

    ptm <- proc.time()
    # pt <- tryCatch({permTest(A=A, evaluate.function=func.list,  #aggiungere TryCatch
    #                randomize.function=randomizeRegions)}, error=function(e){
    #                  print ("error detected")
    #                  as.list(rep(NA,length(Blist)))})

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
    print(paste0(" run in ",time,"  minute"))
    tab<-data.frame()
    for (j in 1:length(pt)){
      if(pt[[j]]$zscore==0){
        zscore.norm<-0
        max.val<-0
        max.z<-0
        max.z.norm<-0
        zscore.std<-0
      }else{

        if(pt[[j]]$zscore>0){ # control NA
          zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
          max.val<-min(length(A),length(Blist[[j]]))
          max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
          max.z.norm<-max.z/sqrt(length(A))
          zscore.std<-pt[[j]]$zscore/max.z
        }
        if(pt[[j]]$zscore<=0){
          zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
          max.val<-0
          max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
          max.z.norm<-max.z/sqrt(length(A))
          zscore.std<-abs(pt[[j]]$zscore)/max.z
        }
      }
      vec<-data.frame(order.id=j,
                      name=new.names[j],
                      n_regions=length(Blist[[j]]),
                      z_score=pt[[j]]$zscore,
                      norm_zscore=zscore.norm,
                      std_zscore=zscore.std,
                      p_value=pt[[j]]$pval,
                      n_overlaps=pt[[j]]$observed,
                      mean_perm_test=mean(pt[[j]]$permuted),
                      sd__perm_test=sd(pt[[j]]$permuted))
      tab<-rbind(vec,tab)
    }
    #colnames(tab)<-c("order.id","name","nº region","z-score","norm.z-score","standard.z-score","p-value","n.overlaps","mean.perm","sd.perm")
    print(tab)
    list.tabs[[i]]<-tab
    #list.pt[[i]]<-pt
    names(list.tabs)[i]<-names(Alist)[i]
    #names(list.pt)[i]<-names(Alist)[i]

  }
  return(list.tabs)
}
