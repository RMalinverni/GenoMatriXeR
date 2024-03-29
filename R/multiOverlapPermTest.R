#' Multiple Permutation test
#'
#' Perform a multiple permutation test between each elements of 2 list of Region
#' set the resut will be store in a multiOvelapPermTest class
#'
#' @usage multiOvelapPermTest(Alist, Blist, sampling=FALSE, fraction=0.15, min_sampling=1000,
#'                            ranFUN="randomizeRegions",universe=NULL, adj_pv_method="BH", 
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
#' @param ranFUN (default = "randomizeRegions") choose the randomization strategy used for the test see  \code{\link{regioneR}}
#' @param evFUN  (default = "numOverlaps) choose the evaluation strategy used for the test see  \code{\link{regioneR}}
#' @param universe (default = NULL) used only when \code{\link{resampleRegions}} function is selected
#' @param adj_pv_method Charachter, the method used for the calculation of the adjusted p-value, 
#' to choose between the options of \code{\link{p.adjust}}. (default = "BH")
#' @param max_pv Numeric, the z-scores associate a p-values higher of this parameter will be transform in 0. (default =0.05)
#' @param pv Charachter, (default = "adj.pv") value choose to set the cutoff
#' @param subEx Numeric, (default = 0) substitute this value to a z-score when the pvalue is higher than max_pv
#' @param genome Charachter, (defalut = "hg19") genome used to compute the randomization
#' @param verbose Boolean, if verbose
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



multiOverlapPermTest<-function(Alist, Blist = NULL, sampling = FALSE, fraction = 0.15, min_sampling = 5000,
                              ranFUN = "randomizeRegions", evFUN = "numOverlaps", universe = NULL, ntimes=100,
                              adj_pv_method = "BH", max_pv = 0.05, pv = "adj.pv", subEx=0, 
                              genome = "hg19", verbose = FALSE, ...){
  
  ranFUN<-as.character(substitute(ranFUN)) 
  
  paramList<-list(Alist=deparse(substitute(Alist)),
              Blist=deparse(substitute(Blist)),
              sampling=deparse(substitute(sampling)),
              fraction=deparse(substitute(fraction)),
              min_sampling=deparse(substitute(fraction)),
              ranFUN=ranFUN,
              universe=deparse(substitute(universe)),
              adj_pv_method=adj_pv_method,
              max_pv=deparse(substitute(max_pv)),
              ntimes=ntimes,
              nc=NULL,
              matOrder=NULL
              )
 
  if( is.null ( Blist ) ){ Blist <- Alist }
  list.tabs <- list()
  list.pt <- list()
  if ( sampling == TRUE ){  
    Alist <- subList (Alist, min_sampling = min_sampling, fraction = fraction)  
  }
  
  
  
  
  if (ranFUN=="resampleRegions") {FUN<-resampleRegions}
  if (ranFUN=="randomizeRegions") {FUN<-randomizeRegions}
  if (ranFUN=="circularRandomizeRegions") {FUN<-circularRandomizeRegions}
  
  print(ranFUN)
  
  for ( i in 1:length(Alist)){
    
    print(names(Alist[i]))
    A<-Alist[[i]]
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=evFUN, param.name="B", values=Blist)
    ptm <- proc.time()
    
    if(ranFUN=="resampleRegions"){
      if (is.null(universe)){
        warning("resampleRegions function need that 'universe' is not NULL, universe was created using all the regions present in Alist")
        uniList<-data.frame()
        for(u in 1:length(Alist)){
          df<-toDataframe(Alist[[u]])[,1:3]
          uniList<-rbind(uniList,df)
        }
        universe<-uniList
      }
    }
    
    pt <- permTest( A=A, evaluate.function=func.list, ntimes=ntimes, #aggiungere TryCatch
                    randomize.function=FUN,universe=universe,count.once=TRUE, genome=genome, ...)
    
  
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
                      n_regionA=length(Alist[[i]]),
                      n_regionB=length(Blist[[j]]),
                      z_score=pt[[j]]$zscore,
                      p_value=pt[[j]]$pval,
                      n_overlaps=pt[[j]]$observed,
                      mean_perm_test=mean(pt[[j]]$permuted),
                      sd_perm_test=sd(pt[[j]]$permuted))
      tab<-rbind(vec,tab)
    }
    tab$norm_zscore<-tab$z_score/sqrt(tab$n_regionA)
    max_zscore<-(tab$n_regionA-tab$mean_perm_test)/tab$sd_perm_test
    tab$std_zscore<-tab$z_score/max_zscore
    tab$adj.p_value<-round(p.adjust(tab$p_value,method=adj_pv_method),digits = 4)
    
    
    if (verbose==TRUE){print(tab)}  # remember to activate only if verbose.... 
    list.tabs[[i]]<-tab
    names(list.tabs)[i]<-names(Alist)[i]
  }
  GMXRobj<-gMXR(parameters=paramList,multiOverlaps=list.tabs ,matrix = list(NULL))

  
  return(GMXRobj)
}
