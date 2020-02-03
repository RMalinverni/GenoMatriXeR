#' Multiple Local-Z score test
#'
#' calculate a a local-z score test btween a Region Set and a listof Region Set
#'
#' @usage lZAssociations(A, Blist,ranFun="randomizeRegions",universe=NULL,
#' window=2000,step=100,...)
#'
#' @param A test Region Set
#'
#' @param Blist test Region Set List
#' @param ranFun c("randomizeRegions","circularRandomizeRegions",
#' "resampleRegions") choose the randomization strategy used for the test
#' @param universe (default = NULL) using only when resamplinRegions function is
#' selected
#' @param sampling boolean if is true the function will use only a sample of
#' each element of Alist to perform the test
#' @param fraction numeric, if sampling==TRUE is the fraction of the region sets
#' used to perform the test
#' @param window (default 2000) window associated to \code{\link{localZScore}}
#' function
#' @param step (default 100) step associated to \code{\link{localZScore}}
#' @param ...  further arguments to be passed to other methods.
#'
#'
#' @seealso \code{\link{localZScore}} \code{\link{regioneR}}
#'
#' @examples  ...
#'
#' @export lZAssociations
#'
lZAssociations <- function(A, Blist,ranFun="randomizeRegions",sampling=FALSE,
                           fraction=0.15, universe=NULL,window=2000,step=100,...){
  
  if(sampling==TRUE){
    A<-A[sample(length(A),round(length(A)*fraction))] #controllare se è vero
  }
  
  function.list <- createFunctionsList(FUN = numOverlaps, param.name = "B", values = Blist)

  if(ranFun=="randomizeRegions"){
    pts <- permTest(A = A,
                    evaluate.function = function.list,
                    randomize.function = randomizeRegions,count.once = TRUE)
  }
  if(ranFun=="circularRandomizeRegions"){
    pts <- permTest(A = A,
                    evaluate.function = function.list,
                    randomize.function = circularRandomizeRegions,
                    count.once = TRUE)
  }
  if(ranFun=="resampleRegions"){
    if (is.null(universe)){
      print("resampleRegions function need that universe parameters in not NULL universe will created using all the regions present in Blist")
      uniList<-data.frame()
      for(u in 1:length(Blist)){
        df<-toDataframe(Blist[[u]])[,1:3]
        uniList<-rbind(uniList,df)
      }
      universe<-uniList
    }
    pts <- permTest(A=A,evaluate.function=function.list,
                   randomize.function=resampleRegions,universe=universe,
                   count.once=TRUE)

  }


  lZs <- list()
  for(i in 1:length(pts)){
    lZs[[i]] <- localZScore(A = A,  pt = pts[[i]], window = window, step = step)
    names(lZs)[i]<-names(pts[i])
  }
  return (list(p_values=pts,local_zscores=lZs))
}
