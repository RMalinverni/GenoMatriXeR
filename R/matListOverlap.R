#' Matrix from ListOverlap Object
#'
#' transform the ListOverlap Oject in a numeric matrix
#'
#' @usage matListOverlap(listOverlapPermTest.obj,zs.type="ranged_zscore",...)
#'
#' @param listOverlapPermTest.obj object from function
#' (\code{\link{listOverlapPermTest}}
#' @param zs.type c("std","norm") choose if create the matrix using standard or
#' normalize Z-score#' @param nameA string name for the plot related to matrix A
#'
#' @seealso ...
#'
#' @examples  ...
#'
#' #@references ...
#' #@importFrom ...
#' @export matListOverlap
#'
matListOverlap<-function(listOverlapPermTest.obj,zs.type='ranged_zscore',...){
  A.obj<-listOverlapPermTest.obj[-1]
  mat<-vector()
  for (i in 1:length(A.obj)){
      mat<-cbind(mat,as.numeric(A.obj[[i]][,'ranged_zscore']))
  }
  colnames(mat)<-names(A.obj)
  rownames(mat)<-A.obj[[1]][,2]
  mat<-as.matrix(mat)
  return(mat)
}
