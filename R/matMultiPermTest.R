#' Matrix from multiPermTest Object
#'
#' transform the ListOverlap Oject in a numeric matrix
#'
#' @usage matMultiPermTest(multiPermTest.obj,zs.type='ranged_zscore',...)
#'
#' @param multiPermTest.obj,  object from function (\code{\link{multiPermTest}}
#' 
#' @param zs.type Character, choose if create the matrix using 'ranged_zscore' or
#' normal_zscore' 
#' 
#' @param nameA string name for the plot related to matrix A
#'
#' @seealso ...
#'
#' @examples  ...
#'
#' #@references ...
#' #@importFrom ...
#' @export matMultiPermTest
#'
matMultiPermTest<-function(mPT,zs.type='ranged_zscore',...){
  A.obj<-mPT@multiOverlaps
  mat<-vector()
  for (i in 1:length(A.obj)){
      mat<-cbind(mat,as.numeric(A.obj[[i]][,zs.type]))
  }
  colnames(mat)<-names(A.obj)
  rownames(mat)<-A.obj[[1]][,2]
  mat<-as.matrix(mat)
  return(mat)
}
