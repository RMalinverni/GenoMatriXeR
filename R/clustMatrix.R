#' Clusterize Matrix
#'
#' apply a clusterization to a numeric matrix
#'
#' @usage clustMatrix(mat,method="complete",mirrored=FALSE,...)
#'
#' @param mat matrix
#' @param method string method for the (\code{\link{hclust}}
#' @param mirrored boolean choose if the matrix will be mirrored
#'
#' @seealso ...
#'
#' @examples  ...
#'
#' #@references ...
#' #@importFrom ...
#' @export clustMatrix
#'
clustMatrix<-function(mat,method="complete",mirrored=FALSE,lZ_tab=FALSE,...){
  if(lZ_tab==TRUE){
    d<-dist(mat)
    hc<-hclust(d)
    ind<-hc$order
    mat1<-mat[ind,]
  }else{

    if(mirrored==FALSE){
      d<-dist(mat)
      hc<-hclust(d,...)
      ind<-hc$order
      d1<-dist(t(mat))
      hc1<-hclust(d1,method = method)
      ind1<-hc1$order
      mat1<-mat[ind,ind1]
    }
    if(mirrored==TRUE){
      ordx<-order(colnames(mat))
      ordy<-order(rownames(mat))
      mat<-mat[ordy,ordx]
      d<-dist(mat)
      hc<-hclust(d,method = method)
      ind<-hc$order

      mat1<-mat[ind,ind]
    }
  }
  return(mat1)
}
