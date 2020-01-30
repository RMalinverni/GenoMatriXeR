#' Plot local Z-Score Matrix
#'
#' Create a image local Z-Score Matrix
#'
#' @usage plotMatrix( mat, ...)
#'
#' @param mat association matrix
#'
#' @param ...  further arguments to be passed to other methods.
#'
#'
#' @seealso \code{\link{plotMatrix}}
#'
#' @examples  ...
#'
#' @export plotLZMatrix
#'
plotLZMatrix<-function(mat,...){
  plotMatrix(t(mat),type="compare",palette="gbr")
}
