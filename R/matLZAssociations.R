#' Matrix from lZAssociation Object
#'
#' transform the lZAssociation Object in a numeric matrix
#'
#' @usage matLZAssociations ( lZobj, A = NULL, ...)
#'
#' @param lZobj object obtained from lZAssociation function
#' @param A (default = NULL) region set
#'
#' @seealso \code{\link{lZAssociation}}
#'
#' @examples  ...
#'
#' #@references ...
#' #@importFrom ...
#' @export matLZAssociations
#'
matLZAssociations <- function(lZobj,A=NULL,...){
  lZs<-lZobj$local_zscores
  tab<-vector()
  coln<-seq(from = 0-lZs[[1]]$window, to=0+lZs[[1]]$window,by=lZs[[1]]$step)
  for(i in 1:length(lZs)){
    tab<-rbind(tab,lZs[[i]][[1]])
  }
  rownames(tab)<-names(lZs)
  colnames(tab)<-coln
  mat<-as.matrix(tab)
  if (!is.null(A)){
    mat<-mat/sqrt(length(toGRanges(A)))
  }
  return(mat)
}
