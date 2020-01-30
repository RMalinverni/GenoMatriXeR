matGlobalAssociations <- function(obj,A=NULL,...){
  lZs<-obj$local_zscores
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
