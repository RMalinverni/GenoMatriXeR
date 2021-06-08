makeLZMatrix<-function(mlZA,normalize=TRUE,matLim=NULL){
  
  if (class(mlZA)!="multiLocalZscore"){     
    stop("the object mlZA need to be an multiLocalZscore object")
  }
  
  vec<-vector(length = length(mlZA$shifts))
  if (normalize==TRUE){
    for(i in 1:length(mlZA$shifed_ZSs)){
      vec<-rbind(vec,mlZA$shifed_ZSs[[i]]/sqrt(mlZA$max_zscores[i]))
    }
  }else{
    for(i in 1:length(mlZA$shifed_ZSs)){
      vec<-rbind(vec,mlZA$shifed_ZSs[[i]])
    } 
  }

  vec<-vec[-1,]
  
  if (is.vector(vec)){
    vec<-t(as.data.frame(vec))
  }
  rownames(vec)<-names(mlZA$shifed_ZSs)
  colnames(vec)<-mlZA$shifts
  # I need to add a matLim integration
  
  return(vec)
}



