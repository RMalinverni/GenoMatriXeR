makeGenomicMatrix<-function(mPT,clusterize=TRUE,hc.method="ward.D",dist.method="euclidean", transform=FALSE,
         nboot=1000,zs.type='norm_zscore', symm_matrix=TRUE,nc=NULL,...){
  
  if (class(mPT)=="GenoMatriXeR"){     #add an error in other case
    mat<-matMultiPermTest(mPT,zs.type=zs.type)
    print("class recognize")
  }
  
  if(transform==TRUE){
    mat<-t(mat)
  }
  
  if (symm_matrix==TRUE & (ncol(mat)!=nrow(mat))){
    symm_matrix<-FALSE
    warning("impossible to create symmetrical matrix, number of matrix columns is different from number of rows")
  }
  
  
  if (clusterize==TRUE){
    fit <- pvclust(mat, method.hclust=hc.method, method.dist=dist.method,
                   nboot = nboot)
    ind<-fit$hclust$order
    if (symm_matrix==TRUE){
      
      newNames<-colnames(mat)[ind] # non mi piace ordinare per nome
      mat<-mat[newNames,newNames]
      fit2 <- NULL
      
    }else{
      
      fit2 <- pvclust(t(mat), method.hclust=hc.method, method.dist=dist.method,
                      nboot = nboot)
      ind2<-fit2$hclust$order
      mat<-mat[ind2,ind]
      
    }  
    
    if (is.null(nc)){
      set.seed(123)
      nc<-pamk(t(mat),krange = 1:(ncol(mat)-1))$nc
    }
  }else{
    fit=NULL
    fit2=NULL
  }
  
  
  clus <- kmeans(t(mat), centers=nc)
  
  mat1<-list(GMat=mat,GFit=fit,GFit2=fit2,GKmean=clus)
  mPT@matrix<-mat1
  #class(mat1)<-"GenomicMatrix"
  
  return(mPT)
  
}