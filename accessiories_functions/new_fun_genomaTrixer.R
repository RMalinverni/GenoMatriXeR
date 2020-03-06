
# choose systematic names for functions
collapse_list<-function(list_GR){  
  GR<-GRanges()
  for (i in 1:length(list_GR)){
    GR<-c(GR,filterChromosomes(list_GR[[i]]))  
  }
  return(joinRegions(GR))
}

listRGBinTable<-function(listGRs,bin=1000,column=NULL,genGR=NULL,genome="hg19",collapse=T){
  
  if (is.null(genGR)){
    gen<-getGenome(genome)
  }else{
    gen<-genGR
  }
  
  tiled_Gen<-unlist(tile(gen,width=bin))
  
  if (collapse==T){
    validRanges<-collapse_list(listGRs)  
    ind<-overlapRegions(tiled_Gen,validRanges,only.boolean = T)
    tiled_Gen<-tiled_Gen[ind]
  }
  
  mta<-vector()
  if (is.null(column)){
    for(i in 1:length(listGRs)){
      newcol<-as.numeric(overlapRegions(tiled_Gen,listGRs[[i]],only.boolean=T))
      mta<-cbind(mta,newcol)
    }
  }else{
    for(i in 1:length(listGRs)){
      newcol<-overlapRegions(tiled_Gen,listGRs[[i]],colB=column)$dd
      newcol<-newcol/mean(newcol)
      mta<-cbind(mta,newcol)
    }
  }
  
  colnames(mta)<-names(listGRs)
  mcols(tiled_Gen)<-mta
  
  return(tiled_Gen)
}


randomizeRegionsPerc<-function(GR,genome,frac=0.2,...){
  nc<-round(length(GR)*frac)
  change<-sample(length(GR),nc)
  GR1<-GR[-change]
  GR2<-randomizeRegions(GR[change],genome=genome)
  GR3<-c(GR1,GR2)
  return(GR3)
}

testAssociations<-function(A,B,
                           seqFraction=seq(0.1,1,0.1),
                           ntimes=100){

  vecZs<-vecNZs<-vecSZs<-vecOE<-vecP<-vecPV<-vector()
  for(i in seq(0.1,1,0.1)){
    ind<-sample(length(A),round(length(A)*i))
    print(length(ind))
    A1=A[ind]
    pt<-permTest(A=A1,
                 B=B,
                 randomize.function = randomizeRegions,
                 evaluate.function = numOverlaps,
                 ntimes=200,
                 count.once=TRUE)
    orig.ev<-pt$numOverlaps$observed
    rand.ev<-pt$numOverlaps$permuted
    zscore <- round((orig.ev - mean(rand.ev, na.rm = TRUE))/stats::sd(rand.ev,na.rm = TRUE), 4)
    maxZscore <- round((length(A1) - mean(rand.ev, na.rm = TRUE))/stats::sd(rand.ev,na.rm = TRUE), 4)
    normZS<-zscore/sqrt(length(A1))
    normZSmax<-maxZscore/sqrt(length(A1))
    stZscore<-normZS/normZSmax
    vecZs<-c(vecZs,zscore)
    vecNZs<-c(vecNZs,normZS)
    vecSZs<-c(vecSZs,stZscore) #cambiare
    vecOE<-c(vecOE,orig.ev)
    vecPV<-c(vecPV,pt$numOverlaps$pval)
    vecP<-c(vecP,mean(rand.ev, na.rm = TRUE))
  }
  
  DF<-DataFrame(frac=seq(0.1,1,0.1),
                Nreg=round(length(A)*seq(0.1,1,0.1)),
                zcore=round(vecZs,digits = 2),
                p_value=vecPV,
                norm_zscore=round(vecNZs,digits = 2),
                stZscore=round(vecSZs,digits = 2),
                n_ov=vecOE,
                menPerm_ov=round(vecP,digits = 2))
  return(DF)
}



cleanNames<-function(listRS,cellName,
                     vecEx=c("Uni","V0","Pcr","sc[0-9]","anb[0-9]","Iggrab","Forskln","Ucd")){
  newNames<-gsub(paste0("^(.*)",cellName),"",names(listRS))
  for(i in vecEx){
    newNames<-gsub(paste0(i,"(.*)$"),"",newNames) 
  }
  names(listRS)<-newNames
  listRS<-listRS[!is.na(names(listRS))]
  names(listRS)<-paste0(1:length(names(listRS)),"-",names(listRS))
  return(listRS)
}



rangedVector <- function(x){   # cambiare nome
  ind<-is.finite(x)
  x[!is.finite(x)]<-max(x[ind])*2 # this is to change- I made only to have a value
  if(sum(x[ind])!=0){
    tr_vec<-as.numeric(x>0)
    tr_vec[tr_vec==0]<-(-1)
    x<-abs(x)
    ranged_x<-(x-0)/(max(x)-0)
    ranged_x<-ranged_x*tr_vec
  }else{ranged_x<-rep(0,length(x))}
  return(ranged_x)
}

subList<-function(Alist,min_sampling,fraction){
  subAlist<-list()
  for (i in 1:length(Alist)){
    A<-Alist[[i]]
    if(min_sampling<length(A)){   
      subN<-round(length(A)*fraction)
      if (subN<min_sampling) {subN<-min_sampling}
      A<-A[sample(length(A),subN)]
      subAlist[[i]]<-A
    }else{subAlist[[i]]<-A}
  }
  names(subAlist)<-names(Alist)
  return(subAlist)
}

funRemove<-function(x,max_pv=0.05,pv="adj.pv",subEx=0,...){     # funzione per azzerare gli zScore che non passano in test 
  
  
  
  if(pv=="adj.pv"){
    x$z_score[x$adj.p_value>max_pv]<-subEx     
    x$norm_zscore[x$adj.p_value>max_pv]<-subEx
    x$ranged_zscore[x$adj.p_value>max_pv]<-subEx
  }
  
  if(pv=="p.value"){
    x$z_score[x$p_value>max_pv]<-subEx     
    x$norm_zscore[x$p_value>max_pv]<-subEx
    x$ranged_zscore[x$p_value>max_pv]<-subEx
  }
  
  return(x)
}

reSizeRegions<-function(regGR,method,width=NULL){
  
  if(method=="tile"){
    
    if (!is.null(width)){
      regGR<-unlist(tile(regGR,width = width))
    }else{
      warning("if you do not provide a value for 'width' will be use as 'width'  the mean of the width of the provided GRanges")
      regGR<-unlist(tile(regGR,width = mean(width(regGR))))
    }
  }
  
  if(method=="homogenize"){
    meanWidth<-mean(width(regGR))
    sdWidth<-sd(width(regGR))
    acceptedReg<-regGR[meanWidth<=meanWidth+(sdWidth*2)]
    removedReg<-regGR[meanWidth>meanWidth+(sdWidth*2)]
    removedReg<-unlist(tile(removedReg,width = meanWidth))
    regGR<-c(acceptedReg,removedReg)
  }
  
  return(regGR)
}

clustGenomicMatrix<-function(mat,hc.method="ward.D",dist.method="euclidean",nboot=1000,...){
  fit <- pvclust(mat, method.hclust="ward.D2",
                 method.dist="euclidean",nboot = nboot,...)
  ind<-fit$hclust$order
  newNames<-colnames(mat)[ind]
  mat<-mat[newNames,newNames]
  
  mat1<-list(GMat=mat,GFit=fit)
  class(mat1)<-"GenomicMatrix"
  return(mat1)
}


createRandomRegionSets<-function(RegionSet,frac=0.5,seedName='',seed=123,nregions=10,genome="hg19",...){
  set.seed(seed)
  lista<-list()
  for(i in 1:nregions){
    if (i==1){
      lista[[i]]<-randomizeRegionsPerc(RegionSet,frac=frac,genome=genome)
    }else{
      lista[[i]]<-randomizeRegionsPerc(lista[[i-1]],frac=frac,genome=genome)
    }
    names(lista)[[i]]<-paste0(seedName,i)
  }
  return(lista)
}


makeGenomicMatrix<-function(mPT,hc.method="ward.D",dist.method="euclidean",
                            nboot=1000,zs.type='ranged_zscore', symm_matrix=TRUE,...){
  
  
  # if (symm_matrix==TRUE & (ncol(mPT)!=nrow(mPT))){
  #   symm_matrix<-FALSE
  #   warning("impossible to create symmetrical matrix, number of matrix columns is different from number of rows")
  # }
  if (class(mPT)=="multiPermTest"){     
    mat<-matMultiPermTest(mPT,zs.type=zs.type)      
    fit <- pvclust(mat, method.hclust=hc.method, method.dist=dist.method,
                   nboot = nboot,...)
    ind<-fit$hclust$order
    if (symm_matrix==TRUE){
      
      newNames<-colnames(mat)[ind]
      mat<-mat[newNames,newNames]
      nc<-pamk(mat)$nc 
      clus <- kmeans(mat, centers=nc)
      
    }
    mat1<-list(GMat=mat,GFit=fit,GKmean=clus)
    class(mat1)<-"GenomicMatrix"
    return(mat1)
  }
  
}


plotGenomeMatrix<-function(GenMat,
                             graph.type, tl.col= "black",tl.srt=45, colMatrix="default",
                             tl.cex = 0.5, pch.col ="black",cl.lim = c(-1,1),
                             nc=NULL,
                             color=TRUE, shade=TRUE, labels=2, lines=0,
                             alpha=.95, lwd=2 , pv="au", border="red") {
    
    graph.type <- match.arg(graph.type, c("matrix", "pvclust", "clusplot","all"))
    
    if (class(GenMat)!="GenomicMatrix"){stop("the input is not a GenoMatrix Object")}
    
    if (!hasArg(GenMat)) {stop("A is missing")}
    
    paletteMatrix<-colorRampPalette(c("#67001F", "#B2182B", "#D6604D",  "#F4A582", 
                                      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                      "#4393C3", "#2166AC", "#053061"))
    
    if(colMatrix=="default") {colMatrix<-rev(paletteMatrix(50))}
    
    if (is.null(nc)){
      set.seed(123)
      nc<-pamk(GenMat$GMat)$nc}
    clus <- kmeans(GenMat$GMat, centers=nc)
    
    if (graph.type=="matrix" | graph.type=="all" ){
      
      corrA<-corrplot(GenMat$GMat, tl.col = tl.col, 
                      tl.srt = tl.srt, is.corr = F,
                      col = colMatrix, tl.cex = tl.cex,
                      pch.col=pch.col, cl.lim =  cl.lim)
      
    }
    
    if (graph.type=="clusplot" | graph.type=="all" ){
    
      clusplot(GenMat$GMat, clus$cluster, color = color, shade = shade, 
               labels = labels, lines = lines)
      
    }
    
    if (graph.type=="pvclust" | graph.type=="all" ){
   
      plot(GenMat$GFit)
      pvrect(GenMat$GFit, alpha=alpha, lwd=lwd ,pv=pv, border = "red")
      
    }
    
  }


createUniverse<-function(Alist){
  uniList<-GRanges()
  for(u in 1:length(Alist)){
    uniList<-c(uniList,Alist[[1]])
  }
  return(uniList)
}

TOT<-c(RSL1,RSL2,RSL3,RSL4)
a<-sample(length(TOT),4)
b<-sample(length(TOT),7)

mpt<-multiPermTest(Alist = TOT[a],
                   Blist = TOT[b],
                   genome = alienGenome,subEx = 0,max_pv=1)
mat<-matMultiPermTest(mpt)

mPT=multiOPT_circular
    
fit <- pvclust(mat, method.hclust="ward.D",method.dist="euclidean",
               nboot = 1000)
ind<-fit$hclust$order

newNames<-colnames(mat)[ind]
mat<-mat[,newNames]

fit1 <- pvclust(t(mat), method.hclust="ward.D",method.dist="euclidean",
               nboot = 1000)

ind1<-fit1$hclust$order
newNames1<-row(mat)[ind1]
mat<-mat[newNames1,]

nc<-pamk(mat)$nc 
clus <- kmeans(mat, centers=nc)
clus1 <- kmeans(t(mat), centers=nc)

a<-clus$cluster*1000
b<-clus1$cluster
mat<-mat1<-vector()
for(i in 1:length(b)){mat<-rbind(mat,a)}
for(i in 1:length(b)){mat1<-cbind(mat1,b)}

unique(as.vector(mat+mat1))
