

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


randomizeRegionsPerc<-function(GR,frac=0.2){
  nc<-round(length(GR)*frac)
  change<-sample(length(GR),nc)
  GR1<-GR[-change]
  GR2<-randomizeRegions(GR[change],genome = "hg19")
  GR3<-c(GR1,GR2)
  return(GR3)
}

