
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


randomizeRegionsPerc<-function(GR,frac=0.2){
  nc<-round(length(GR)*frac)
  change<-sample(length(GR),nc)
  GR1<-GR[-change]
  GR2<-randomizeRegions(GR[change],genome = "hg19")
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
  }else{ranged_x<-vec}
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

funRemove<-function(x,max_pv=0.05){     # funzione per azzerare gli zScore che non passano in test 
  x$z_score[x$adj.p_value>max_pv]<-0        # da portare fuori
  x$norm_zscore[x$adj.p_value>max_pv]<-0
  x$ranged_zscore[x$adj.p_value>max_pv]<-0
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

