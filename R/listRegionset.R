listRegionset<-function(cell.line,type,folder){
  a<-read.delim(paste(folder,"files.txt",sep=""))
  a<-as.character(a[,1])
  a<-a[grep(type,a)]
  a<-a[grep(cell.line,a)]
  a<-paste(folder,a,sep="")
  return(a)
}



