downloadRegionsObject<-function(list){
  test<-list()
  for(i in 1:length(list)){
    download.file(list[i],"test.gz")
    test[i]<-toGRanges(read.delim(gzfile("test.gz")))
  }
  names(test)<-list
  system("rm test.gz")
  return(test)
}
