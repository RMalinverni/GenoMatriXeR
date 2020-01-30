library(regioneR)
library(AnnotationHub)
hg19<-filterChromosomes(getGenome("hg19"))
hg19
ah <- AnnotationHub()
qh <- query(ah, c("hg19", "genes","Gencode v32"))
qh$coordinate_1_based


smurfenize<-toGRanges(data.frame(chr=c("chr1","chr7","chr13"),
                     start=c(155110773,97336314,79844644),
                     end=c(155110773,97336314,79844644)+10e6))

chrs<-seqlevels(smurfenize)
for(i in 1:length(chrs)){
  genes[seqnames(genes)==chrs[i]]
  
}

genes<-qh[[which(qh$title == "Annotated genes for Gencode v32 on hg19 coordinates")]]

subtractRegions(smurfenize,genes)
a<-overlapRegions(genes,smurfenize,type="AinB",colA="all")




toSmurfGenome