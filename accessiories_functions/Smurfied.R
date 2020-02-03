library(regioneR)
library(AnnotationHub)
load("/home/malli/work/data/regionSets/HepG2_list.BroadPeak.RData")
load("/home/malli/work/data/regionSets/HepG2_CS.RData")

smurfenize<-toGRanges(data.frame(chr=c("chr1","chr1","chr7","chr13"),
                                 start=c(155110773,50000000,97336314,79844644),
                                 end=c(155110773,50000000,97336314,79844644)+30e6))


a<-HepG2_CS[HepG2_CS$name=="Active TSS"]
b<-CellPeakList$`E118-H3K27ac.broadPeak.gz`

a1<-joinRegions(a[overlapRegions(a,smurfenize,only.boolean = T)])
b1<-b[overlapRegions(a,smurfenize,only.boolean = T)]
ov<-overlapPermTest(a1,b1,ntimes=230,genome=smurfenize)
ov1<-overlapPermTest(a,b,ntimes=230,genome=smurfenize)###???


lz<-localZScore(A = a1,B=b1,pt=ov,window = 10e4,step = 500)




hg19<-filterChromosomes(getGenome("hg19"))
hg19
ah <- AnnotationHub()
qh <- query(ah, c("hg19", "genes","Gencode v32"))

qh$coordinate_1_based
genes<-qh[[which(qh$title == "Annotated genes for Gencode v32 on hg19 coordinates")]]

smurfenize<-toGRanges(data.frame(chr=c("chr1","chr1","chr7","chr13"),
                     start=c(155110773,50000000,97336314,79844644),
                     end=c(155110773,50000000,97336314,79844644)+10e6))

GR<-subtractRegions(hg19,smurfenize)




chrs<-seqlevels(smurfenize)
for(i in 1:length(chrs)){
  overlapRegions(genes[seqnames(genes)==chrs[i]],smurfenize,
                 type="AinB",colA='Gene')
  
}



subtractRegions(smurfenize,genes)
a<-overlapRegions(genes,smurfenize,type="AinB",colA="all")




toSmurfGenome