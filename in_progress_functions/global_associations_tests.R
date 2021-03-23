source("/home/malli/GenoMatriXeR/accessiories_functions/new_fun_genomaTrixer.R")
source("/home/malli/GenoMatriXeR/accessiories_functions/rquery_cormat-copy.r")

listGR<-function(GR,colname){
  nc<-which(colnames(mcols(GR))==colname)
  colnames(mcols(GR))[nc]<-"sel"
  uniHM<-unique(GR$sel)
  listHM<-list()
  for(i in 1:length(uniHM)){
    listHM[[i]]<- GR %>% filter(sel==uniHM[i])
  }
  names(listHM)<-uniHM
  return(listHM)
}


load("/home/malli/work/data/regionSets/HepG2_list.BroadPeak.RData")
load("/home/malli/work/data/regionSets/TF_ENCODE/TF_HepG2_peaks.RData")
load("/home/malli/work/data/regionSets/HepG2_CS.RData")


library(pvclust)
library(cluster)
library(fpc)
library(corrplot)

macro_HepG2<-toGRanges(read.delim("/home/malli/work/insulation/peaks/peaks_HepG2_epic.bed"))
AL<-c(HM_HepG2_peaks,macroH2A=macro_HepG2)
BL<-listGR(HepG2_CS,"name")

TF_HepG2_peaks<-cleanNames(TF_HepG2_peaks,"Hepg2")

mOpT<-multiOverlapPermTest(Alist = AL,   # Inf problem (remember)
                     Blist = TF_HepG2_peaks,     
                     sampling = TRUE,
                     fraction =  0.1,
                     genome="hg19",
                     ranFun = "circularRandomizeRegions",
                     ntimes=50,
                     mc.cores = 1,
                     subEx=0,per.chromosome=TRUE) 

a<-makeGenomicMatrix(mOpT,symm_matrix = FALSE,zs.type = 'norm_zscore',transform = TRUE)
plotGenomeMatrix(a,graph.type = "matrix",nc=7,tl.cex = 1)


options(digits=3)
set.seed(20000)
face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
pk1 <- pamk(face,krange=1:5,criterion="asw",critout=TRUE)
pk2 <- pamk(face,krange=1:5,criterion="multiasw",ns=2,critout=TRUE)
# "multiasw" is better for larger data sets, use larger ns then.
pk3 <- pamk(face,krange=1:5,criterion="ch",critout=TRUE)


