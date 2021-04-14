library("regioneR")
library("pvclust")
library(cluster)
library(fpc)
library(corrplot)
library(RColorBrewer)

load("/mnt/data1/GenoMartrixeR_folder/regionSets/TF_ENCODE/TF_HepG2_peaks.RData")
setwd("/mnt/data1/RProjects/GenoMatriXeR/R")
lf<-list.files()
sapply(lf,source)
source("/mnt/data1/RProjects/GenoMatriXeR/accessiories_functions/new_fun_genomaTrixer.R")

TF_HepG2_peaks<-cleanNames(TF_HepG2_peaks,cellName = "Hepg2")
set.seed(123)
smp<-sample(1:length(TF_HepG2_peaks),10)
mPt_smp<-multiOverlapPermTest(Alist = TF_HepG2_peaks,
                          sampling = TRUE, fraction = 0.15, min_sampling = 2000,
                          ranFun = "resampleRegions",evFUN = "numOverlaps",
                          mc.cores = 12
                          )
mPt_smp_circ<-multiOverlapPermTest(Alist = TF_HepG2_peaks,
                              sampling = TRUE, fraction = 0.15, min_sampling = 2000,
                              ranFun = "circularRandomizeRegions",evFUN = "numOverlaps",
                              mc.cores = 12
)
mPt_smp_rand<-multiOverlapPermTest(Alist = TF_HepG2_peaks,
                              sampling = TRUE, fraction = 0.15, min_sampling = 2000,
                              ranFun = "randomizeRegions",evFUN = "numOverlaps",
                              mc.cores = 12,genome=posGen
)

mPt_smp_rand_small<-multiOverlapPermTest(Alist = TF_HepG2_peaks,
                                   sampling = TRUE, fraction = 0.15, min_sampling = 2000,
                                   ranFun = "randomizeRegions",evFUN = "numOverlaps",
                                   mc.cores = 12,genome=posGen2
)
mPt_smp_rand_xx<-multiOverlapPermTest(Alist = TF_HepG2_peaks,
                                         sampling = TRUE, fraction = 0.15, min_sampling = 2000,
                                         ranFun = "randomizeRegions",evFUN = "numOverlaps",
                                         mc.cores = 12,
)


save(mPt_smp_rand,file="/mnt/data1/GenoMartrixeR_folder/regionSets/TF_ENCODE/TF_Hepg2_mPt_randReg.RData")
save(mPt_smp_rand_small,file="/mnt/data1/GenoMartrixeR_folder/regionSets/TF_ENCODE/TF_Hepg2_mPt_randReg_small.RData")
save(mPt_smp_rand_total,file="/mnt/data1/GenoMartrixeR_folder/regionSets/TF_ENCODE/TF_Hepg2_mPt_randReg_total.RData")

mPt_smp_rand<-makeGenomicMatrix(mPt_smp_rand,zs.type = 'norm_zscore')
plotGenomeMatrix_bis(mPt_smp_rand,cl.lim = NULL,maxNZS=100,nc=5)

genome<-getGenomeAndMask("hg19")
validRanges<-collapse_list(TF_HepG2_peaks)
validRanges<-joinRegions(validRanges)
#validRanges<-extendRegions(validRanges,extend.start=1e3, extend.end=1e3)
negGen<-subtractRegions(genome$genome,validRanges)
posGen2<-subtractRegions(genome$genome,negGen)
posGen2<-subtractRegions(posGen2,genome$mask)
kp<-plotKaryotype()
kpPlotRegions(kp,posGen)  
  


kp<-plotKaryotype()
kpPlotRegions(kp,subtractRegions(genome$genome,genome$mask))
  
  
mPt_smp<-makeGenomicMatrix(mPt_smp,zs.type = 'std_zscore')
mPt_smp_circ<-makeGenomicMatrix(mPt_smp_circ,zs.type = 'norm_zscore')



mPt_norm<-makeGenomicMatrix(mPt2,zs.type = 'norm_zscore')
mPt_std<-makeGenomicMatrix(mPt2,zs.type = 'std_zscore')
plotGenomeMatrix(mPt_ranged,nc = 5)
mPt_ranged<-makeGenomicMatrix(mPt,zs.type = 'ranged_zscore')
save(mPt_smp_circ,file="/mnt/data1/GenoMartrixeR_folder/regionSets/TF_ENCODE/TF_Hepg2_mPt_circRand.RData")
#load("/mnt/data1/GenoMartrixeR_folder/regionSets/TF_ENCODE/TF_Hepg2_mPt.RData")
plotGenomeMatrix(mPt_norm,cl.lim = c(0,100),nc=5)
mPt<-makeGenomicMatrix(mPt,zs.type = 'std_zscore')


GM_norm<-mPt_norm@matrix$GMat
GM_std<-mPt_std@matrix$GMat
colXX<-c(rev(brewer.pal(9, "PuBuGn")),brewer.pal(9, "YlOrRd"))

mPtest<-mPt_smp
maxNZS<-4
GM2<-mPtest@matrix$GMat
mat_pv<-mPtest@matrix$GMat_pv
mat_pv<-mat_pv[rownames(GM2),colnames(GM2)]
GM2[mat_pv>0.05]<-0
GM2[GM2>=maxNZS]<-maxNZS
GMX<-GM2/abs(GM2)
GMX[is.nan(GMX)]<-1
GM2<-abs(GM2)
GM2<-log1p(GM2)
GM2<-GM2/max(GM2)
GM2<-GM2*GMX


corrA<-corrplot::corrplot(corr=t(GM3),tl.col= "black",tl.srt=45,  method = "circle",
                is.corr = FALSE,type="full",order="original", col=colXX,
                tl.cex = 0.5, pch.col ="black",main="test")

mPt_smp_circ@matrix$GMat

plotGenomeMatrix_bis(mPt_smp_rand_small,cl.lim = NULL,maxNZS=100,nc=5)

GM3<-GM2
GM3[GM3>0]<-1
GM3[GM3<0]<-(-1)

corrplot::corrplot(GM,col=colXX,is.corr = FALSE)
memoise::forget()


mPT=mPt_smp_circ
graph.type="all"
main=""
tl.col= "black"
tl.srt=45
colMatrix="default"
tl.cex = 0.5
pch.col ="black"
cl.lim = c(-1,1)
nc=NULL
color=TRUE
shade=TRUE
labels=2
lines=0
alpha=.95
lwd=2 
pv="au" 
border="red"
cex=0.7
maxNZS=20
