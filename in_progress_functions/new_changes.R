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

mPt_smp<-makeGenomicMatrix(mPt_smp,zs.type = 'norm_zscore')

mPt_norm<-makeGenomicMatrix(mPt2,zs.type = 'norm_zscore')
mPt_std<-makeGenomicMatrix(mPt2,zs.type = 'std_zscore')
plotGenomeMatrix(mPt_ranged,nc = 5)
mPt_ranged<-makeGenomicMatrix(mPt,zs.type = 'ranged_zscore')
save(mPt_smp,file="/mnt/data1/GenoMartrixeR_folder/regionSets/TF_ENCODE/TF_Hepg2_mPt_resample.RData")
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


corrA<-corrplot::corrplot(corr=t(GM2),tl.col= "black",tl.srt=45,  method = "circle",
                is.corr = TRUE,type="full",order="original", col=colXX,
                tl.cex = 0.5, pch.col ="black")
corrplot::corrplot(GM,col=colXX,is.corr = FALSE)
memoise::forget()
