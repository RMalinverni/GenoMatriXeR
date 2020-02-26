library(regioneR)
library(AnnotationHub)
library(lattice)
library(forcats)
library(RColorBrewer) 
library(miceadds)
library(corrplot)
library(ggplot2)


ah<-AnnotationHub()

# qhs<-query(ah,c("ENCODE","broadPeak","E118-"))
# HM_HepG2_peaks<-list()
# for(i in 1:length(qhs$title)){
#   HM_HepG2_peaks[[i]] <- qhs[[qhs$ah_id[i]]]
# }
paletteMatrix<-colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                   "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                   "#4393C3", "#2166AC", "#053061"))
paletteLZMatrix<-colorRampPalette(c("red", "black","green"))


setwd("/home/malli/GenoMatriXeR/R")
for(i in list.files()){source(i)}
source("/home/malli/GenoMatriXeR/accessiories_functions/new_fun_genomaTrixer.R")
source("/home/malli/GenoMatriXeR/accessiories_functions/rquery_cormat-copy.r")

####################Calculate matrix (pt and lz for HepG2)
lovPT_HepG2<-listOverlapPermTest(Alist = list_HepG2_CS,
                                 Blist = CellPeakList,
                                 sampling=TRUE)

mat_HepG2<-matListOverlap(lovPT_HepG2)
mat_HepG2<-clustMatrix(mat_HepG2)
plotMatrix(mat_HepG2)

llzPT_HepG2<-lZAssociations(A=list_HepG2_CS$`7_Enh`,sampling = TRUE,
                      Blist=CellPeakList,
                      window=10000,
                      step=500)

lzmat_HepG2<-matLZAssociations(llzPT_HepG2)
lzmat_HepG2<-clustMatrix(lzmat_HepG2,lZ_tab = TRUE)


plot(lzmat_K562[1,],type='n',ylim=c(min(lzmat_K562),max(lzmat_K562)))
for(i in 1:nrow(lzmat_K562)){
  lines(lzmat_K562[i,])
}

max(lzmat_K562)

landscapePlot(lzmat=lzmat_CTCF_K562,rotation=1,
              title="CTCF K562", lim = c(-7e3,7e3),
              subtitle=" TFs - K562",
              xlab="bp",ylab="z-score")

corrplot(mat, tl.col="black", 
         tl.srt=45,is.corr = F,
         col=rev(paletteMatrix(50)),
         tl.cex = 0.5,
         pch.col="black")


corrplot(as.matrix(RedMat), tl.col="black", 
         tl.srt=45,is.corr = F,
        col=rev(paletteLZMatrix(50)),
         tl.cex = 0.5,
         pch.col="black",bg = "black")


llzPT_K562<-lZAssociations(A=HM_peaks_K562$`E123-H3K4me1.broadPeak.gz`,
                           Blist=list_K562_CS,
                           sampling=TRUE,
                           per.chromosome=TRUE,
                           fraction = 0.1,
                           window=10000,
                           step=500)

mat1<-matLZAssociations(llzPT_K562)


corrplot(clustMatrix(mat1,lZ_tab = T), tl.col="black", 
         tl.srt=45,is.corr = F,
         col=rev(paletteLZMatrix(50)),
         tl.cex = 0.5,
         pch.col="black",bg = "black")

#####################################################33333

load("/home/malli/work/data/regionSets/TF_ENCODE/TF_HepG2_peaks.RData")
load("/home/malli/work/data/tables_GenoMatriXeR/tab_HepG2_GMXR_CRR_F010_nt25.RData")




TF_HepG2_peaks<-cleanNames(TF_HepG2_peaks,"Hepg2")

mat<-matListOverlap(PTList_HepG2)
mat<-clustMatrix(mat,mirrored=TRUE)



corrplot(mat, tl.col="black", 
         tl.srt=45,is.corr = F,
         col=rev(paletteMatrix(50)),
         tl.cex = 0.5,
         pch.col="black")


mat['58-Mxi1',"41-Brca1a300"]     #58- son 20371 il 10% 2037
mat["41-Brca1a300",'58-Mxi1']     #41- son 1497  il 10% son 149



A=TF_HepG2_peaks$`41-Brca1a300`
pt<-permTest(A=A,
         B=TF_HepG2_peaks$`58-Mxi1`,
         randomize.function = randomizeRegions,
         evaluate.function = numOverlaps,
         ntimes=25,
         count.once=TRUE)

ANzscore=pt$numOverlaps$zscore/sqrt(length(A))

A=TF_HepG2_peaks$`41-Brca1a300`
ind<-sample(1:length(A),round(length(A))*0.1)
A=A[ind]
pt<-permTest(A=A,
             B=TF_HepG2_peaks$`58-Mxi1`,
             randomize.function = randomizeRegions,
             evaluate.function = numOverlaps,
             ntimes=25,
             count.once=TRUE)

orig.ev<-pt$numOverlaps$observed
rand.ev<-pt$numOverlaps$permuted
zscore <- round((orig.ev - mean(rand.ev, na.rm = TRUE))/stats::sd(rand.ev,na.rm = TRUE), 4)
maxZscore <- round((length(A) - mean(rand.ev, na.rm = TRUE))/stats::sd(rand.ev,na.rm = TRUE), 4)
normZS<-zscore/sqrt(length(A))
normZSmax<-maxZscore/sqrt(length(A))
stZscore<-normZS/normZSmax

ANzscore=pt$numOverlaps$zscore/sqrt(length(ind))

ind<-sample(length(A),round(length(A)*0.10))
A=TF_HepG2_peaks$`58-Mxi1`
B=TF_HepG2_peaks$`41-Brca1a300`



A=TF_HepG2_peaks$`58-Mxi1`
B=TF_HepG2_peaks$`41-Brca1a300`
C=TF_HepG2_peaks$`74-Cmyc`
DF_1<-testAssociations(A=A,B=B,ntimes = 200)
DF_X<-testAssociations(A=B,B=B,ntimes = 200)


hist(width(A))




DF_2<-testAssociations(A=B,B=A)
DF_3<-testAssociations(A=A,B=C)
DF_4<-testAssociations(A=C,B=A)



A=TF_HepG2_peaks$`58-Mxi1`
B=TF_HepG2_peaks$`41-Brca1a300`
pt<-permTest(A=A,
             B=B,
             randomize.function = circularRandomizeRegions,
             evaluate.function = numOverlaps,
             ntimes=ntimes,
             mc.cores=4,
             count.once=TRUE)





plot(DF_1$Nreg,DF_1$zcore,type="l")
lines(DF_1$Nreg,DF_2$zcore,col="red")
plot(DF_3$Nreg,DF_3$zcore,type="l")
lines(DF_3$Nreg,DF_4$zcore,col="red")

plot(DF_1$Nreg,DF_1$norm_zscore,type="l")
lines(DF_1$Nreg,DF_2$zcore,col="red")

plot(DF_1$Nreg,DF_1$n_ov,type="l")
lines(DF_1$Nreg,DF_1$zcore,col="red")
plot(DF$Nreg,DF$norm_zscore,type="l")
plot(DF$Nreg,DF$stZscore,type="l",ylim=c(-1,1))
plot(DF$Nreg,DF$n_ov,type="l")


x<-lOpT$`73-Znf274`
funRemove<-function(x){
  x$z_score[x$p_value>0.05]<-0
  x$norm_zscore[x$p_value>0.05]<-0
  return(x)
}
lOpT<-lapply(lOpT,funRemove)
mat<-matListOverlap(lOpT,zs.type = "norm")
mat[is.infinite(mat)]<-0
mat<-clustMatrix(mat,mirrored = TRUE)

rangedVector <- function(x){
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

set.seed(136)
samp20<-sample(77,20)

fake_genome<-collapse_list(TF_HepG2_peaks[samp20])
antiGenome<-subtractRegions(genome,fake_genome)
fake_genome<-filterChromosomes(extendRegions(fake_genome,extend.start = 1e3,extend.end = 1e3))
genome<-filterChromosomes(getGenome("hg19"))

antiRegion<-createRandomRegions(nregions=20000, length.mean=300, length.sd=150, genome=antiGenome, mask=NULL, non.overlapping=TRUE)

sum(width(fake_genome))

gg<-getGenome(antiGenome)
overlapRegions(antiRegion,antiGenome)
antiRegion<-antiRegion[!overlapRegions(antiRegion,fake_genome,only.boolean = T)]

matCorr<-listRGBinTable(TF_HepG2_peaks[samp20])
mat1<-rquery.cormat(as.matrix(mcols(matCorr)),type = "full")

seqlevels(gg)[1]



testPeaks<-TF_HepG2_peaks[samp20]
testPeaks$TestAntiPeaks<-antiRegion
lOpT_withTest<-listOverlapPermTest2(Alist =testPeaks,   #new rule problem (check function comment)
                           Blist = testPeaks,
                           sampling = T,
                           genome=fake_genome,
                           ranFun = "randomizeRegions",
                           mc.cores = 4) 
load("/home/malli/work/data/regionSets/HepG2_list.BroadPeak.RData")

lOpT_HM<-listOverlapPermTest2(Alist =totalList,   #new rule problem (check function comment)
                                    Blist = totalList,
                                    sampling = T,
                                    genome="hg19",
                                    ranFun = "randomizeRegions",
                                    mc.cores = 4) 

mat<-matListOverlap(lOpT_HM)
mat<-clustMatrix(mat,mirrored = TRUE)

load("/home/malli/work/data/tables_GenoMatriXeR/Matrices_HM_CvsP.RData")
names(CellPeakList)
HepG2_HM<-cleanNames(CellPeakList,vecEx = c("broadPeak"),cellName = "-")
TF_HepG2_peaks<-cleanNames(TF_HepG2_peaks,cellName = "Hepg2")
totalList<-c(TF_HepG2_peaks[samp20],HepG2_HM)

corrplot(Matrix_Hepg2_HM_Permutation, tl.col="black",
         tl.srt=45,
         col=rev(paletteMatrix(50)),
         tl.cex = 1,
         pch.col="black")


# Matrix_Hepg2_HM_Permutation<-mat
# Matrix_Hepg2_HM_Correlation<-as.matrix(mcols(mat1))
# save(Matrix_Hepg2_HM_Permutation,
#      Matrix_Hepg2_HM_Correlation,
#      file="/home/malli/work/data/tables_GenoMatriXeR/Matrices_HM_CvsP.RData")


mat1<-listRGBinTable(totalList)
rquery.cormat(Matrix_Hepg2_HM_Correlation,type="full")
rquery.cormat(as.matrix(mcols(mat1)),type="full")

vec<-c(32,0,Inf,34,0,1)
rangedVector(vec)

load('/home/malli/work/data/regionSets/HepG2_list.BroadPeak.RData')
load('/home/malli/work/data/regionSets/K562_list.BroadPeak.RData')



HM_HepG2_peaks<-cleanNames(HM_HepG2_peaks,vecEx = c(".broad"),cellName = "-")
HM_K562_peaks<-cleanNames(HM_K562_peaks,vecEx = c(".broad"),cellName = "-")

HM_K562_peaks<-HM_K562_peaks[-(grep("H3K9me1",names(HM_K562_peaks)))]

setwd("/home/malli/work/data/regionSets/")
load("PTlists.RData")


# PTList_HepG2_HM<-listOverlapPermTest2(Alist = HM_HepG2_peaks,
#                                       Blist = HM_HepG2_peaks,
#                                       ranFun = "randomizeRegions",
#                                       sampling = TRUE,
#                                       fraction = 0.15)
# 
# PTList_K562_HM<-listOverlapPermTest2(Alist = HM_K562_peaks,
#                                     Blist = HM_K562_peaks,
#                                     ranFun = "randomizeRegions",
#                                     sampling = TRUE,
#                                     fraction = 0.15)
  PTList_K562_HM<-listOverlapPermTest2(Alist = HM_HepG2_peaks,
                                     Blist = HM_HepG2_peaks,
                                     ranFun = "resampleRegions",
                                     sampling = TRUE,
                                     fraction = 0.15, verbose = TRUE)

 PTList_K562_HM<-listOverlapPermTest2(Alist = HM_K562_peaks,
                                     Blist = HM_K562_peaks,
                                     ranFun = "resampleRegions",
                                     sampling = TRUE,
                                     fraction = 0.15)


 
 
mat_HepG2<-clustMatrix(matListOverlap(PTList_HepG2_HM),mirrored = TRUE)
mat_K562<-clustMatrix(matListOverlap(PTList_K562_HM),mirrored = TRUE)




colnames(a)[order(colnames(a))]
colnames(b)[order(colnames(b))]
colnames(mat_HepG2)<-gsub("^(.*)-","",colnames(mat_HepG2)) 
rownames(mat_HepG2)<-gsub("^(.*)-","",rownames(mat_HepG2)) 
colnames(mat_K562)<-gsub("^(.*)-","",colnames(mat_K562)) 
rownames(mat_K562)<-gsub("^(.*)-","",rownames(mat_K562)) 

ind<-order(colnames(mat_HepG2))
mat_HepG2<-mat_HepG2[ind,ind]
ind<-order(colnames(mat_K562))
mat_K562<-mat_K562[ind,ind]


barplot(unlist(lapply(HM_K562_peaks,function(a){max(width(a))})))
barplot(unlist(lapply(HM_HepG2_peaks,function(a){max(width(a))})))





a<-corrplot(mat_HepG2, tl.col="black", 
         tl.srt=45,is.corr = F,
         col=rev(paletteMatrix(50)),
         tl.cex = 0.8,
         pch.col="black",
         cl.lim = c(-1,1))

b<-corrplot(mat_K562, tl.col="black", 
            tl.srt=45,is.corr = F,
            col=rev(paletteMatrix(50)),
            tl.cex = 0.8,
            pch.col="black",
            cl.lim = c(-1,1),title = "K562-Histone Modifications")

corrplot(clustMatrix(a-b,mirrored = TRUE), tl.col="black", 
            tl.srt=45,is.corr = F,
            col=rev(paletteMatrix(50)),
            tl.cex = 0.8,
            pch.col="black",cl.ratio=0.1,cl.lim = c(-1,1))


corrMatA<-listRGBinTable(HM_HepG2_peaks)
mat1<-rquery.cormat(as.matrix(mcols(corrMatA)),type = "full")



corrMatB<-listRGBinTable(HM_K562_peaks)
mat2<-rquery.cormat(as.matrix(mcols(corrMatB)),type = "full")








