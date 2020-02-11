library(regioneR)
library(AnnotationHub)
library(lattice)
library(forcats)
library(RColorBrewer) 
library(miceadds)
library(corrplot)
library(ggplot2)


ah<-AnnotationHub()
query(ah,c("K562","TF","Uni"))
query(ah,c("Gm12878","TF","Uni"))


q_ah<-query(ah,c("HepG2","TF","Uni"))
TF_HepG2_peaks<-list()
for(i in 1:length(names(q_ah))){
  TF_HepG2_peaks[[i]]<-q_ah[[names(q_ah)[i]]]
}
names(TF_HepG2_peaks)<-q_ah$title
save(TF_HepG2_peaks,file="/home/malli/work/data/regionSets/TF_HepG2_peaks.RData")

q_ah<-query(ah,c("K562","TF","Uni"))
TF_K562_peaks<-list()
for(i in 1:length(names(q_ah))){
  TF_K562_peaks[[i]]<-q_ah[[names(q_ah)[i]]]
}
names(TF_K562_peaks)<-q_ah$title
save(TF_K562_peaks,file="/home/malli/work/data/regionSets/TF_K562_peaks.RData")

q_ah<-query(ah,c("HCT116","TF","Uni"))
TF_HCT116_peaks<-list()
for(i in 1:length(names(q_ah))){
  print(paste0("downloading: ",i,"/",length(names(q_ah))))
  TF_HCT116_peaks[[i]]<-q_ah[[names(q_ah)[i]]]
}

names(TF_HCT116_peaks)<-q_ah$title
save(TF_HCT116_peaks,file="/home/malli/work/data/regionSets/TF_HCT116_peaks.RData")

gen_ah<-query(ah,c("GEO"))


source("http://www.sthda.com/upload/rquery_cormat.r")
paletteMatrix<-colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                   "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                   "#4393C3", "#2166AC", "#053061"))
paletteLZMatrix<-colorRampPalette(c("red", "black","green"))


setwd("/home/malli/GenoMatriXeR/R")
for(i in list.files()){source(i)}

load("/home/malli/work/data/regionSets/HepG2_list.BroadPeak.RData")
HM_peaks_HepG2<-CellPeakList
load("/home/malli/work/data/regionSets/K562_list.BroadPeak.RData")
HM_peaks_K562<-GRangesList(CellPeakList)
load("/home/malli/work/data/regionSets/HepG2_CS.RData")
load("/home/malli/work/data/regionSets/K562_CS.RData")
load("/home/malli/work/GenoMatrixeR_package/RData/listTF_K562.ENCODE.RData")
load("/home/malli/work/GenoMatrixeR_package/RData/tab_circular_TF_K562.ENCODE.RData")



K562_TF_circular_tab

list_HepG2_CS<-list()
names_CS<-unique(HepG2_CS$abbr)
for(i in 1:length(names_CS)){
  list_HepG2_CS[[i]]<-HepG2_CS[HepG2_CS$abbr==names_CS[i]]
}
names(list_HepG2_CS)<-names_CS


list_K562_CS<-list()
names_CS<-unique(K562_CS$abbr)
for(i in 1:length(names_CS)){
  list_K562_CS[[i]]<-K562_CS[K562_CS$abbr==names_CS[i]]
}
names(list_K562_CS)<-names_CS



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
####################Calculate matrix (pt and lz for K562)
smp<-sample(length(listTF_K562),10)
lovPT_K562<-listOverlapPermTest(Alist = listTF_K562[smp],
                                Blist = listTF_K562[smp],
                                sampling=TRUE,
                                mc.cores = 4,
                                ranFun = "circularRandomizeRegions",
                                fraction = 0.2)

mat_K562<-matListOverlap(lovPT_K562,zs.type = "std")
mat_K562<-clustMatrix(mat_K562,mirrored = T)
plotMatrix(mat_K562,type="normal")





llzPT_K562<-lZAssociations(A=HM_peaks_K562$`E123-H3K4me1.broadPeak.gz`,
                            Blist=K562_CS,
                            sampling=TRUE,
                            per.chromosome=TRUE,
                            fraction = 0.1,
                            window=10000,
                            step=500)

lzmat_K562<-matLZAssociations(llzPT_K562)
lzmat_K562<-clustMatrix(lzmat_K562,lZ_tab = TRUE)
plotLZMatrix(lzmat_K562)

lzmat_Corest_K562<-matLZAssociations(lZ_tab_Corest_TF_K562)
lzmat_Corest_K562<-clustMatrix(lzmat_Corest_K562,lZ_tab = TRUE)

rotation<-300
my_palette<-colorRampPalette(palette_vec)
library(ggplot2)

plot(lzmat_K562[1,],type='n',ylim=c(min(lzmat_K562),max(lzmat_K562)))
for(i in 1:nrow(lzmat_K562)){
  lines(lzmat_K562[i,])
}

max(lzmat_K562)

landscapePlot(lzmat=lzmat_CTCF_K562,rotation=1,
              title="CTCF K562", lim = c(-7e3,7e3),
              subtitle=" TFs - K562",
              xlab="bp",ylab="z-score")




tB<-listRGBinTable(listTF_K562[smp])
corrP<-rquery.cormat(as.matrix((mcols(tB))),graphType="heatmap")
corrP<-rquery.cormat(as.matrix((mcols(tB))))


corrplot(cormat, type=ifelse(type[1]=="flatten", "lower", type[1]),
         tl.col="black", tl.srt=45,col=col,...)


matGK562_TF_circular_tab

mat<-matListOverlap(K562_TF_circular_tab)
mat[mat>1]<-1
mat<-mat[1:100,3:50]
mat<-clustMatrix(mat)


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

llzPT_K562_2<-lZAssociations(A=HM_peaks_K562$`E123-H3K36me3.broadPeak.gz`,
                           Blist=list_K562_CS,
                           sampling=TRUE,
                           per.chromosome=TRUE,
                           fraction = 0.1,
                           window=10000,
                           step=500)


mat1<-matLZAssociations(llzPT_K562)
mat2<-matLZAssociations(llzPT_K562_2)

mat1<-lzmat_Corest_K562
mat2<-lzmat_CTCF_K562
mat1<-mat1[order(rownames(mat1)),]
mat2<-mat1[order(rownames(mat2)),]
matdiff<-mat1-mat2
matdiff<-clustMatrix(matdiff,lZ_tab = TRUE)


corrplot(clustMatrix(mat1,lZ_tab = T), tl.col="black", 
         tl.srt=45,is.corr = F,
         col=rev(paletteLZMatrix(50)),
         tl.cex = 0.5,
         pch.col="black",bg = "black")

corrplot(clustMatrix(mat2,lZ_tab = T), tl.col="black", 
         tl.srt=45,is.corr = F,
         col=rev(paletteLZMatrix(50)),
         tl.cex = 0.5,
         pch.col="black",bg = "black")

corrplot(matdiff, tl.col="black", 
         tl.srt=45,is.corr = F,
         col=rev(paletteLZMatrix(50)),
         tl.cex = 0.5,
         pch.col="black",bg = "black")


dim(mcols(a))
cor(mcols(a[,6])[[1]],mcols(a[,2])[[1]])
heatmap(mat_K562, col=paletteMatrix(200), symm=TRUE)



