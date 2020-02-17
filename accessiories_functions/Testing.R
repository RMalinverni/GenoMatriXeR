library(regioneR)
library(AnnotationHub)
library(lattice)
library(forcats)
library(RColorBrewer) 
library(miceadds)
library(corrplot)
library(ggplot2)


ah<-AnnotationHub()


# q_ah<-query(ah,c("HepG2","TF","Uni"))
# TF_HepG2_peaks<-list()
# for(i in 1:length(names(q_ah))){
#   TF_HepG2_peaks[[i]]<-q_ah[[names(q_ah)[i]]]
# }
# names(TF_HepG2_peaks)<-q_ah$title
# save(TF_HepG2_peaks,file="/home/malli/work/data/regionSets/TF_HepG2_peaks.RData")

paletteMatrix<-colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                   "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                   "#4393C3", "#2166AC", "#053061"))
paletteLZMatrix<-colorRampPalette(c("red", "black","green"))


setwd("/home/malli/GenoMatriXeR/R")
for(i in list.files()){source(i)}

# list_HepG2_CS<-list()
# names_CS<-unique(HepG2_CS$abbr)
# for(i in 1:length(names_CS)){
#   list_HepG2_CS[[i]]<-HepG2_CS[HepG2_CS$abbr==names_CS[i]]
# }
# names(list_HepG2_CS)<-names_CS


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


cleanNames<-function(listRS,cellName,
                    vecEx=c("Uni","sc[0-9]","V0","Pcr","anb[0-9]","Iggrab","Forskln","Ucd")){
                    newNames<-gsub(paste0("^(.*)",cellName),"",names(listRS))
                    for(i in vecEx){
                      newNames<-gsub(paste0(i,"(.*)$"),"",newNames) 
                    }
                    names(listRS)<-newNames
                    listRS<-listRS[!is.na(names(listRS))]
                    names(listRS)<-paste0(1:length(names(listRS)),"-",names(listRS))
                    return(listRS)
}

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

testAssociations<-function(A,B,
                           seqFraction=seq(0.1,1,0.1),
                           ntimes=25){

  vecZs<-vecNZs<-vecSZs<-vecNZs1<-vecSZs1<-vecOE<-vecP<-vecPV<-vector()
  for(i in seq(0.1,1,0.1)){
    ind<-sample(length(A),round(length(A)*i))
    print(length(ind))
    A1=A[ind]
    pt<-permTest(A=A1,
                 B=B,
                 randomize.function = randomizeRegions,
                 evaluate.function = numOverlaps,
                 ntimes=ntimes,
                 mc.cores=4,
                 count.once=TRUE)
    orig.ev<-pt$numOverlaps$observed
    rand.ev<-pt$numOverlaps$permuted
    zscore <- round((orig.ev - mean(rand.ev, na.rm = TRUE))/stats::sd(rand.ev,na.rm = TRUE), 4)
    maxZscore <- round((length(A1) - mean(rand.ev, na.rm = TRUE))/stats::sd(rand.ev,na.rm = TRUE), 4)
    normZS<-zscore/sqrt(length(A1))
    normZSmax<-maxZscore/sqrt(length(A1))
    stZscore<-normZS/normZSmax
    
    maxZscore1 <- round((min(length(A1),length(B)) - mean(rand.ev, na.rm = TRUE))/stats::sd(rand.ev,na.rm = TRUE), 4)
    normZS1<-zscore/sqrt(min(length(A1),length(B)))
    normZSmax1<-maxZscore/sqrt(min(length(A1),length(B)))
    stZscore1<-normZS1/normZSmax1
    
    vecZs<-c(vecZs,zscore)
    vecNZs<-c(vecNZs,normZS)
    vecSZs<-c(vecSZs,stZscore)
    vecNZs1<-c(vecNZs1,normZS1)
    vecSZs1<-c(vecSZs1,stZscore1)

    vecOE<-c(vecOE,orig.ev)
    vecPV<-c(vecPV,pt$numOverlaps$pval)
    vecP<-c(vecP,mean(rand.ev, na.rm = TRUE))
  }
  
  DF<-DataFrame(frac=seq(0.1,1,0.1),
                Nreg=round(length(A)*seq(0.1,1,0.1)),
                zcore=round(vecZs,digits = 2),
                p_value=vecPV,
                norm_zscore=round(vecNZs,digits = 2),
                norm_zscore1=round(vecNZs1,digits = 2),
                stZscore=round(vecSZs,digits = 2),
                stZscore1=round(vecSZs1,digits = 2),
                n_ov=vecOE,
                menPerm_ov=round(vecP,digits = 2))
  return(DF)
}

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

samp20<-sample(77,20)
lOpT<-listOverlapPermTest2(Alist = TF_HepG2_peaks[samp20],
                          Blist = TF_HepG2_peaks[samp20],
                          sampling = T,
                          ranFun = "randomizeRegions",
                          mc.cores = 4)

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




mat<-matListOverlap(lOpT)
mat<-clustMatrix(mat,mirrored = TRUE)
corrplot(mat, tl.col="black", 
         tl.srt=45,is.corr = F,
         col=rev(paletteMatrix(50)),
         tl.cex = 0.5,
         pch.col="black")

vec<-c(32,0,Inf,34,0,1)
rangedVector(vec)
