

---
title:"Test GenoMatriXeR (GMXR)"
output: html_document
---

```{r}
library(regioneR)
library(AnnotationHub)
library(RColorBrewer) 
library(corrplot)
library(ggplot2)
library(rmarkdown)  

  
  
paletteMatrix<-colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                                  "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                  "#4393C3", "#2166AC", "#053061"))
paletteLZMatrix<-colorRampPalette(c("red", "black","green"))


setwd("/home/malli/GenoMatriXeR/R")
for(i in list.files()){source(i)}
source("/home/malli/GenoMatriXeR/accessiories_functions/new_fun_genomaTrixer.R")
source("/home/malli/GenoMatriXeR/accessiories_functions/rquery_cormat-copy.r")
#load list of GRangers HM for HepG2
load("/home/malli/work/data/regionSets/HepG2_list.BroadPeak.RData")
#load list of GRangers TF for HepG2
load("/home/malli/work/data/regionSets/TF_ENCODE/TF_HepG2_peaks.RData")

lOpT_HM<-listOverlapPermTest2(Alist = HM_HepG2_peaks,   # Inf problem (remember)
                              Blist = HM_HepG2_peaks,
                              sampling = T,
                              fraction =  0.01,
                              genome="hg19",
                              ranFun = "randomizeRegions",
                              mc.cores = 4) 

```



