library(regioneR)
library(AnnotationHub)
load("/home/malli/work/data/regionSets/HepG2_list.BroadPeak.RData")
HM_peaks_HepG2<-CellPeakList
load("/home/malli/work/data/regionSets/K562_list.BroadPeak.RData")
HM_peaks_K562<-CellPeakList
load("/home/malli/work/data/regionSets/HepG2_CS.RData")
load("/home/malli/work/data/regionSets/K562_CS.RData")


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

llzPT_HepG2<-lZAssociations(A=list_HepG2_CS$`7_Enh`,
                      Blist=CellPeakList,
                      window=10000,
                      step=500)

lzmat_HepG2<-matLZAssociations(llzPT_HepG2)
lzmat_HepG2<-clustMatrix(lzmat_HepG2,lZ_tab = TRUE)
####################Calculate matrix (pt and lz for K562)

lovPT_K562<-listOverlapPermTest(Alist = list_K562_CS,
                                 Blist = CellPeakList,
                                 sampling=TRUE)

mat_K562<-matListOverlap(lovPT_K562)
mat_K562<-clustMatrix(mat_K562)
plotMatrix(mat_K562)

llzPT_K562_s<-lZAssociations(A=list_K562_CS$`7_Enh`,
                            Blist=CellPeakList,
                            sampling=TRUE,
                            fraction=0.05,
                            window=10000,
                            step=500)

lzmat_K562<-matLZAssociations(llzPT_K562_s)
lzmat_K562<-clustMatrix(lzmat_K562,lZ_tab = TRUE)




plotLZMatrix(lzmat_K562)
plotLZMatrix(lzmat_K562,name='aaaa')

