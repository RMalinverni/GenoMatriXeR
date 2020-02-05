library(regioneR)
library(AnnotationHub)
library(lattice)
library(forcats)
library(RColorBrewer) 
library(miceadds)
source.all("/home/malli/GenoMatriXeR/R/*.R")
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

llzPT_HepG2<-lZAssociations(A=list_HepG2_CS$`7_Enh`,sampling = TRUE,
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

llzPT_K562<-lZAssociations(A=HM_peaks_HepG2$`E118-H3K4me3.broadPeak.gz`,
                            Blist=list_K562_CS,
                            sampling=TRUE,
                            per.chromosome=TRUE,
                            fraction = 0.1,
                            window=10000,
                            step=500)

lzmat_K562<-matLZAssociations(llzPT_K562)
lzmat_K562<-clustMatrix(lzmat_K562,lZ_tab = TRUE)
plotLZMatrix(lzmat_K562)


lzmat<-lzmat_HepG2
rotation<-300
my_palette<-colorRampPalette(palette_vec)
library(ggplot2)


landscapePlot(lzmat=lzmat_HepG2,rotation=-300)

landscapePlot<-function(lzmat,rotation=0,
                        palette_vec=c("brown","darkgreen","lightgreen","yellow","lightblue"),
                        theme_vec=c("lightblue","black","lightblue","white")){  






#ind<-order(as.numeric(apply(lzmat[,"0"],1,FUN=max)))
ind<-order(lzmat[,"0"])
lzmat<-lzmat[rev(ind),]
rownames(lzmat)<-paste0(letters[1:nrow(lzmat)],"_",rownames(lzmat))
rownames(lzmat)<-factor(rownames(lzmat), levels = rownames(lzmat)[order(rownames(lzmat))])
rownames(lzmat)<-as.factor(rownames(lzmat))
values<-as.vector(lzmat)
names<-rep(rownames(lzmat),ncol(lzmat))
steps<-vector()
for(i in 1:ncol(lzmat)){
  a<-rep(colnames(lzmat)[i],nrow(lzmat))
  steps<-c(steps,a)
  }


df<-data.frame(names=names,steps=steps,values=values)
df$steps<-as.numeric(as.character(df$steps))
lll<-length(unique(df$names))
changed<-seq(-(lll/2),lll/2,by=1)*rotation
for(i in 1:length(unique(df$names))){
  ndf<-as.character(unique(df$names)[i])
  df$steps[df$names==ndf]<-df$steps[df$names==ndf]+changed[i]
}



p<-ggplot(df, aes(x=as.numeric(as.character(df$steps)), 
          y=values,fill=names,group=as.factor(names))) + 
          geom_area(color="white") +
          coord_cartesian(xlim=c(-5000,5000)) +
          scale_fill_manual(values=my_palette(nlevels(as.factor(names))))+
          theme(panel.background = element_rect(fill = theme_vec[1],
                                                colour = theme_vec[2],
                                                size = 0.5, 
                                                linetype = "solid"),
                panel.grid.major = element_line(size = 0.5, 
                                                linetype = 'solid', 
                                                colour = theme_vec[3]), 
                panel.grid.minor = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = theme_vec[4])
          )

p
}

#scale_fill_brewer(palette="Spectral") +