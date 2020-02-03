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

llzPT_K562<-lZAssociations(A=list_K562_CS$`7_Enh`,
                            Blist=CellPeakList,
                            window=10000,
                            step=500)

lzmat_K562<-matLZAssociations(llzPT_K562)
lzmat_K562<-clustMatrix(lzmat_K562,lZ_tab = TRUE)
plotLZMatrix(lzmat_K562)


lzmat<-lzmat_HepG2
plot(lzmat[,1],type="n",ylim=c(min(lzmat),max(lzmat)),xlim=c(1,ncol(lzmat)))
for(i in 1:nrow(lzmat)){
  lines(lzmat[i,])
}
library(ggplot2)



df<-data.frame(X=as.numeric(colnames(lzmat_K562)),
               Y=lzmat_K562[1,])





p <- ggplot(data=df, aes(x=X, y=Y)) + 
      geom_area(fill = "lightblue") +
      geom_line(color = "blue") 

time <- as.numeric(rep(seq(1,7),each=7))  # x Axis
value <- runif(49, 10, 100)               # y Axis
group <- rep(LETTERS[1:7],times=7)        # group, one shape per group
data <- data.frame(time, value, group)

values<-as.vector(lzmat_K562)
names<-rep(rownames(lzmat_K562),ncol(lzmat_K562))
steps<-vector()
for(i in 1:ncol(lzmat_K562)){
  a<-rep(colnames(lzmat_K562)[i],nrow(lzmat_K562))
  steps<-c(steps,a)
  }
steps<-rep(colnames(lzmat_K562),nrow(lzmat_K562))

df<-data.frame(names=names,steps=steps,values=values)
p<-ggplot(df, aes(x=steps, y=values)) + 
  geom_area(aes(fill = names)) 

p<-ggplot(df, aes(x=steps, y=values,group=names)) + 
  geom_line() 

df[seq(1,100,by=12),]




as.vector(lzmat_K562)
df<-as.data.frame(lzmat_K562)
p<-ggplot(df,aes(x=as.numeric(colnames(df)), y='E123-H3K4me1.broadPeak.gz')) +
  geom_area()

# stacked area chart
ggplot(data, aes(x=time, y=value, fill=group)) + 
  geom_area()

