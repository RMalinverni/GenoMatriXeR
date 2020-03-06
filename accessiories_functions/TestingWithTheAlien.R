library(regioneR)
library(AnnotationHub)
library(lattice)
library(forcats)
library(RColorBrewer) 
#library(miceadds)
library(corrplot)
library(ggplot2)
library(cluster)
library(clustree)
library(pvclust)
library(fpc)
#library("Seurat")


alienGenome<-toGRanges(data.frame(chr=c("aln1","aln2","aln3"),start=c(1,1,1),end=c(5e6,4e6,4e6)))
set.seed(336)
RS1<-createRandomRegions(nregions = 230,length.mean = 240,length.sd = 100,genome = alienGenome)
set.seed(431)
RS2<-createRandomRegions(nregions = 210,length.mean = 130,length.sd = 200,genome = alienGenome)
set.seed(231)
RS3<-createRandomRegions(nregions = 250,length.mean = 270,length.sd = 10,genome = alienGenome)
set.seed(131)
RS4<-createRandomRegions(nregions = 220,length.mean = 310,length.sd = 3,genome = alienGenome)



createRandomRegionSets<-function(RegionSet,frac=0.5,seedName='',seed=123,nregions=10,genome="hg19",...){
  set.seed(seed)
  lista<-list()
  for(i in 1:nregions){
    if (i==1){
      lista[[i]]<-randomizeRegionsPerc(RegionSet,frac=frac,genome=genome)
    }else{
      lista[[i]]<-randomizeRegionsPerc(lista[[i-1]],frac=frac,genome=genome)
    }
    names(lista)[[i]]<-paste0(seedName,i)
  }
  return(lista)
}

RSL1<-createRandomRegionSets(RS1,seedName = "RSL1_",nregions = 5,genome=alienGenome,frac=0.2)
RSL2<-createRandomRegionSets(RS2,seedName = "RSL2_",nregions = 10,genome=alienGenome,frac=0.7)
RSL3<-createRandomRegionSets(RS3,seedName = "RSL3_",nregions = 7,genome=alienGenome,frac=0.3)
RSL4<-createRandomRegionSets(RS4,seedName = "RSL4_",nregions = 3,genome=alienGenome,frac=0.1)


multiOPT<-multiOverlapPermTest(Alist=c(RSL1,RSL2,RSL3,RSL4),sampling = FALSE,ranFun = "resampleRegions",
                               verbose=TRUE,genome=alienGenome,ntimes=100,
                               max_pv = 0.05,pv="p.value",subEx=0)

multiOPT_circular<-multiPermTest(Alist=c(RSL1,RSL2,RSL3,RSL4),sampling = FALSE,ranFun = "circularRandomizeRegions",
                               verbose=TRUE,genome=alienGenome,ntimes=100,
                               max_pv = 0.05,pv="p.value",subEx=0)


mat_circular<-makeGenomicMatrix(mPT=multiOPT_circular)
plotGenomeMatrix(mat_circular,graph.type = "all",cl.lim=c(0,1))


makeGenomicMatrix<-function(multiOPTL.obj,hc.method="ward.D",dist.method="euclidean",
                            nboot=1000,zs.type='ranged_zscore',...){
  
  
  if (class(multiOPTL.obj)=="multiOverlapPermTestList"){     #this work fro simmetric matrix
      mat<-matListOverlap(multiOPTL.obj,zs.type=zs.type)      
      fit <- pvclust(mat, method.hclust="ward.D2", method.dist="euclidean",
                     nboot = nboot,...)
      ind<-fit$hclust$order
      newNames<-colnames(mat)[ind]
      mat<-mat[newNames,newNames]
      nc<-pamk(mat)$nc 
      clus <- kmeans(mat, centers=nc)
      mat1<-list(GMat=mat,GFit=fit,GKmean=clus)
  class(mat1)<-"GenomicMatrix"
  return(mat1)
  }
  
}

class(mat_circular)




corrA<-corrplot(mat_circular$GMat, tl.col="black", 
            tl.srt=45,is.corr = F,
            col=rev(paletteMatrix(150)),
            tl.cex = 0.8,
            pch.col="black",
            #cl.lim = c(-1,1),
            )

corrB<-corrplot(mat_resalmping$GMat, tl.col="black", 
                tl.srt=45,is.corr = F,
                col=rev(paletteMatrix(150)),
                tl.cex = 0.8,
                pch.col="black",
                #cl.lim = c(-1,1),
)


clusplot(mat$GMat, clus$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

plot(mat$GFit)
pvrect(mat$GFit, alpha=.90,lwd=3,pv ="au",border = "green")
pvrect(mat$GFit, alpha=.95,lwd=3,pv ="au",border = "blue") 
pvrect(mat$GFit, alpha=.99,lwd=3,pv ="au",border = "red") 
 
nc<-pamk(mat_circular$GMat)$nc 
clus <- kmeans(mat_circular$GMat, centers=nc)
clusplot(mat_circular$GMat, clus$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)


nc<-pamk(mat_resalmping$GMat)$nc 
clus <- kmeans(mat_resalmping$GMat, centers=nc)
clusplot(mat_resalmping$GMat, clus$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)


plot(mat$GFit)
pvrect(mat$GFit, alpha=.95,lwd=3,pv ="bp",border = "green")
pvrect(mat$GFit, alpha=.98,lwd=3,pv ="bp") 

d<-dist(mat)
hc<-hclust(d)

d <- dist(mat, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D2")

groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
plot(fit) 
rect.hclust(fit, k=5, border="red") 
ind<-fit$order
newNames<-colnames(mat)[ind]

mat[newNames,rev(newNames)]

colnames(mat[ind,rev(rev(ind))])
rownames(mat[ind,rev(rev(ind))])

fit<-kmeans(mat$GMat,2)
plot(mat$GFit)
pvrect(fit, alpha=.89) 
clus <- kmeans(mat$GMat, centers=2)
clusplot(mat$GMat, clus$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)


library(clustree)
kmeans(mat$GMat,)

K1<-kmeans(mat$GMat, 1)$cluster
K2<-kmeans(mat$GMat, 2)$cluster
K3<-kmeans(mat$GMat, 3)$cluster
K4<-kmeans(mat$GMat, 4)$cluster
K5<-kmeans(mat$GMat, 5)$cluster
K6<-kmeans(mat$GMat, 6)$cluster
K10<-kmeans(mat$GMat, 10)$cluster
PC1<-prcomp(mat$GMat)$rotation[,1]
PC2<-prcomp(mat$GMat)$rotation[,2]
objK<-data.frame(K1=K1,K2=K2,K3=K3,K4=K4,K5=K5,K6=K6,PC1,PC2,group=substr(names(K1),1,4))

pamk(mat$GMat)$nc 

substr(names(K1),1,4)
library("Seurat")
clustree(objK, prefix = "K")
clustree_overlay(objK, prefix = "K", x_value = "PC1", y_value = "PC2",groups="group",label_nodes = TRUE)
overlay_list <- clustree_overlay(objK, prefix = "K", x_value = "PC1",
                                 y_value = "PC2", plot_sides = TRUE)
overlay_list$x_side
overlay_list$y_side

as.data.frame(fit$cluster)
seurat <- CreateSeuratObject(counts = mat$GMat,
                             meta.data = sc_example$seurat_clusters)
# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster)
fit <- kmeans(mat$GMat, 3)
a<-clusplot(mat$GMat, fit$cluster, color=TRUE, shade=FALSE,
         labels=3, lines=0)
plot(mat$GFit)
pvrect(mat$GFit, alpha=.99) 



# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(mat$GMat, fit$cluster) 






b<-corrplot(mat$GMat, tl.col="black", 
            tl.srt=45,is.corr = F,
            col=rev(paletteMatrix(50)),
            tl.cex = 0.8,
            pch.col="black",
            cl.lim = c(-1,1))


library(pvclust)


