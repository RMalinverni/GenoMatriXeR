#' Plot GenomeMatrix
#'
#' Create plots from a "GenomicMatrix" object
#'
#' @usage plotMatrix( (GenMat, graph.type = "all", main = "", tl.col = "black", tl.srt = 45, 
#' colMatrix = "default", tl.cex = 0.5, pch.col = "black", cl.lim = c(-1,1),
#' nc = NULL, color = TRUE, shade = TRUE, labels = 2, lines = 0,
#' alpha = .95, lwd = 2 , pv = "au", border = "red",cex = 0.7,...))
#'
#' @param 
#'
#' @param 
#' 
#' @param 
#' @param 
#' @param ...  further arguments to be passed to other methods.
#'
#'
#' @seealso ...
#'
#' @examples  ...
#'
#' #@references ...
#' #@importFrom ...
#' @export plotGenomeMatrix
#'
plotGenomeMatrix_bis<-function(mPT, graph.type="all",main="",
                           tl.col= "black",tl.srt=45, colMatrix="default",
                           tl.cex = 0.5, pch.col ="black",cl.lim = c(-1,1),
                           nc=NULL, color=TRUE, shade=TRUE, labels=2, lines=0,
                           alpha=.95, lwd=2 , pv="au", border="red",cex=0.7,...) {
  
  graph.type <- match.arg(graph.type, c("matrix", "pvclust", "clusplot","all"))
  
  if (class(mPT)!="GenoMatriXeR"){stop("the input is not a GenoMatrix Object")}
  
  if (!hasArg(mPT)) {stop("A is missing")}
  
  paletteMatrix<-c(rev(brewer.pal(9, "PuBuGn")),brewer.pal(9, "YlOrRd"))
  
  if (!is.null(mPT@matrix)){ 
    GM<-mPT@matrix$GMat 
    }else{stop("the matrix slot in the GenoMatriXeR is NULL")}
  
  if(colMatrix=="default") {colMatrix<-rev(paletteMatrix(50))}
  
  if (is.null(nc)){
    
    set.seed(123)
    nc<-pamk(t(GM),krange = 1:(ncol(GM)-1))$nc
    
  }
  
  clus <- kmeans(t(GM), centers=nc)
  
  if (graph.type=="matrix" | graph.type=="all" ){
    ind<-mPT@matrix$GFit$hclust$order
    
    maxNZS<-4
    GM2<-mPT@matrix$GMat
    mat_pv<-mPT@matrix$GMat_pv
    mat_pv<-mat_pv[rownames(GM2),colnames(GM2)]
    GM2[mat_pv>0.05]<-0
    GM2[GM2>=maxNZS]<-maxNZS
    GMX<-GM2/abs(GM2)
    GMX[is.nan(GMX)]<-1
    GM2<-abs(GM2)
    #GM2<-log1p(GM2)
    GM2<-GM2/max(GM2)
    GM<-GM2*GMX
    
    corrA<-corrplot(GM, tl.col = tl.col, 
                    tl.srt = tl.srt, is.corr = TRUE,
                    col = colMatrix, 
                    tl.cex = tl.cex,
                    pch.col=pch.col, cl.lim =  cl.lim)
    
  }
  
  if (graph.type=="clusplot" | graph.type=="all" ){
    
    dimMat<-dim(GM)
    if (dimMat[2]>=dimMat[1]){
      clusMat<-t(GM)
      colX<-which(colSums(clusMat)==0 & colMeans(clusMat)==0)
      rowX<-which(rowSums(clusMat)==0 & rowMeans(clusMat)==0)
      if (!is.integer0(colX)) { clusMat<-clusMat[ , -colX ] }
      if (!is.integer0(rowX)) { clusMat<-clusMat[ -rowX , ] }
      
      clusplot(pam(clusMat, nc), color = color, shade = shade,cex=cex,
               labels = labels, lines = lines,main=paste0(main," method: PAM n.cluster = ",nc))
    }else{
      warning("is impossible to calculate a clusplot using a matrix with more rows than column")
    }
  }
  
  if (graph.type=="pvclust" | graph.type=="all" ){
    
    plot(mPT@matrix$GFit, main = paste0(main," method: PAM n.cluster = ",nc))
    pvrect(mPT@matrix$GFit, alpha=alpha, lwd=lwd ,pv=pv, border = border)
    
  }
  
}
