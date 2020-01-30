#' Plot Matrix
#'
#' Create a image of association matrix
#'
#' @usage plotMatrix( mat, palette = "bwr", type = "normal", name = "", ...)
#'
#' @param mat association matrix
#'
#' @param palette c( "bwr", "gbr") association's color palette blue-white-red or
#' green--black-red
#' @param type c( "normal", "compare")
#' @param name string name of the plot
#' @param ...  further arguments to be passed to other methods.
#'
#'
#' @seealso ...
#'
#' @examples  ...
#'
#' #@references ...
#' #@importFrom ...
#' @export plotMatrix
#'
plotMatrix<-function(mat,palette="bwr",type="normal",cexX=1,cexY=1,name="",...){
  atx<-seq(-1.1,1.1,length.out=100)
  if (palette=="bwr"){rgb.palette <- colorRampPalette(c("blue","white","red"),
                                                      space = "rgb")
  }
  if (palette=="gbr"){rgb.palette <- colorRampPalette(c("green","black","red"),
                                                      space = "rgb")
  }
  trellis.par.set(regions=list(col=rgb.palette(100)))
  if (type=="normal"){
    mat[mat>1]<-0.99
    mat[mat<(-1)]<-(-0.99)
    pl<-levelplot(mat,at=atx, scales=list(x=list(rot=45,cex=cexX),y=list(cex=cexY)),main=name)
  }
  if (type=="compare"){
    pl<-levelplot(mat, scales=list(x=list(rot=45,cex=cexX),y=list(cex=cexY)),main=name)
  }
  return(pl)
}

