#' landscapePlot
#'
#' draw a plot for multiple associations starting from a matrix of local z-scare multi associations 
#'
#' @usage landscapePlot<-function(lzmat,rotation=1,
#'          palette_vec=c("brown","darkgreen","lightgreen","yellow","lightblue"),
#'          theme_vec=c("lightblue","black","lightblue","white"))
#' 
#'
#' @param lzmat matrix of mutiple local z-scare associations
#'
#' @param rotation sliding of the central axis (in base pairs)
#' @param palette_vec (default = c("brown","darkgreen","lightgreen","yellow","lightblue"))is a colorRampPaette or the color of the landscape
#' 
#' @param theme_vec (default = c("lightblue","black","lightblue","white"))palette for the background elements of the plot
#' 
#'
#' @seealso \code{\link{localZScore}} \code{\link{regioneR}}
#'
#' @examples  ...
#'
#' @export lZAssociations
#'



landscapePlot<-function(lzmat,rotation=1,
                        palette_vec=c("brown","darkgreen","lightgreen","yellow","lightblue"),
                        theme_vec=c("lightblue","black","lightblue","white")){  
  
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