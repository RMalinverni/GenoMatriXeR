#' Redux and Compare Matrix
#'
#' compare 2 matrices and redux it choosing only the shared elements
#'
#' @usage reduxAndCompareMatrix(mat.A,mat.B,nameA="",nameB="",mirrored=FALSE,
#'method="complete",...)
#'
#' @param mat.A association matrix
#' @param mat.B association matrix
#' @param nameA string name for the plot related to matrix A
#' @param nameB string name for the plot related to matrix B
#' @param mirrored boolean choose if the matrix will be mirrored
#' @param method for the (\code{\link{hclust}}
#'
#'
#' @seealso ...
#'
#' @examples  ...
#'
#' #@references ...
#' #@importFrom ...
#' @export reduxAndCompareMatrix
#'
reduxAndCompareMatrix<-function(mat.A,mat.B,nameA="",nameB="",mirrored=FALSE,
                                method="complete",...){

  colnames(mat.A)<-names.A<-sub(".*-", "", colnames(mat.A))
  colnames(mat.B)<-names.B<-sub(".*-", "", colnames(mat.B))
  rownames(mat.A)<-sub(".*-", "",rownames(mat.A))
  rownames(mat.B)<-sub(".*-", "",rownames(mat.B))

  mat1<-mat.A[which(rownames(mat.A)%in%rownames(mat.B)),
              which(colnames(mat.A)%in%colnames(mat.B))]
  mat1<-mat1[!duplicated(rownames(mat1)),!duplicated(colnames(mat1))]
  mat2<-mat.B[which(rownames(mat.B)%in%rownames(mat.A)),
              which(colnames(mat.B)%in%colnames(mat.A))]
  mat2<-mat2[!duplicated(rownames(mat2)),!duplicated(colnames(mat2))]
  mat1<-mat1[order(rownames(mat1)),order(colnames(mat1))]
  mat2<-mat2[order(rownames(mat2)),order(colnames(mat2))]


  mat1<-clustMatrix(mat1,method=method,mirrored=mirrored)
  mat2<-mat2[match(rownames(mat1),rownames(mat2)),match(colnames(mat1),
                                                        colnames(mat2))]
  mat3<-mat1-mat2


  mat.tot<-list(mat1,mat2,mat3)
  names(mat.tot)<-c(nameA,nameB,paste0(nameA,"-",nameB))

  return(mat.tot)

}
