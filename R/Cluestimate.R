
#' 	Estimate the number of cell clusters in a single-cell dataset
#'
#' @param procdata The preprocessed single cell data, which is a numeric matrix or data.frame. Rows represent genes/features cluster and columns represent single cells.
#' @param K_max Maximum number of clusters. If NULL , then it will be set as 20.
#'
#' @return The optimal number of clustering clusters.
#' @export
#'
#' @examples data(lpsdata)
#' procdata = preprocess(lpsdata, takenormalize=TRUE, clusternum=350)
#' K_opt = Cluestimate(procdata, K_max=20)
Cluestimate <- function(procdata, K_max=20){
  data=as.matrix(procdata[[1]])
  PC_num=fa.parallel(data,fa='pc',n.iter = 100,show.legend = F,main = 'Scree plot with parallel analysis')
  if (dim(data)[1]>dim(data)[2]){
    car.pr1= principal(data,nfactors = PC_num[["ncomp"]],rotate = 'varimax')
    b=round(unclass(car.pr1$weights),PC_num[["ncomp"]])
  }else{
    car.pr1= princomp(data)
    b=car.pr1$loadings[,1:PC_num[["ncomp"]]]
  }
  K_opt=fviz_nbclust(b, kmeans,k.max = K_max, method = "silhouette")
  K_opt=which(K_opt[["data"]][["y"]]==max(K_opt[["data"]][["y"]]))
  return(K_opt)
}
