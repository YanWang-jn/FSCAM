#' Single-cell clustering algorithm based on FSCAM method
#'
#' @param procdata The preprocessed single cell data, which is a numeric matrix or data.frame. Rows represent genes/features cluster and columns represent single cells.
#' @param K The number of cell clusters. If NULL , then it will be estimated by function Cluestimate.
#' @param T Running times of FSCAM feature selection algorithm. If NULL , then it will be set as 20.
#'
#' @return The name vector of the selected features and K types of cell-labels for cells in input dataset.
#' @export
#'
#' @examples data(lpsdata)
#' procdata = preprocess(lpsdata, takenormalize=TRUE, clusternum=350)
#' DEGs = FSCAM(procdata,K)
#' features_labels = SCC_FSCAM(procdata,K=NULL, T=20)
SCC_FSCAM <- function(procdata,K=NULL, T=20) {
  procdata_c=list(procdata[[1]],procdata[[2]])
  procdata_uc=as.matrix(procdata[[3]])
  if (is.null(K)){
    K=Cluestimate(procdata_c)
  }
  rescult=matrix(numeric(0), 0,0)
  for (i in 1:T){
    try({
      rescult=c(rescult,FSCAM(procdata_c,K))
    })
  }
  precise_g=matrix(numeric(0), 0,0)
  for (k in 1:(length(rescult))){
    precise_g=c(precise_g,rescult[[k]])
  }
  precise_g_score=table(precise_g)
  precise_g_score=sort(precise_g_score,decreasing = TRUE)
  precise_select_g=precise_g_score[(precise_g_score!=1)]
  a=data.frame(procdata_uc[names(precise_select_g),])
  d = as.dist(1-cor(a,method="spearman"))
  kp<- pamk(d,krang=K,diss = TRUE)
  feature_set=names(precise_select_g)
  cell_label=kp$clustering
  return(list(feature_set,cell_label))
}
