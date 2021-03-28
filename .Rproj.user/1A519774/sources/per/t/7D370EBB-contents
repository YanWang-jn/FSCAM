#' FSCAM: CAM-based feature selection for clustering scRNA-seq
#'
#' @param procdata The preprocessed single cell data, which is a numeric matrix or data.frame. Rows represent genes/features cluster centers and columns represent single cells.
#' @param K The number of cell clusters. If NULL , then it will be estimated by function Cluestimate.
#'
#' @return The name vector of the selected DE genes.
#' @export
#'
#' @examples data(lpsdata)
#' procdata = preprocess(lpsdata, takenormalize=TRUE, clusternum=350)
#' DEGs =FSCAM(procdata,K)
FSCAM <- function(procdata, K) {
  gene_cluster=as.matrix(procdata[[2]])
  data=as.matrix(procdata[[1]])
  genename=rownames(gene_cluster)
  if(is.null(K)){
    K= Cluestimate(data)
  }
  landmark=kmeans(t(data),K)
  landmark_cluster=as.matrix(landmark$cluster)
  templandmark=aggregate(as.matrix(t(data),nrow=dim(procdata)[2],ncol=dim(data)[1]) ,by=list(landmark_cluster),FUN=mean)
  templandmark=t(templandmark[,-1])
  marker_gene=matrix(numeric(0), 0,0)
  for(j in 1:ceiling(K/2)){
    sample1=1:K
    sample1=sample(sample1,2)
    templandmark1=templandmark[,sample1]
    J <- dim(templandmark1)[1]
    cluster1=as.matrix(templandmark1)
    convex <- convhulln(rbind(templandmark1,0))
    corner <- matrix(numeric(0), 0,0)
    for (i in 1:2){
      corner <- union(corner,convex[,i])
    }

    for (j in 1:length(corner)){
      if (corner[j]==(J+1)){
        break
      }
    }
    corner <- corner[-j]  # throw away the origin point
    J_out <- length(corner)

    ##### estimate A and S ###########
    cat("\nEstimating A and S ... \n")
    cornerResult <- FOA_measure_convex(t(templandmark1[corner,]),2)
    ind <- cornerResult[[1]]
    measure_corner=corner[ind==1]
    for(i in 1:dim(matrix(measure_corner))[1])
    {
      marker_gene=c(marker_gene,genename[which(gene_cluster==measure_corner[i])])
    }
  }
  return(marker_gene)
}
