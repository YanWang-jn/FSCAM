#' preprocess
#'
#' @param data The raw single_cell data, which is a numeric matrix or data.frame. Rows represent genes/features and columns represent single cells.
#' @param takenormalize Logical value indicating whether to take normalize. If FALSE no normalizing will be performed.
#' @param clusternum The number of clusters for doing cluster, typically 5% of number of all genes. The clustering will be done after all the gene filtering and normalizing. If NULL no clustering will be performed.
#'
#' @return Preprocessed cluster-level expression matrix and cluster to which each gene belongs
#' @export
#'
#' @examples data(lpsdata)
#' procdata = preprocess(lpsdata, takenormalize=TRUE, clusternum=350)
preprocess <- function(data, takenormalize=TRUE, clusternum=350) {
  data=data[apply(data,1,sum)!=0,]
  gene_quantile=quantile(apply(data,1,sum),probs=c(0.05,0.95))
  data=data[apply(data,1,sum)>gene_quantile[1],]
  data=data[apply(data,1,sum)<gene_quantile[2],]
  Nsample=dim(data)[2]
  if (takenormalize==TRUE){
    data_max=apply(data,1,max)
    data_min=apply(data,1,min)
    data1=(data-data_min)/(data_max-data_min)
  }else{
    data1=data
  }
  Ngene=dim(data1)[1]
  genename=rownames(data1)
  cluster <- kmeans(data1,clusternum,iter.max=100)
  for (i in 1:50){
    tmp <- kmeans(data1,clusternum,iter.max=100)
    if (cluster$tot.withinss>tmp$tot.withinss){
      cluster <- tmp
    }
  }
  small_cluster <- matrix(numeric(0),0,0)
  for (k in 1:clusternum){
    if (cluster$size[k]<0.1*Ngene/clusternum)
      small_cluster <- c(small_cluster,k)
  }
  cluster1=cluster
  genename=rownames(data1)
  if (length(small_cluster)==0){
    cluster1 <- cluster1
  } else {
    cat("\nfind small cluster !!!... \n")
    for(i in 1:length(small_cluster)){
      genename=genename[-which(cluster1$cluster==small_cluster[i])]
      cluster1$cluster=cluster1$cluster[-which(cluster1$cluster==small_cluster[i])]
      cluster1$cluster[which(cluster1$cluster>small_cluster[i])]=cluster1$cluster[which(cluster1$cluster>small_cluster[i])]-1
      small_cluster=small_cluster-1
    }
  }
  gene_cluster=cluster1$cluster
  data2=data1[genename,]
  data3=aggregate(as.matrix(data2,nrow=dim(data2)[1],ncol=dim(data2)[2]) ,by=list(gene_cluster),FUN=mean)
  data3=data3[,2:dim(data3)[2]]
  return(list(data3,gene_cluster,data2))
}
