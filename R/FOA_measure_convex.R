# This software is under license BSD 3.

#' Estimate vertices of convex sets
#'
#' @param X Convex set, which is a numeric matrix or data.frame. Rows represent genes/features cluster and columns represent single cells.
#' @param N The number of vertices
#'
#' @return Estimated vertices number
#' @export
#'
#' @examples Y=FOA_measure_convex(X,N)
FOA_measure_convex <- function(X,N){

  M <- nrow(X); L <- ncol(X)
  Npop=50
  r=0.3
  Nfeature=N
  error <- matrix(0,Npop,1)
  Nrange=10
  temperror <- matrix(0,Nrange,1)
  sample1=c(rep(1,Nfeature),rep(0,L-Nfeature))
  pop=matrix(nrow = L,ncol = Npop)
  for (i in 1:Npop) {
    pop[,i]=sample(sample1,length(sample1))
  }

  for (p in 1:Npop) {
    A <- X[,pop[,p]==1]
    Others <- X[,pop[,p]==0]
    Ae <- rbind(1e-5*A, matrix(1,1,sum(pop[,p]==1)))
    Oe <- rbind(1e-5*Others, matrix(1,1,sum(pop[,p]==0)))
    alpha <- matrix(0,sum(pop[,p]==1),sum(pop[,p]==0))
    for (i in 1:sum(pop[,p]==0)) {
      alpha[,i] <- nnls2(Ae, Oe[,i])
    }

    error[p] <- norm(Others - A %*% alpha, 'f')^2
  }


  val <- min(error)
  ind <- which.min(error)
  best_smell=val
  best_fly=pop[,ind]
  history_best_fly=best_fly
  history_best_smell=best_smell

  m=1
  t=0
  tmax=20
  while(t<=tmax){
    t=t+1
    for (j in 1:Npop) {
      temppop=rep(pop[,j],Nrange)
      temppop=matrix(temppop,nrow = L,ncol = Nrange)
      for (i in 1:Nrange) {
        temppop[sample(c(which(temppop[,i]==0)),1),i]=1
        temppop[sample(c(which(temppop[,i]==1)),1),i]=0
      }

      for (p in 1:Nrange) {
        A <- X[,temppop[,p]==1]
        Others <- X[,temppop[,p]==0]
        Ae <- rbind(1e-5*A, matrix(1,1,sum(temppop[,p]==1)))
        Oe <- rbind(1e-5*Others, matrix(1,1,sum(temppop[,p]==0)))
        alpha <- matrix(0,sum(temppop[,p]==1),sum(temppop[,p]==0))
        for (i in 1:sum(temppop[,p]==0)) {
          alpha[,i] <- nnls2(Ae, Oe[,i])
        }

        temperror[p] <- norm(Others - A %*% alpha, 'f')^2
      }

      tempval <- min(temperror)
      tempind <- which.min(temperror)
      pop[,j]=temppop[,tempind]
    }


    for (p in 1:Npop) {
      A <- X[,pop[,p]==1]
      Others <- X[,pop[,p]==0]
      Ae <- rbind(1e-5*A, matrix(1,1,sum(pop[,p]==1)))
      Oe <- rbind(1e-5*Others, matrix(1,1,sum(pop[,p]==0)))
      alpha <- matrix(0,sum(pop[,p]==1),sum(pop[,p]==0))
      for (i in 1:sum(pop[,p]==0)) {
        alpha[,i] <- nnls2(Ae, Oe[,i])
      }

      error[p] <- norm(Others - A %*% alpha, 'f')^2
    }


    val <- min(error)
    ind <- which.min(error)


    if(history_best_smell[m]>val){
      m=m+1
      best_smell=val
      best_fly=pop[,ind]
      history_best_fly=cbind(history_best_fly,best_fly)
      history_best_smell=cbind(history_best_smell,best_smell)
      for (j in 1:Npop){
        simulate_set1=sample(which(best_fly==1),ceiling(N*r))
        increase1=dim(matrix(which(pop[simulate_set1,j]==0)))[1]
        simulate_set0=sample(setdiff(which(pop[,j]==1),simulate_set1),increase1)
        pop[simulate_set0,j]=0
        pop[simulate_set1,j]=1
      }
    }else{
      best_smell=val
      best_fly=pop[,ind]
      if(m==1){
        history_best_fly=cbind(history_best_fly,history_best_fly)
      }else{
        history_best_fly=cbind(history_best_fly,history_best_fly[,m])
      }
      history_best_smell=cbind(history_best_smell,history_best_smell[m])
      m=m+1
      for (j in 1:Npop){
        simulate_set1=sample(which(best_fly==1),ceiling(N*r))
        decrease1=dim(matrix(which(pop[simulate_set1,j]==0)))[1]
        simulate_set0=sample(setdiff(which(pop[,j]==1),simulate_set1),decrease1)
        pop[simulate_set0,j]=0
        pop[simulate_set1,j]=1
      }
    }

  }
  cornerind <- history_best_fly[,m]
  return(list(cornerind))

}
