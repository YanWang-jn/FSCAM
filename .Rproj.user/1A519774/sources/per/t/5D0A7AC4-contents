# This software is under license BSD 3.
# This function implements the Lawson and Hanson method for solving the least
# squares problems with nonnegativity constraints.
#
# Usage: x <- nnls(C,d)
# returns the vector x that minimizes Norm(d-C*x)
# subject to x >= 0;


#' Title
#'
#' @param C The independent variables
#' @param d The dependent variable
#'
#' @return the Regression coefficient vector x that minimizes Norm(d-C*x)
#' @export
#'
#' @examples x=nnls(C,d)
nnls2 <- function(C,d) {
# Set initial parameters.
# ----------------------------------------------------------------
  M   <- nrow(C)
  N   <- ncol(C)
  tol <- 10 * (.Machine$double.eps) * norm(C,'1') * max(M,N)
  x0  <- 1:N
  x1  <- rep(0,N)

  x0_pos<- x0
  x     <- x1

  re <- d - C %*% x
  w     <- t(C) %*% re

  iter  <- 0
  itmax <- 3 * N
# ----------------------------------------------------------------
  CP <- matrix(0,M,N)
  while (any(x0 != 0) && any(w[x0_pos] > tol)) {
    t <- which.max(w[x0_pos])
    t <- x0_pos[t]

    x0[t] <- 0
    x1[t] <- t
    x0_pos   <- which(x0 != 0)
    x1_pos   <- which(x1 != 0)

    Ns  	<- length(x0_pos)
    CP[,x0_pos] <- matrix(0,M,Ns)
    CP[,x1_pos] <- C[,x1_pos]
    x_iter      <- ginv(CP) %*% d
    x_iter[x0_pos] <- rep(0,Ns)

    while (any(x_iter[x1_pos] <= tol)) {
      iter <- iter + 1
      if (iter > itmax) {
        x = x_iter
        return(x)
      }

      QQ    <- which((x_iter <= tol) & (x1 != 0))
      alpha <- min(x[QQ] / (x[QQ] - x_iter[QQ]))
      x     <- x + alpha * (x_iter - x)
      index    <- which((abs(x) < tol) & (x1 != 0))

      x0[index] <- index
      x1[index] <- rep(0,length(index))
      x1_pos    <- which(x1 != 0)
      x0_pos    <- which(x0 != 0)

      Ns   	  <- length(x0_pos)
      CP[,x0_pos] <- matrix(0,M,Ns)
      CP[,x1_pos] <- C[,x1_pos]
      x_iter      <- ginv(CP) %*% d
      x_iter[x0_pos] <- rep(0,Ns)
    }
    x     <- x_iter
    re <- d - C %*% x
    w     <- t(C) %*% re
  }

  return(x)
}
