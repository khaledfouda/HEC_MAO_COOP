AIS_fit <- function(Y, lambda, max.rank = 10, maxit=100, tol=1e-6, thresh=1e-3, trace.it=FALSE){
   
   yfill <- Y
   ymiss <- Y == 0
   m <- dim(Y)[1]
   n <- dim(Y)[2]
   R <- matrix(rnorm(n), ncol = 1)
   
   Q <- Compute_basis(Y, R, maxit=5, tol=tol)
   V0 <- V1 <- svd(t(Q) %*% Y, nu=0)$v
   
   iter <- 0
   ratio <- 1
   
   while((ratio > thresh)&(iter < maxit)){
      iter <- iter + 1
      lambdai <- abs(max.rank - lambda) * (decay ^ (i - 1)) + lambda

      
   }
   
}