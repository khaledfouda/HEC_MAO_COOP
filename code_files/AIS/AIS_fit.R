AIS_fit <- function(Y, lambda, maxR = 10, maxit=100, tol=1e-6, thresh=1e-3,decay=0.8, trace.it=FALSE){
   
   yfill <- Y
   ymiss <- Y == 0
   m <- dim(Y)[1]
   n <- dim(Y)[2]
   R <- matrix(rnorm(n), ncol = 1)
   
   lambdaMax <- svd(D, nu = 1, nv = 0)$d[1]
   Q <- Compute_basis(Y, R, maxit=5, tol=tol)
   V0 <- V1 <- svd(t(Q) %*% Y, nu=0)$v
   
   iter <- 0
   ratio <- 1
   
   while((ratio > thresh)&(iter < maxit)){
      #svtY_old = svtY
      iter <- iter + 1
      lambdai <- abs(lambdaMax - lambda) * (decay ^ (iter - 1)) + lambda
      theta = (iter-1)/(iter+2)
      Z = M + theta * (M - M.old)
      yfill[ymiss] = M[ymiss]
      #-----------------------------------------------------------------
      # Compute R using V1 and V0
      R <- V0 - V1 %*% t(V1) %*% V0
      R <- colSums(R^2) > 0
      R <- cbind(V1, V0[,R])
      R <- R[, 1:min(ncol(R), maxR)]
      #---------------------------------------
      svtY <- SVT_Approx(yfill, R, lambdai, maxit)
      M.old <- M
      M <- svtY$u %*% (svtY$d * t(svt$v))
      #------------------------------------------------------------
      # Update V1 and V0
      V0 <- V1
      V1 = svtY$v
      #-------------------
      # find another way for the ratio as U and U.old won't have the same dimension
      ratio=Frob(svtY_old$u,svtY_old$d,svtY_old$v,svtY$u,svtY$d,svtY$v)
      
   }
   
}