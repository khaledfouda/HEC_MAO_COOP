# This functions generates either model M1 or M2 from Yousef's simulation. Set model=1 for M1 or model=2 for M2.
generate_simulation_data_ysf <- function(model=1, n1 =300,  n2 = 300, m1 = 5, m2 = 10, missing_prob = 0.7,coll=FALSE,
                                         seed=2023){
   set.seed(seed=seed)
   X <- matrix(rnorm(n1*m1), ncol = m1)
   Z <- matrix(rnorm(n2*m2), ncol = m2)
   E <- matrix(rnorm(n1*n2), ncol = n2)
   # normalize X
   #X <- normalize_matrix(X)
   beta.x <- matrix(runif(m1*n2), ncol=n2)
   beta.z <- matrix(runif(m2*n1), ncol=n1)
   B.x <- matrix(runif(n1*m2), ncol=m2)
   B.z <- matrix(runif(m2*n2), ncol=n2)
   B <- B.x %*% B.z
   
   # if collinearity is needed, make the correlation between the first two columns in X and Z between (.99,1)
   # the actual correlation value is very close to 95%
   if(coll == TRUE){
      X[,2]  <- X[,1] + rnorm(n1, mean = 0, sd = 0.001)
      Z[,2]  <- Z[,1] + rnorm(n2, mean = 0, sd = 0.001)
   }
   
   
   P_X = X %*% solve(t(X) %*% X) %*% t(X)
   P_bar_X = diag(1,n1) - P_X
   P_Z = Z %*% solve(t(Z) %*% Z) %*% t(Z)
   P_bar_Z = diag(1,n2) - P_Z
   
   W <- matrix( rbinom(n1*n2, 1, (1 - missing_prob) ) , nrow = n1)
   
   if(model == 1){
      A <- (X %*% beta.x) + t(Z %*% beta.z) + P_bar_X %*% B %*% P_bar_Z 
      Y <- (A+E) * W
      rank <- qr(A)$rank
      return(list(A=A, W=W, X=X, Z=Z, Y=Y, beta.x=beta.x, beta.z=beta.z, B=B, rank=rank))
   }else if (model == 2){
      A <- (X %*% beta.x) + P_bar_X %*% B 
      Y <- (A + E) * W
      rank <- qr(A)$rank
      return(list(A=A, W=W, X=X, Y=Y, beta.x=beta.x, B=B, rank=rank))
   }
   else{
      stop("Error: Unrecognized model.")
   }
   #-----------------------------------------------------------------------------
   #---------------------------------------------------------------------
   
}
