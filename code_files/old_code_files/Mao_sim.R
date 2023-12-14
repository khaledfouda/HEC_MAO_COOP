normalize_matrix <- function(X) { # NOT USED FOR NOW
   # Normalize rows
   X_row_normalized <- t(apply(X, 1, function(row) {
      (row - mean(row)) / sd(row)
   }))
   
   # Normalize columns
   X_col_normalized <- apply(X_row_normalized, 2, function(col) {
      (col - mean(col)) / sd(col)
   })
   
   return(X_col_normalized)
}


# This function generates Mao's simulation as defined in (M0) above. Dimensions, missing probability are given as parameters. 
generate_simulation_data_mao <- function(n1 =400,  n2 = 400, m = 20, r = 10, MAR=TRUE, seed=2023){
   #' Input: 
   #'      n1, n2: are the dimensions of the A, Y, and B matrices
   #'      m: number of covariates
   #'      r: Second dimension of U and V where  B = P_bar_X U V^T
   #'      MAR: If True, missing at random method is employed to compute the missing probability and W
   #'      seed: random seed  
   set.seed(seed=seed)
   X <- matrix(rnorm(n1*m), ncol = m) #%>% normalize_matrix()
   beta <- matrix(rnorm(m*n2), ncol=n2)
   U <- matrix(rnorm(n1*r),ncol=r)
   V <- matrix(rnorm(n2*r),ncol=r)
   P_X = X %*% solve(t(X) %*% X) %*% t(X)
   P_bar_X = diag(1,n1) - P_X 
   B = P_bar_X %*% U %*% t(V)
   A <- X %*% beta + B
   rank <- qr(A)$rank # = m + r
   #-----------------------------------------------------------------------------
   # Simulation Theta using missing-at-random MAR method
   gamma <- rbind( matrix(rnorm(1, -1.5, 0.1),1), matrix(rnorm(3, 0.3, 0.1),3))
   # only the first 3 columns of X are used
   # theta is the inverse of the probability of missingness following a logistic model 
   inclusion_prob = 1 / (1 + exp( -(cbind(1,X[,1:3]) %*% gamma)) ) #  all values should be close to 0.2
   theta = 1 / inclusion_prob
   # Wij=1 if A is observed following the probabilities of missingness defined above.
   W <- matrix( rbinom(n1*n2, 1, inclusion_prob ) , nrow = n1)
   # sum(W==0)/(n1*n2) should be equal to 0.2
   #----------------------------------------------------------
   # Does fully observed Y = A (ie,  without noise?)? In that case ignore the code below.
   #----------------------------------------------------------------------
   # Computing epsilon as iid zero mean Gaussian with variance chosen such that the signal-to-noise ratio (SNR) is 1
   signal_A <- sum((A - mean(A))^2) / (n1 * n2 - 1)
   sigma_epsilon <- sqrt(signal_A)  # Since SNR = 1
   epsilon <- matrix(rnorm(n1 * n2, mean = 0, sd = sigma_epsilon), n1, n2)
   #--------------------------------------------------------------------------------------------
   # Y is a corrupted, partially-observed version of A. Yij=0 for missing data [maybe consider making it NA]
   Y <- (A + epsilon) * W
   #---------------------------------------------------------------------
   return(list(A=A, W=W, X=X, Y=Y, beta=beta, B=B, theta=theta, gamma=gamma, inclusion_prob=inclusion_prob, rank=rank))
}
