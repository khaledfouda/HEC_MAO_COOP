#' The following function, $Mao.fit()$, estimates $\hat\beta$ and $\hat B$ as above with fixed (given) hyperparameters.
#'  I will try to use the same notations as above to avoid confusion.


# EDIT: These functions return the 1 / (probability of inclusion) NOT missingness.
theta_default <- function(X, W, ...){
   # using logistic regression as indicated in (a1)
   n1 = dim(W)[1]
   n2 = dim(W)[2]
   theta_hat = matrix(NA, n1, n2)
   for(j in 1:n2){
      model_data = data.frame(cbind(W[,j], X))
      model_fit = glm(X1~., family=binomial(), data=model_data)
      theta_hat[, j] = 1 / predict(model_fit, type="response")
   }
   return(theta_hat)
}
theta_random <- function(W, ...){
   # A theta estimation function that selects the proportion of missing data within the same column
   # using formula (a2)
   n1 = dim(W)[1]
   n2 = dim(W)[2]
   theta_hat = matrix(NA, n1, n2)
   for(j in 1:n2){
      theta_hat[,j] = n1 / sum(W[,j]==1) 
   }
   return(theta_hat)
}
theta_simple <- function(W, ...){
   # A theta estimation function that selects the proportion of missing data in the matrix
   # using formula (a3)
   n1 = dim(W)[1]
   n2 = dim(W)[2]
   theta_hat = matrix( (n1*n2) / sum(W==1), n1, n2)
   return(theta_hat)
}


Mao.fit <- function(Y, X, W, lambda.1, lambda.2, alpha, n1n2_optimized=TRUE, theta_estimator=theta_default){
   #
   #' ----------------------------------------------
   #' Input: Y: corrupted, partially observed A, (Y is assumed to be the product of Y*W)
   #'        W: Wij=1 if Aij is observed (ie, Yij !=0), and 0 otherwise
   #'         X: covariate matrix
   #'         lambda.1, lambda.2, alpha: hyperparameters
   #' ----------------------------------------------
   #' output: list of  A, Beta_hat, B_hat
   #' ----------------------------------------------
   n1 = dim(Y)[1]
   n2 = dim(Y)[2]
   m  = dim(X)[2]
   # The following two lines are as shown in (c) and (d)
   X.X = t(X) %*% X
   P_X = X %*% solve(X.X) %*% t(X)
   P_bar_X = diag(1,n1) - P_X 
   
   if(n1n2_optimized == TRUE){
      # we define the factor that will be used later:
      n1n2 = svd(X.X)$d[1] 
   }else{
      n1n2 = n1 * n2
   }
   
   # The following part estimates theta (missingness probabilities)
   theta_hat = theta_estimator(W=W, X=X)
   # the following is the product of W * theta_hat * Y
   W_theta_Y = Y * theta_hat # * W 
   
   # beta hat as (8)
   beta_hat = solve(X.X + n1n2 * lambda.1 * diag(1, m)) %*% t(X) %*% W_theta_Y
   # SVD decomposition to be used in (b)
   svdd = svd(P_bar_X %*% W_theta_Y)
   if(n1n2_optimized == TRUE){
      # evaluation of  (b)
      n1n2 = svdd$d[1]
   }else{
      n1n2 = n1 * n2
   }
   T_c_D = svdd$u %*% (pmax(svdd$d - alpha*n1n2*lambda.2, 0) * t(svdd$v))
   # B hat as in (11)
   B_hat = (1 + 2 * (1-alpha) * n1n2 * lambda.2 )^(-1) * T_c_D
   # computing the rank of B [Copied from Mao's code; Don't understand how it works.]
   # EQUIVALENT to  qr(B_hat)$rank + m   or   qr(A_hat)$rank
   # B is a low rank matrix
   rank = sum(pmax(svdd$d - n1n2 * lambda.2 * alpha, 0) > 0) + m
   
   # Estimate the matrix as given in the model at the top
   A_hat = X %*% beta_hat + B_hat
   
   return(list(A_hat = A_hat, B_hat = B_hat, beta_hat = beta_hat, rank = rank))
}
