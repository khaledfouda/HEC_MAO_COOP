compute_test_estimates <- function(test_indices, Y.minus.B, X, beta.estim){
   estimates <- numeric(length(test_indices[, 1]))
   for (idx in seq_along(estimates)) {
      i <- test_indices[idx, 1]
      j <- test_indices[idx, 2]
      estimates[idx] <- Y.minus.B[i, j] + sum(X[i, ] * beta.estim[, j])
   }
   return(estimates)
}


optimize_lambda1 <- function(W, Y.minus.B, X, lambda1.grid, Y.valid, trace=FALSE, acc=4){
   #' W: Validation mask where Wij=0 if the cell belongs to the validation set and 1 if belongs to the training or test
   #' fiti: training object containing the SVD decomposition of the low-rank matrix (UDV) 
   #'  lambda1.grid: grid of possible values for lambda1
   #'  Y.valid: A vector of the true values of Y for the validation set.
   #'  ------------------
   #'  Returns: best.lambda1: optimal value
   #'           Xbeta: the term X^T beta corresponding to the optimal lambda
   X.X = t(X) %*% X
   XY <- t(X) %*% Y.minus.B
   validation_indices <- which(W == 0, arr.ind = TRUE)
   best_score = Inf
   best_lambda1 = NULL
   best_beta = NULL
   nr <- nrow(X.X)
   n1n2 = nrow(W)*nrow(W)#nr*nr#svd(X.X)$d[1]
   
   for(lambda1 in lambda1.grid){
      beta_partial = solve(X.X + n1n2*lambda1* diag(1, nr))
      beta.estim = beta_partial %*% XY
      soft_estim = Y.minus.B  + X %*% beta.estim
      #soft_estim = compute_test_estimates(validation_indices, Y.minus.B, X, beta.estim)
      err = round(test_error(soft_estim[W==0], Y.valid),acc)
      if(err < best_score){
         best_score = err
         best_lambda1 = lambda1
         beta_partial = beta_partial
         best_beta = beta.estim
         #print(paste(err, lambda1))
      }
   }
   list(lambda1 = best_lambda1, score=best_score, beta_partial =beta_partial%*% t(X), beta=best_beta)
}


simpute.cov.cv.v2 <- function(Y, X, W, A, 
                              # lambda 2 for the regularization on the low-rank matrices A & B
                              lambda2.factor=1/4, lambda2.init=NA, n.lambda2=20,
                              # lambda 1 for the regularization on the covariate coefficients
                              lambda1.optimize = FALSE, lambda1.grid = seq(0,400,length.out=20),
                              # trace prints at every iteration and print.best prints at the end.
                              trace=FALSE, print.best=TRUE, 
                              # threshold and tolerance for convergence
                              thresh=1e-5,  tol=5,
                              # rank limit parameters. minimum rank is "rank.init" and maximum is "rank.limit".
                              rank.init=10, rank.limit=50, rank.step=2,
                              # type of fitting function: "als" or "svd" only. (considering adding original?)
                              type="als"){
   
   # choose fitting function
   stopifnot(type %in% c("svd", "als"))
   if(type == "svd"){
      fit.function <- simpute.svd.cov
   }else
      fit.function <- simpute.als.cov
   
   # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
   #Y[Y==0] = NA
   
   # initial lambda is selected here as well as the sequence of lambdas to fit.
   lam0 <- ifelse(is.na(lambda2.init), lambda0.cov(Y, X) * lambda2.factor, lambda2.init) 
   lamseq <- seq(from=lam0, to=0, length=n.lambda2)
   
   # variable initialization
   fits <- as.list(lamseq)
   ranks <- as.integer(lamseq)
   rank.max <- rank.init
   warm <- NULL
   best_estimates = NA
   best_fit <- list(error=Inf, rank_A=NA, rank_B=NA, lambda2=NA, rank.max=NA)
   counter <- 1
   Y.minus.B = Y
   Y.valid = A[W==0]
   X.X = t(X) %*% X
   beta_partial = solve(X.X) %*% t(X)
   best_lambda1 = 0
   #--------------------------------------------------
   # we start by fidning the optimal lambda 1 
   
   #-----------------------------------------------------
   # for each lambda in the grid do ...>
   for(i in seq(along=lamseq)) {
      # find the optimal lambda before going to every fit function
      #fiti = optimize_lambda1(W, Y.minus.B, X, lambda1.grid, Y.valid, trace=FALSE)
      #best_lambda1 = fiti$lambda1
      #beta_partial = fiti$beta_partial
      #-----------------------------------------------------
      fiti <- fit.function(Y, X, beta_partial, thresh=thresh, lambda = lamseq[i], J=rank.max, warm.start = warm)
      
      # compute rank.max for next iteration
      rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
      rank.max <- min(rank+rank.step, rank.limit)
      
      # get test estimates and test error
      v=as.matrix(fiti$v)
      vd=v*outer(rep(1,nrow(v)),fiti$d)
      B_hat = fiti$u %*% t(vd)
      soft_estim = B_hat  + X %*% fiti$beta.estim
      Y.minus.B = Y - B_hat
      err = test_error(soft_estim[W==0], A[W==0])
      #----------------------------
      warm <- fiti # warm start for next 
      if(trace==TRUE)
         cat(sprintf("%2d lambda1=%9.5g, lambda2=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
                     i, best_lambda1, lamseq[i], rank.max, rank, err))
      #-------------------------
      # register best fir
      if(err < best_fit$error){
         best_fit$error = err
         best_fit$rank_A = qr(soft_estim)$rank
         best_fit$rank_B = rank
         best_fit$lambda1 = best_lambda1
         best_fit$lambda2 = lamseq[i]
         best_fit$rank.max = rank.max
         best_estimates = soft_estim
         best_beta = fiti$beta.estim
         best_B = fiti$u %*% t(vd)
         counter=1
      }else counter = counter + 1
      if(counter >= tol){
         cat(sprintf("Performance didn't improve for the last %d iterations.", counter))
         break
      }
   }
   #-----------
   # step 2: Optimizing for lambda.1 after finding the optimal lambda.2
   fiti = optimize_lambda1(W, Y.minus.B, X, lambda1.grid, Y.valid, trace=FALSE)
   best_fit$lambda1 = fiti$lambda1
   best_beta = fiti$beta
   #----------------------------------------------------------------------
   
   if(print.best==TRUE) print(best_fit)
   best_fit$B_hat = best_B
   best_fit$beta_hat = best_beta
   best_fit$A_hat = best_B + X %*% best_beta
   best_fit$rank_A = qr(best_fit$A_hat)$rank
   return(best_fit)
}
#--------------------------------------------------------------------
