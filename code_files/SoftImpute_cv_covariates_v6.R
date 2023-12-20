
simpute.cov.kfold.lambda1 <- function(Y, X, W, lambda2=NA, J=NA,
                              trace=FALSE, print.best=TRUE, thresh=1e-5, n_folds=3,
                           type="als", lambda1.grid=NA, n1n2=1, warm=NULL){
   
   stopifnot(type %in% c("svd", "als"))
   if(type == "svd"){
      fit.function <- simpute.svd.cov
   }else
      fit.function <- simpute.als.cov
   
   n <- ncol(X)
   stopifnot(n1n2 %in% 1:3)
   # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
   
   #ymiss <- W == 0
   #rank.max <- rank.init
   #counter <- 1
   X.X = t(X) %*% X
   if(n1n2 == 2){
      n1n2 = svd(X)$d[1]
   }else if(n1n2 == 3){
      n1n2 = nrow(Y) * ncol(Y)
   }
   
   
   best_lambda1 = NA
   best_error = Inf
   
   fiti = warm
   
   folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, W)
   fold_data <- lapply(1:n_folds, function(i) {
      W_fold = folds[[i]]
      ymiss = W_fold==0 & W==1
      list(ymiss=ymiss, W_fold=W_fold, Y_train = Y * W_fold, Y_valid = Y[ymiss])
   })
   
   for(i in 1:length(lambda1.grid)){
      #print('hi')
      beta_partial = solve(X.X +  diag(n1n2*lambda1.grid[i], n)) %*% t(X)
      
      err = 0
      for(f in 1:n_folds){
         data = fold_data[[f]]
         fiti <- fit.function(data$Y_train, X, beta_partial, thresh=thresh, lambda = lambda2, J=J, warm.start = fiti)
         test_estim = (fiti$u %*%(fiti$d*t(fiti$v)))[data$ymiss] + (X %*% fiti$beta.estim)[data$ymiss]
         err = err + test_error(test_estim, data$Y_valid)
      }
      err = err / n_folds
      
      
      if(err < best_error){
         best_lambda1 = lambda1.grid[i]
         best_index = i
         best_error = err
      }
   if(trace==TRUE)
      cat(sprintf("%2d lambda1=%9.5g, lambda2=%9.5g, rank.max = %d, error = %.5f\n",
                  i, best_lambda1, lambda2, J, best_error))
   }
   
   beta_partial = solve(X.X +  diag(n1n2*lambda1.grid[best_index], n)) %*% t(X)
   fiti <- fit.function(Y, X, beta_partial, thresh=thresh, lambda = lambda2, J=J, warm.start = fiti)
   
   
   results = list()
   #results$fit = fiti
   results$lambda1 = best_lambda1
   results$lambda2 = lambda2
   results$B_hat = fiti$u %*%(fiti$d*t(fiti$v))
   results$A_hat = results$B_hat + X %*% fiti$beta.estim
   results$beta_hat = fiti$beta.estim
   results$rank_A = qr(results$A_hat)$rank
   results$J = J
   return(results)
}
#--------------------------------------------------------------------
simpute.cov.kfold <- function(Y, X, W, lambda.factor=1/4, lambda.init=NA, n.lambda=20,
                           trace=FALSE, print.best=TRUE, tol=5, thresh=1e-5, n_folds=3,
                           rank.init=3, rank.limit=50, rank.step=2,
                           type="als", lambda1=0, n1n2=1, warm=NULL){
   
   stopifnot(type %in% c("svd", "als"))
   if(type == "svd"){
      fit.function <- simpute.svd.cov
   }else
      fit.function <- simpute.als.cov
   
   stopifnot(n1n2 %in% 1:3)
   
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov(Y, X) * lambda.factor, lambda.init) 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   
   fits <- as.list(lamseq)
   ranks <- as.integer(lamseq)
   
   
   folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, W)
   fold_data <- lapply(1:n_folds, function(i) {
      W_fold = folds[[i]]
      ymiss = W_fold==0 & W==1
      list(ymiss=ymiss, W_fold=W_fold, Y_train = Y * W_fold, Y_valid = Y[ymiss])
   })
   
   rank.max <- rank.init
   #warm <- NULL
   best_estimates = NA
   best_fit <- list(error=Inf, rank_A=NA, rank_B=NA, lambda=NA, rank.max=NA)
   counter <- 1
   best_error <- Inf
   best_rank <- NA
   best_lambda <- NA
   X.X = t(X) %*% X
   if(n1n2 == 2){
      n1n2 = svd(X)$d[1]
   }else if(n1n2 == 3){
      n1n2 = nrow(Y) * ncol(Y)
   }
   beta_partial = solve(X.X +  diag(n1n2*lambda1, ncol(X))) %*% t(X)
   for(i in seq(along=lamseq)) {
   
      fiti <- fit.function(Y, X, beta_partial, thresh=thresh, lambda = lamseq[i], J=rank.max, warm.start = warm)
      
      err <- rank <-  0
      for(f in 1:n_folds){
         data = fold_data[[f]]
         fiti <- fit.function(data$Y_train, X, beta_partial, thresh=thresh, lambda = lamseq[i], J=rank.max, warm.start = fiti)
         test_estim = (fiti$u %*%(fiti$d*t(fiti$v)))[data$ymiss] + (X %*% fiti$beta.estim)[data$ymiss]
         err = err + test_error(test_estim, data$Y_valid)
         rank <- rank + sum(round(fiti$d, 4) > 0) # number of positive sing.values
      }
      err = err / n_folds
      rank = as.integer(rank / n_folds)
      
      # get test estimates and test error
      #v=as.matrix(fiti$v)
      #vd=v*outer(rep(1,nrow(v)),fiti$d)
      #soft_estim = fiti$u %*% t(vd)  + X %*% fiti$beta.estim
      #err = test_error(soft_estim[W==0], Y.valid)
      #----------------------------
      #warm <- fiti # warm start for next 
      if(trace==TRUE)
         cat(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
                     i, lamseq[i], rank.max, rank, err))
      #-------------------------
      # register best fir
      if(err < best_error){
         best_error = err
         best_lambda = lamseq[i]
         best_rank = rank.max
         counter=1
         #best_fit$error = err
         #best_fit$rank_A = qr(soft_estim)$rank
         #best_fit$rank_B = rank
         #best_fit$lambda = lamseq[i]
         #best_fit$rank.max = rank.max
         #best_estimates = soft_estim
         #best_beta = fiti$beta.estim
         #best_B = fiti$u %*% t(vd)
      }else counter = counter + 1
      if(counter >= tol){
         cat(sprintf("Performance didn't improve for the last %d iterations.", counter))
         break
      }
      # compute rank.max for next iteration
      rank.max <- min(rank+rank.step, rank.limit)
   }
   if(print.best==TRUE) cat(sprintf("lambda=%9.5g, rank.max = %d, error = %.5f\n",
                                    best_lambda, best_rank, best_error))
   
   fiti <- fit.function(Y, X, beta_partial, thresh=thresh, lambda = best_lambda, J=best_rank, warm.start = fiti)
   
   
   results = list()
   #results$fit = fiti
   results$lambda1 = lambda1
   results$lambda2 = best_lambda
   results$B_hat = fiti$u %*%(fiti$d*t(fiti$v))
   results$A_hat = results$B_hat + X %*% fiti$beta.estim
   results$beta_hat = fiti$beta.estim
   results$rank_A = qr(results$A_hat)$rank
   results$J = best_rank
   
   return(results)
}
