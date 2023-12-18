
lambda0.cov <- function(Y, X){
   
   n1 <- dim(Y)[1]
   #n2 <- dim(Y)[2]
   #m1 <- dim(X)[2]
   ynas <- Y==0 #is.na(Y)
   
   # The following two lines are as shown in (c) and (d)
   X.X = t(X) %*% X
   P_X = X %*% solve(X.X) %*% t(X)
   P_bar_X = diag(1,n1) - P_X
   #theta_hat = theta_estimator(W=W, X=X)
   beta_partial = solve(X.X) %*% t(X)
   
   yfill <- Y
   yfill[ynas] <- 0
   beta.estim <- beta_partial %*% yfill
   Xbeta <- X %*% beta.estim
   yplus <- yfill - Xbeta 
   return(svd(yplus)$d[1])
}

simpute.cov.cv.L2 <- function(Y, X, W, Y.valid, lambda.factor=1/4, lambda.init=NA, n.lambda=20,
                              trace=FALSE, print.best=TRUE, tol=5, thresh=1e-5,
                           rank.init=10, rank.limit=50, rank.step=2,
                           type="als", lambda1=0, lambda1.grid=NA, n1n2=1, warm=NULL){
   
   stopifnot(type %in% c("svd", "als"))
   if(type == "svd"){
      fit.function <- simpute.svd.cov
   }else
      fit.function <- simpute.als.cov
   
   n <- ncol(X)
   stopifnot(n1n2 %in% 1:3)
   # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
   #Y[Y==0] = NA
   #xs <- as(Y, "Incomplete")
   
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov(Y, X) * lambda.factor, lambda.init) 
   #lam0 <- 40 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   
   fits <- as.list(lamseq)
   ranks <- as.integer(lamseq)
   
   
   rank.max <- rank.init
   #warm <- NULL
   best_estimates = NA
   best_fit <- list(error=Inf, rank_A=NA, rank_B=NA, lambda=NA, rank.max=NA)
   counter <- 1
   X.X = t(X) %*% X
   if(n1n2 == 2){
      n1n2 = svd(X)$d[1]
   }else if(n1n2 == 3){
      n1n2 = nrow(Y) * ncol(Y)
   }
   
   beta_partials = list()
   beta_partial_init = solve(X.X) %*% t(X)
   for(i in 1:length(lambda1.grid))
      beta_partials[[i]] = solve(X.X +  diag(n1n2*lambda1.grid[i], n)) %*% t(X)
   
   best_overall_lambda1 = NA
   best_overall_lambda2 = NA
   best_overall_error = Inf
   best_overall_J = NA
   best_overall_lambda1.index = NA
   
   
   for(i in seq(along=lamseq)) {
      #print('hi')
      fiti <- fit.function(Y, X, beta_partial_init, thresh=thresh, lambda = lamseq[i], J=rank.max, warm.start = warm)
      #print("hi2")
      # compute rank.max for next iteration
      rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
      rank.max <- min(rank+rank.step, rank.limit)
      soft_estim = fiti$u %*%(fiti$d*t(fiti$v)) + X %*% fiti$beta.estim
      err = test_error(soft_estim[W==0], Y.valid)
      warm <- fiti # warm start for next 
      
      best_lambda1 = 0
      best_lambda1_error = err
      
      for(j in seq(along=lambda1.grid)){
         fiti <- fit.function(Y, X, beta_partials[[j]], thresh=thresh, lambda = lamseq[i], J=rank.max, warm.start = warm)
         soft_estim = fiti$u %*%(fiti$d*t(fiti$v)) + X %*% fiti$beta.estim
         err = test_error(soft_estim[W==0], Y.valid)
         
         if(err < best_lambda1_error){
            best_lambda1 = lambda1.grid[j]
            best_lambda1.index = j
            best_lambda1_error = err
            rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
            rank.max <- min(rank+rank.step, rank.limit)
            warm <- fiti
         }
      }
      if(best_lambda1_error < best_overall_error){
         best_overall_lambda1 = best_lambda1
         best_overall_lambda2 = lamseq[i]
         best_overall_error = best_lambda1_error
         best_overall_J = rank.max
         best_overall_lambda1.index = best_lambda1.index
         counter = 1
      }else counter = counter + 1

      
      
      # get test estimates and test error
      #v=as.matrix(fiti$v)
      #vd=v*outer(rep(1,nrow(v)),fiti$d)
      #soft_estim = fiti$u %*% t(vd)  + X %*% fiti$beta.estim
      #----------------------------
      if(trace==TRUE)
         cat(sprintf("%2d lambda1=%9.5g, lambda2=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
                     i, best_lambda1, lamseq[i], rank.max, rank, best_lambda1_error))
      #-------------------------
      # register best fir
      # if(err < best_fit$error){
      #    best_fit$error = err
      #    best_fit$rank_A = qr(soft_estim)$rank
      #    best_fit$rank_B = rank
      #    best_fit$lambda = lamseq[i]
      #    best_fit$rank.max = rank.max
      #    best_estimates = soft_estim
      #    best_beta = fiti$beta.estim
      #    best_B = fiti$u %*% t(vd)
      #    counter=1
      # }else counter = counter + 1
      if(counter >= tol){
         cat(sprintf("Performance didn't improve for the last %d iterations.", counter))
         break
      }
   }
   if(print.best==TRUE) 
      cat(sprintf("BEST lambda1=%9.5g, lambda2=%9.5g, rank.max = %d, error = %.5f\n",
                  best_overall_lambda1, best_overall_lambda2, best_overall_J, best_overall_error))
   # run final fit with the given parameter:
   if(best_overall_lambda1 == 0){
      best_beta_partial = beta_partial_init
   }else best_beta_partial = beta_partials[[best_overall_lambda1.index]]
   
   fiti <- fit.function(Y, X, best_beta_partial, thresh=thresh, lambda = best_overall_lambda2, J=best_overall_J,
                        warm.start = warm)
   
   
   results = list()
   results$B_hat = fiti$u %*%(fiti$d*t(fiti$v))
   results$A_hat = results$B_hat + X %*% fiti$beta.estim
   results$beta_hat = fiti$beta.estim
   #best_fit$last.fit = fiti
   return(results)
}
#--------------------------------------------------------------------
# softimpute with validation and no covariates using original method.

simpute.orig <- function(Y, W, A, n.lambda=20,
                               trace=FALSE, print.best=TRUE, tol=5, thresh=1e-5, rank.init=10, rank.limit=50, rank.step=2){
   # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
   Y[Y==0] = NA
   #xs <- as(Y, "Incomplete")
   lam0 <- lambda0(Y)
   #lam0 <- 40 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   
   fits <- as.list(lamseq)
   ranks <- as.integer(lamseq)
   
   
   rank.max <- rank.init
   warm <- NULL
   best_estimates = NA
   best_fit <- list(error=Inf, rank_A=NA, rank_B=NA, lambda=NA, rank.max=NA)
   counter <- 1
   
   for(i in seq(along=lamseq)) {
      
      fiti <- softImpute(Y, lambda=lamseq[i], rank.max=rank.max, warm=warm, thresh=thresh)
      
      # compute rank.max for next iteration
      rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
      rank.max <- min(rank+rank.step, rank.limit)
      
      # get test estimates and test error
      soft_estim = complete(Y, fiti)
      err = test_error(soft_estim[W==0], A[W==0])
      #----------------------------
      warm <- fiti # warm start for next 
      if(trace==TRUE)
         cat(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
                     i, lamseq[i], rank.max, rank, err))
      #-------------------------
      # register best fir
      if(err < best_fit$error){
         best_fit$error = err
         best_fit$rank_A = qr(soft_estim)$rank
         best_fit$rank_B = rank
         best_fit$lambda = lamseq[i]
         best_fit$rank.max = rank.max
         best_estimates = soft_estim
         counter=1
      }else counter = counter + 1
      if(counter >= tol){
         cat(sprintf("Performance didn't improve for the last %d iterations.", counter))
         break
      }
   }
   if(print.best==TRUE) print(best_fit)
   best_fit$A_hat = best_estimates
   best_fit$B_hat = NA
   best_fit$beta_hat = NA
   return(best_fit)
}
#-------------------------------------------------------------------
softImputeALS_L2 <- function(Y.train, Y.valid, W.valid, X, lambda1.grid=seq(0,20,length.out=10),
                             n1n2=1, no_cores=NA, max_cores=20, rank.limit=30){
   
   if(is.na(no_cores)) no_cores = length(lambda1.grid) + 1
   no_cores = min(max_cores, no_cores)
   
   if(no_cores > 1){
      model_output <- mclapply(lambda1.grid, function(lambda) {
         simpute.cov.cv(Y.train, X, W.valid, Y.valid,
                        trace=FALSE, rank.limit = rank.limit, lambda1=lambda,n1n2 = n1n2,print.best = FALSE)
      }, mc.cores = no_cores)
   }else{
      model_output = list()
      for(i in 1:length(lambda1.grid)){
         model_output[[i]] <- simpute.cov.cv(Y.train, X, W.valid, Y.valid,
                                             trace=FALSE, rank.limit = rank.limit,
                                             lambda1=lambda1.grid[i],n1n2 = n1n2, print.best = FALSE)
      }
   }
   
   valid_errors <- unlist(lapply(model_output, function(d) d$error))
   best_index = which.min(valid_errors)
   best_fit <- model_output[[best_index]]
   list(best_score=valid_errors[best_index], best_fit=best_fit, lambda1 = lambda1.grid[best_index])
}





