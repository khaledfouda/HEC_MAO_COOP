
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
   
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov(Y, X) * lambda.factor, lambda.init) 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   
   #fits <- as.list(lamseq)
   #ranks <- as.integer(lamseq)
   #warm <- NULL
   #best_estimates = NA
   #best_fit <- list(error=Inf, rank_A=NA, rank_B=NA, lambda=NA, rank.max=NA)
   
   ymiss <- W == 0
   rank.max <- rank.init
   counter <- 1
   X.X = t(X) %*% X
   if(n1n2 == 2){
      n1n2 = svd(X)$d[1]
   }else if(n1n2 == 3){
      n1n2 = nrow(Y) * ncol(Y)
   }
   
   beta_partials = list()
   #beta_partial_init = solve(X.X) %*% t(X)
   for(i in 1:length(lambda1.grid))
      beta_partials[[i]] = solve(X.X +  diag(n1n2*lambda1.grid[i], n)) %*% t(X)
   
   best_overall_lambda1 = NA
   best_overall_lambda2 = NA
   best_overall_error = Inf
   best_overall_J = NA
   best_overall_lambda1.index = NA
   best_overall_fit = NULL
   rank_maxes = rep(rank.init, length(lambda1.grid))
   
   
   for(i in seq(along=lamseq)) {
      #print('hi')
      fiti <- fit.function(Y, X, beta_partials[[1]], thresh=thresh, lambda = lamseq[i], J=rank_maxes[1], warm.start = warm)
      
      # compute rank.max for next iteration
      rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
      #rank.max <- min(rank+rank.step, rank.limit)
      rank_maxes[1] <- min(rank+rank.step, rank.limit)
      test_estim = (fiti$u %*%(fiti$d*t(fiti$v)))[ymiss] + (X %*% fiti$beta.estim)[ymiss]
      err = test_error(test_estim, Y.valid)
      warm <- fiti # warm start for next 
      
      best_lambda1 = 0
      best_lambda1_error = err
      best_lambda1.index = -1
      #best_lambda1_fit = fiti
      
      for(j in 2:length(lambda1.grid)){
         fiti <- fit.function(Y, X, beta_partials[[j]], thresh=thresh, lambda = lamseq[i], J=rank_maxes[j], warm.start = warm)
         test_estim = (fiti$u %*%(fiti$d*t(fiti$v)))[ymiss] + (X %*% fiti$beta.estim)[ymiss]
         err = test_error(test_estim, Y.valid)
         
         if(err < best_lambda1_error){
            best_lambda1 = lambda1.grid[j]
            best_lambda1.index = j
            best_lambda1_error = err
            rank.max <- rank_maxes[j] #min(rank+rank.step, rank.limit)
            #warm <- fiti
            #best_lambda1_fit <- fiti
         }
         rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
         rank_maxes[j] <- min(rank+rank.step, rank.limit)
      }
      if(best_lambda1_error < best_overall_error){
         best_overall_lambda1 = best_lambda1
         best_overall_lambda2 = lamseq[i]
         best_overall_error = best_lambda1_error
         best_overall_J = rank.max
         best_overall_lambda1.index = best_lambda1.index
         #best_overall_fit = best_lambda1_fit
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
   best_overall_fit <- fiti
   results$fit = fiti#best_overall_fit
   results$lambda1 = best_overall_lambda1
   results$lambda2 = best_overall_lambda2
   results$B_hat = best_overall_fit$u %*%(best_overall_fit$d*t(best_overall_fit$v))
   results$A_hat = results$B_hat + X %*% best_overall_fit$beta.estim
   results$beta_hat = best_overall_fit$beta.estim
   results$rank_A = qr(results$A_hat)$rank
   return(results)
}
#--------------------------------------------------------------------
