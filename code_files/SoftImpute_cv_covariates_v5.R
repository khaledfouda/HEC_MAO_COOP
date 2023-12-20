
simpute.cov.cv.lambda1 <- function(Y, X, W, Y.valid, lambda2=NA, J=NA,
                              trace=FALSE, print.best=TRUE, thresh=1e-5,
                           type="als", lambda1.grid=NA, n1n2=1, warm=NULL){
   
   stopifnot(type %in% c("svd", "als"))
   if(type == "svd"){
      fit.function <- simpute.svd.cov
   }else
      fit.function <- simpute.als.cov
   
   n <- ncol(X)
   stopifnot(n1n2 %in% 1:3)
   # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
   
   ymiss <- W == 0
   #rank.max <- rank.init
   counter <- 1
   X.X = t(X) %*% X
   if(n1n2 == 2){
      n1n2 = svd(X)$d[1]
   }else if(n1n2 == 3){
      n1n2 = nrow(Y) * ncol(Y)
   }
   
   
   best_lambda1 = NA
   best_error = Inf
   
   fiti = warm
   
   for(i in 1:length(lambda1.grid)){
      #print('hi')
      beta_partial = solve(X.X +  diag(n1n2*lambda1.grid[i], n)) %*% t(X)
      fiti <- fit.function(Y, X, beta_partial, thresh=thresh, lambda = lambda2, J=J, warm.start = fiti)
      test_estim = (fiti$u %*%(fiti$d*t(fiti$v)))[ymiss] + (X %*% fiti$beta.estim)[ymiss]
      err = test_error(test_estim, Y.valid)
         
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
