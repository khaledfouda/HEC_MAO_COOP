CAMC_cv_holdout <- function(Y_train,
                            X,
                            W_valid,
                            Y_valid,
                            trace = TRUE,
                            rank.limit = 30,
                            print.best = TRUE,
                            lambda.1_grid = seq(0, 3, length = 20),
                            rank.step = 2,
                            n.lambda = 30,
                            type = "als",
                            quiet = FALSE,
                            tol = 2) {
   start_time = Sys.time()
   fit1 <-
      CAMC_cv_holdout_lambda2(
         Y_train,
         X,
         W_valid,
         Y_valid,
         trace = trace,
         rank.limit = rank.limit,
         print.best = print.best,
         rank.step = rank.step,
         type = type,
         lambda1 = 0,
         tol = tol,
         n.lambda = n.lambda,
         quiet = quiet
      )
   fit2 <-
      CAMC_cv_holdout_lambda1(
         Y_train,
         X,
         W_valid,
         Y_valid,
         fit1$lambda,
         fit1$rank.max,
         print.best = print.best,
         trace = trace,
         lambda1.grid = lambda.1_grid ,
         n1n2 = 1,
         warm = NULL
      )
   time_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
   list(fit1 = fit1,
        fit2 = fit2,
        time_seconds = time_seconds)
}



#-=--------------------------------------------------------------------
CAMC_cv_holdout_lambda2 <-
   function(Y,
            X,
            W,
            Y.valid,
            lambda.factor = 1 / 4,
            lambda.init = NA,
            n.lambda = 20,
            trace = FALSE,
            print.best = TRUE,
            tol = 5,
            thresh = 1e-6,
            rank.init = 10,
            rank.limit = 50,
            rank.step = 2,
            type = "als",
            lambda1 = 0,
            n1n2 = 1,
            warm = NULL,
            quiet = FALSE,
            maxit = 300) {
      stopifnot(type %in% c("svd", "als"))
      if (type == "svd") {
         fit.function <- CAMC_fit_svd
      } else
         fit.function <- CAMC_fit_als
      
      stopifnot(n1n2 %in% 1:3)
      # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
      #Y[Y==0] = NA
      #xs <- as(Y, "Incomplete")
      
      lam0 <-
         ifelse(is.na(lambda.init),
                lambda0.cov(Y, X) * lambda.factor,
                lambda.init)
      #lam0 <- 40
      lamseq <- seq(from = lam0,
                    to = 0,
                    length = n.lambda)
      
      
      rank.max <- rank.init
      #warm <- NULL
      best_estimates = NA
      best_fit <-
         list(
            error = Inf,
            rank_O = NA,
            rank_M = NA,
            lambda = NA,
            rank.max = NA
         )
      counter <- 1
      X.X = t(X) %*% X
      # I dropped this part so n1n2=1
      if (n1n2 == 2) {
         n1n2 = svd(X)$d[1]
      } else if (n1n2 == 3) {
         n1n2 = nrow(Y) * ncol(Y)
      }
      #-----
      beta_partial = ginv(X.X +  diag(n1n2 * lambda1, ncol(X))) %*% t(X)
      for (i in seq(along = lamseq)) {
         fiti <-
            fit.function(
      Y,
      X,
      beta_partial,
      thresh = thresh,
      lambda = lamseq[i],
      J = rank.max,
      warm.start = warm,
      maxit = maxit
            )
         
         
         # get test estimates and test error
         v = as.matrix(fiti$v)
         vd = v * outer(rep(1, nrow(v)), fiti$d)
         soft_estim = fiti$u %*% t(vd)  + X %*% fiti$beta.estim
         err = test_error(soft_estim[W == 0], Y.valid)
         rank <-
            sum(round(fiti$d, 4) > 0) # number of positive sing.values
         #----------------------------
         warm <- fiti # warm start for next
         if (trace == TRUE)
            print(
               sprintf(
                  "%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
                  i,
                  lamseq[i],
                  rank.max,
                  rank,
                  err
               )
            )
         #-------------------------
         # register best fir
         if (err < best_fit$error) {
            best_fit$error = err
            best_fit$rank_O = qr(soft_estim)$rank
            best_fit$rank_M = rank
            best_fit$lambda = lamseq[i]
            best_fit$rank.max = rank.max
            best_estimates = soft_estim
            best_beta = fiti$beta.estim
            best_M = fiti$u %*% t(vd)
            counter = 1
         } else
            counter = counter + 1
         if (counter >= tol) {
            if (quiet == FALSE)
               print(sprintf(
                  "Performance didn't improve for the last %d iterations.",
                  counter
               ))
            break
         }
         # compute rank.max for next iteration
         rank.max <- min(rank + rank.step, rank.limit)
      }
      if (print.best == TRUE)
         print(best_fit)
      best_fit$estimates = best_estimates
      best_fit$M = best_M
      best_fit$beta = best_beta
      best_fit$last.fit = fiti
      return(best_fit)
   }
#--------------------------------------------------------------------
CAMC_cv_holdout_lambda1 <-
   function(Y,
            X,
            W,
            Y.valid,
            lambda2 = NA,
            J = NA,
            trace = FALSE,
            print.best = TRUE,
            thresh = 1e-5,
            type = "als",
            lambda1.grid = NA,
            n1n2 = 1,
            warm = NULL) {
      stopifnot(type %in% c("svd", "als"))
      if (type == "svd") {
         fit.function <- CAMC_fit_svd
      } else
         fit.function <- CAMC_fit_als
      
      n <- ncol(X)
      stopifnot(n1n2 %in% 1:3)
      # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
      
      ymiss <- W == 0
      #rank.max <- rank.init
      counter <- 1
      X.X = t(X) %*% X
      if (n1n2 == 2) {
         n1n2 = svd(X)$d[1]
      } else if (n1n2 == 3) {
         n1n2 = nrow(Y) * ncol(Y)
      }
      
      
      best_lambda1 = NA
      best_error = Inf
      
      fiti = warm
      
      for (i in 1:length(lambda1.grid)) {
         #print('hi')
         beta_partial = ginv(X.X +  diag(n1n2 * lambda1.grid[i], n)) %*% t(X)
         fiti <-
            fit.function(
      Y,
      X,
      beta_partial,
      thresh = thresh,
      lambda = lambda2,
      J = J,
      warm.start = fiti
            )
         test_estim = (fiti$u %*% (fiti$d * t(fiti$v)))[ymiss] + (X %*% fiti$beta.estim)[ymiss]
         err = test_error(test_estim, Y.valid)
         
         if (err < best_error) {
            best_lambda1 = lambda1.grid[i]
            best_index = i
            best_error = err
         }
         if (trace == TRUE)
            cat(
               sprintf(
                  "%2d lambda1=%9.5g, lambda2=%9.5g, rank.max = %d, error = %.5f\n",
                  i,
                  best_lambda1,
                  lambda2,
                  J,
                  best_error
               )
            )
      }
      
      beta_partial = ginv(X.X +  diag(n1n2 * lambda1.grid[best_index], n)) %*% t(X)
      fiti <-
         fit.function(
      Y,
      X,
      beta_partial,
      thresh = thresh,
      lambda = lambda2,
      J = J,
      warm.start = fiti
         )
      
      
      results = list()
      #results$fit = fiti
      results$lambda1 = best_lambda1
      results$lambda2 = lambda2
      results$M = fiti$u %*% (fiti$d * t(fiti$v))
      results$estimates = results$M + X %*% fiti$beta.estim
      results$beta = fiti$beta.estim
      results$rank_O = qr(results$estimates)$rank
      results$J = J
      return(results)
   }
