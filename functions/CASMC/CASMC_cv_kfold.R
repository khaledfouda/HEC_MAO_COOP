CASMC_cv_kfold <- function(Y, X_r, W,n_folds=5, lambda.factor=1/4, lambda.init=NA, 
                                n.lambda=30,
                                            trace=FALSE, thresh=1e-6, maxit=100,
                                            rank.init=3, rank.limit=20, rank.step=2,
                                            warm=NULL, tol=1, print.best=TRUE, seed=NULL){

   if(!is.null(seed)) set.seed(seed)
   #----------------------------------------------------
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov_splr(Y, X_r$svdH) * lambda.factor, lambda.init)
   lamseq <- seq(from=lam0, to=1, length=n.lambda)
   #-------------
   # Initialize warm.start for the second model
   svdX = fast.svd(X_r$X)
   Ux = svdX$u
   Vx = svdX$d * t(svdX$v)
   X0 = ginv(t(Vx)%*%Vx) %*% t(Vx)
   warm.start.beta = list()
   warm.start.beta$X1 = X0 %*% t(Ux)
   warm.start.beta$X2 = X0 %*% Vx
   Xinv = ginv(X_r$X)
   #---------------------------------------------------------
   # to update the observed part of xbeta after each fold fit
   observed_indices = which(W==1)
   xbeta.observed = rep(NA, length(observed_indices))
   #-----------------------------------------------------------------------
   # prepare the folds
   folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, W)
   fold_data <- lapply(1:n_folds, function(i) {
      W_fold = folds[[i]]
      valid_ind = W_fold==0 & W==1
      
      obs_ind = which(observed_indices %in% which(W_fold==1))
      
      Y_train = Y * W_fold
      Y_train[W_fold==0] = NA
      Y_train = as(Y_train, "Incomplete")

      Y_valid = Y[valid_ind]

      W_fold[valid_ind] = 1
      W_fold[!valid_ind] = NA
      W_fold <- as(W_fold, "Incomplete")
      virow = W_fold@i
      vpcol = W_fold@p
      W_fold <- NULL

      list(Y_train = Y_train, Y_valid = Y_valid,
           virow = virow, vpcol=vpcol, obs_ind = obs_ind)
   })
   #---------------------------------------------------------------------------
   n <- dim(Y)
   m <- n[2]
   n <- n[1]
   Y[Y == 0] = NA
   Y <- as(Y, "Incomplete")
   xbeta.sparse = Y
   #---------------------------------------------------------------------------
   rank.max <- rank.init
   fiti = NULL
   counter = 0
   best_fit = list(error = Inf)
   niter1 <- niter2 <- 0
   #---------------------------------------------------------------------
   for(i in seq(along=lamseq)) {
      # initial fit to whole data
      if(i != 1) fiti$xbeta.obs <- xbeta.observed 
         # or, initialize it using the secondmodel
         #fiti$xbeta.obs <- suvC(Xv, t(fitx$d * t(fitx$u)), Y@i, Y@p)
      fiti <-  CASMC_fit(y=Y, svdH=X_r$svdH,  trace=F, J=rank.max,
                                    thresh=thresh, lambda=lamseq[i], init = "naive",
                                    final.svd = T,maxit = maxit, warm.start = fiti)
      xbeta.sparse@x <- fiti$xbeta.obs
      #---------
      # prepare warm.start.beta:
      # use naive model to initialize.
      if(i == 1){
         B = t( Xinv %*% naive_MC(as.matrix(xbeta.sparse))) # B = (X^-1 Y)'
         warm.start.beta$Bsvd = fast.svd(B)
      }else warm.start.beta$Bsvd = fitx
      #---------------------------
      # fit second model:
      fitx = SZIRCI(xbeta.sparse, X_r$X, X_r$rank, final.trim = F, thresh=thresh,
                                       warm.start = warm.start.beta, trace.it = F,maxit=maxit)
      Xv = X_r$X %*% fitx$v
      #-------------------------------------------------------------------
      if(trace){
         niter1 <- niter1 + fiti$n_iter
         niter2 <- fitx$n_iter
      }
      #-----------------------------------------------
      err <- rank <- 0
      if(i == 1){
      xbeta.observed <- fiti$xbeta.obs
      }else xbeta.observed <- (fiti$xbeta.obs + xbeta.observed)/2
      #-----------------------------------------------------
      for(fold in 1:n_folds){
         data = fold_data[[fold]]
         #-----
         fiti$xbeta.obs <- suvC(Xv, t(fitx$d * t(fitx$u)), data$Y_train@i, data$Y_train@p)
         fiti <-  CASMC_fit(y=data$Y_train, svdH=X_r$svdH,  trace=F, J=rank.max,
                                       thresh=thresh, lambda=lamseq[i], init = "naive",
                                       final.svd = T,maxit = maxit, warm.start = fiti)
         xbeta.observed[data$obs_ind] <- (xbeta.observed[data$obs_ind] + fiti$xbeta.obs)/2
         #--------------------------------------------------------------
         # predicting validation set and xbetas for next fit:
         Xbvalid = suvC(Xv, t(fitx$d * t(fitx$u)), data$virow, data$vpcol)
         Mvalid = suvC(fiti$u, t(fiti$d * t(fiti$v)), data$virow, data$vpcol)
         #--------------------------------------------
         err = err + test_error(Mvalid+Xbvalid, data$Y_valid)
         rank <- rank + sum(round(fiti$d, 4) > 0) # number of positive sing.values
         if(trace){
            niter1 <- niter1 + fiti$n_iter
         }
         #-----------------------------------------------------------------------------------
      }
      err = err / n_folds
      rank = as.integer(rank / n_folds)
      #------------------------------------------------
      if(trace){
         print(sprintf(paste0("%2d lambda=%9.5g, rank.max = %d  ==>",
                              " rank = %d, error = %.5f, niter = %d + %d"),
                       i, lamseq[i], rank.max, rank, err, niter1,niter2))
         niter1 <- niter2 <- 0
      }
      #-----------------------------------------------------------------------
      if(err < best_fit$error){
         best_fit$error = err
         best_fit$lambda = lamseq[i]
         best_fit$rank.max = rank.max
         best_fit$rank = rank
         best_fit$fit1 = fiti
         best_fit$fit2 = fitx
         best_fit$iter = i
         counter = 0
      }else
         counter = counter + 1
      if(counter >= tol){
         if(trace | print.best)
            print(sprintf("Performance didn't improve for the last %d iterations.", counter))
         break
      }
      #-------------------------------------------------------------------
      rank.max <- min(rank+rank.step, rank.limit)
      #----------------------------------------------------------------

   }
   if(print.best) print(sprintf("lambda=%9.5g, rank.max = %d, error = %.5f",
                                      best_fit$lambda, best_fit$rank.max, best_fit$error))
# 
#    # one last fit!
   best_fit$fit1$xbeta.obs <-xbeta.observed
   # or, initiate with model 2 as:
   # suvC(X_r$X %*% best_fit$fit2$v,t(best_fit$fit2$d * t(best_fit$fit2$u)), Y@i, Y@p)
   best_fit$fit1 <-  CASMC_fit(y=Y, svdH=X_r$svdH,  trace=F, J=best_fit$rank.max,
                                          thresh=thresh, lambda=best_fit$lambda, init = "naive",
                                          final.svd = T,maxit = maxit, warm.start = best_fit$fit1)
   xbeta.sparse@x <- best_fit$fit1$xbeta.obs
   B = t( Xinv %*% naive_MC(as.matrix(xbeta.sparse))) # B = (X^-1 Y)'
   warm.start.beta$Bsvd = fast.svd(B)
   best_fit$fit2 = SZIRCI(xbeta.sparse, X_r$X, X_r$rank, final.trim = F, thresh=thresh,
                                             warm.start = warm.start.beta, trace.it = F,maxit=maxit)


   return(best_fit)
}

