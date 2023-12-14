#' Hyperparameter optimization for $\lambda_1,\lambda_2,\alpha$ is done using k-fold (k=5 by default) and grid search.
#'  Mao's paper optimizes for each parameter separately while fixing the other two.
#'   The explored/recommended the range of (0-2) for $\lambda_1$, (0.1,0.9) for $\lambda_2$, and (0.992,1) for $\alpha$.

prepare_fold_data <- function(Y, X, W, A, n1n2_optimized, theta_estimator) {
  n1 = dim(Y)[1]
  n2 = dim(Y)[2]
  m  = dim(X)[2]
  
  # The following two lines are as shown in (c) and (d)
  X.X = t(X) %*% X
  P_X = X %*% solve(X.X) %*% t(X)
  P_bar_X = diag(1,n1) - P_X 

  theta_hat = theta_estimator(W=W, X=X)
  
  #---------
  # The following are partial parts of equations 8 and 11 that don't involve the hyperparameters.
   # this is useful to avoid unneccessary matrix multiplications.
   #----------
   if(n1n2_optimized == TRUE){
      # this one is for equation 8, the product n1n2 is replace with the Eigen value
      n1n2Im = svd(X.X)$d[1]  * diag(1, m)  
   }else{
      n1n2Im = n1 * n2  * diag(1, m)
   }
   # the following is the product of W * theta_hat * Y
   W_theta_Y = Y * theta_hat
   X.W.theta.Y = t(X) %*% W_theta_Y
   svdd = svd(P_bar_X %*% W_theta_Y)
   if(n1n2_optimized == TRUE){
      # this one is for equation 11, the product is also replace with the Eigen value of the SVD
      n1n2 = svdd$d[1]
   }else{
      n1n2 = n1 * n2
   }
   
   return(list(A.test=A[W==0], W=W, X=X, X.X=X.X, n1n2Im=n1n2Im, n1n2=n1n2,
               X.W.theta.Y = X.W.theta.Y, svdd=svdd))
}

Mao.fit_optimized <- function(data, lambda.1, lambda.2, alpha){
   
   beta_hat = solve( data$X.X + data$n1n2Im * lambda.1) %*% data$X.W.theta.Y
   T_c_D = data$svdd$u %*% (pmax(data$svdd$d - alpha*data$n1n2*lambda.2, 0) * t(data$svdd$v))
   # B hat as in (11)
   B_hat = (1 + 2 * (1-alpha) * data$n1n2 * lambda.2 )^(-1) * T_c_D
   # Estimate the matrix as given in the model at the top
   A_hat = data$X %*% beta_hat + B_hat
   
   return(A_hat[data$W == 0])
}

# Computing the test error as given by Mao in page 205
test_error <- function(A.hat.test, A.test){
   return(sum( (A.hat.test-A.test)^2 )/ sum((A.test-mean(A.test))^2))
}

Mao.cv <- function(A, X, Y, W, n_folds=5, lambda.1_grid = seq(0,1,length=30),
                   lambda.2_grid = seq(0.9, 0.1, length=30),
                   alpha_grid = seq(0.992, 1, length=20), seed=2023, numCores=16, n1n2_optimized=TRUE,
                   theta_estimator=theta_default){
   
   #' -------------------------------------------------------------------
   #' Input :
   #' A : Complete (True) A matrix as in the model above of size n1 by n2
   #' X :  Covariate matrix of size  n1 by m
   #' W : Binary matrix representing the mask. wij=1 if yij is observed. size similar to A
   #' The rest are cross validation parameters
   #' --------------------------------------------------------------------
   #' Output:
   #' list of best parameters and best score (minimum average MSE across folds)
   #' --------------------------------------------------------------------
   
   set.seed(seed = seed)
   indices = sample(cut(seq(1, nrow(A)), breaks=n_folds, labels=FALSE))
   best_score = Inf
   best_params = list(alpha = NA, lambda.1 = NA, lambda.2 = NA)
   
   fold_data = lapply(1:n_folds, function(i) {
      train_indices = which(indices != i, arr.ind = TRUE)
      Y_train = Y[train_indices,]
      X_train = X[train_indices,]
      W_train = W[train_indices,]
      A_train = A[train_indices,]
      prepare_fold_data(Y_train, X_train, W_train, A_train, n1n2_optimized, theta_estimator)
   })
   
   # Original for loop to be executed on one node 
   # *******************************************
   # for(alpha in alpha_grid){
   #   for(lambda.1 in lambda.1_grid){
   #     for(lambda.2 in lambda.2_grid){
   #       
   #       scores = numeric(n_folds)
   #       for(i in 1:n_folds){
   #         data = fold_data[[i]]
   #         #train_indices = which(indices != i, arr.ind=TRUE)
   #         
   #         # compute the estimates with a modified fit function
   #         A_hat = Mao.fit_optimized(data, lambda.1, lambda.2, alpha)
   #         #print(dim(A[test_indices,]))
   #         #print(dim(A_hat))
   #         #print(dim(data$X))
   #         # Evaluate model performance using MSE
   #         # IMPORTANT: As this is not a predictive model, the MSE is calculated on the training data
   #         # that is, the test fold isn't used at all.
   #         scores[i] = mean((data$A - A_hat)^2)
   #       }
   #       avg_score = mean(scores)
   #       
   #       if(avg_score < best_score){
   #         best_score = avg_score
   #         best_params = list(alpha=alpha, lambda.1=lambda.1, lambda.2=lambda.2)
   #       }
   #     }
   #   }
   # }
   # ************************************************************
   if(numCores == 1){
      
      # fixing optimal values of lambda 1 and alpha and optimizing for alpha separately
      lambda.1 = 0
      alpha = 1
      best_score = Inf
      for(lambda.2 in lambda.2_grid){
         score = 0
         for(i in 1:n_folds){
            data = fold_data[[i]]
            # compute the estimates with a modified fit function
            A_hat_test = Mao.fit_optimized(data, lambda.1, lambda.2, alpha)
            # Evaluate model performance using MSE
            # IMPORTANT: As this is not a predictive model, the MSE is calculated on the training data
            # that is, the test fold isn't used at all.
            #scores[i] = mean((data$A.test - A_hat_test)^2)
            # -- EDIT: Using Mao's formula in page 205 to compute the test error
            score = score + test_error(A_hat_test, data$A.test)
         }
         score = score / n_folds
         
         if(score < best_score){
            best_score = score
            best_params$lambda.2 = lambda.2
         }
      }
      # fixing optimal values of lambda 2 and lambda 1 and optimizing for alpha separately
      lambda.2 = best_params$lambda.2
      lambda.1 = 0
      best_score = Inf
      for(alpha in alpha_grid){
         score = 0
         for(i in 1:n_folds){
            data = fold_data[[i]]
            # compute the estimates with a modified fit function
            A_hat_test = Mao.fit_optimized(data, lambda.1, lambda.2, alpha)
            # Evaluate model performance using MSE
            # IMPORTANT: As this is not a predictive model, the MSE is calculated on the training data
            # that is, the test fold isn't used at all.
            #scores[i] = mean((data$A.test - A_hat_test)^2)
            # -- EDIT: Using Mao's formula in page 205 to compute the test error
            score = score + test_error(A_hat_test, data$A.test)
         }
         score = score / n_folds
         
         if(score < best_score){
            best_score = score
            best_params$alpha = alpha
         }
      }
   
   }else{
   # prepare the cluster
   cl <- makeCluster(numCores) 
   registerDoParallel(cl)
   # fixing lambda 1 at 0 and optimizing for lambda 2 and alpha using a grid
   # Export the Mao.fit_optimized function and any other necessary objects to each worker
   clusterExport(cl, varlist = c("Mao.fit_optimized","test_error"))
   results <- foreach(alpha = alpha_grid, .combine = rbind) %:%
      foreach(lambda.2 = lambda.2_grid, .combine = rbind) %dopar% {
         lambda.1 = 0
         score = 0
         for (i in 1:n_folds) {
            data = fold_data[[i]]
            A_hat_test = Mao.fit_optimized(data, lambda.1, lambda.2, alpha)
            # scores[i] = mean((data$A.test - A_hat_test)^2)
            # -- EDIT: Using Mao's formula in page 205 to compute the test error
            score = score + test_error(A_hat_test, data$A.test)
         }
         #score = score / n_folds
         c(alpha, lambda.2, score)
      }
   # Process results to find the best parameters
   # Edited on Dec 1st to pick the minimum score with highest lambda.2 value.
   min_score <- min(results[, 3])
   # Subset to only include results with the minimum score
   min_results <- results[results[, 3] == min_score, , drop = FALSE] # drop to keep it as df
   # In case of multiple results with the same score, find the one with the highest lambda.2
   if (nrow(min_results) > 1) {
      best_result <- min_results[which.max(min_results[, 2]), ]
   } else {
      best_result <- min_results  # If only one row, it's already the best result
   }
   #best_result <- results[which.min(results[, 3]), ] # old line
   # Extract the best parameters
   best_params <- list(alpha = best_result[1], lambda.1 = 0, lambda.2 = best_result[2])
   best_score <- best_result[3] / n_folds
   # close the cluster
   stopCluster(cl)
   }
   #--------------------------------------------
   # fixing optimal values of lambda 2 and alpha and optimizing for lambda 1 separately
   lambda.2 = best_params$lambda.2
   alpha = best_params$alpha
   best_score = Inf
   for(lambda.1 in lambda.1_grid){
      score = 0
      for(i in 1:n_folds){
         data = fold_data[[i]]
         # compute the estimates with a modified fit function
         A_hat_test = Mao.fit_optimized(data, lambda.1, lambda.2, alpha)
         # Evaluate model performance using MSE
         # IMPORTANT: As this is not a predictive model, the MSE is calculated on the training data
         # that is, the test fold isn't used at all.
         #scores[i] = mean((data$A.test - A_hat_test)^2)
         # -- EDIT: Using Mao's formula in page 205 to compute the test error
         score = score + test_error(A_hat_test, data$A.test)
      }
      score = score / n_folds
      
      if(score < best_score){
         best_score = score
         best_params$lambda.1 = lambda.1
      }
   }
   
   #---------------------------------------------------
   return(list(best_parameters = best_params, best_score = best_score))
   
}
