#' Hyperparameter optimization for $\lambda_1,\lambda_2,\alpha$ is done using k-fold (k=5 by default) and grid search.
#'  Mao's paper optimizes for each parameter separately while fixing the other two.
#'   The explored/recommended the range of (0-2) for $\lambda_1$, (0.1,0.9) for $\lambda_2$, and (0.992,1) for $\alpha$.

k_fold_cells <- function(n_rows, n_cols, n_folds, W, min_per_row=15, max_iter=100) {
   # Create a data frame of all matrix indices
   indices <- expand.grid(row = 1:n_rows, col = 1:n_cols)
   # we only consider non-missing data (ie, with Wij=1)
   indices <- indices[W==1,]
   # Shuffle indices (both rows and columns are shuffled. later, we will reshuffle the columns)
   indices <- indices[sample(1:nrow(indices)), ]
   
   # Assign each observed index to one of k groups, ensuring each row has at least 5 samples
   indices <- indices %>%
      group_by(row) %>%
      # the following is to shuffle within each row
      do(sample_n(., size = nrow(.))) %>% 
      mutate(fold = rep(1:n_folds, length.out = n())) %>%
      ungroup()
   
   # Ensure that each row has at least min_per_row samples
   # counter = 0
   # best_occ = -1 
   # while( best_occ < min_per_row) {
   #    temp = sample(indices$fold)
   #    min_occ = min(table(indices$row, temp))
   #    if(min_occ > best_occ) { 
   #       best_occ = min_occ
   #       indices$fold <- temp
   #    }
   #    #indices$fold <- sample(indices$fold)
   #    counter = counter + 1
   #    if(counter>=max_iter)
   #    {
   #       print(paste0("failed to keep ",min_per_row, " obs per row after ", counter, " iterations. ",
   #                    "Minimum number of occurences is: ",min(table(indices$row, indices$fold)),
   #                    ". Attempting ", min_per_row-5, " obs per row..."))
   #       min_per_row = min_per_row - 5
   #       counter = 0
   #    } 
   # }
   
   # Assign each index to one of k groups
   #indices$fold <- rep(1:n_folds, length.out = nrow(indices))
   # Create a list to hold each fold
   folds <- vector("list", n_folds)
   for (i in 1:n_folds) {
      # Create a mask for the test cells in this fold
      # 1 -> train | 0 -> valid | 0 -> test (missing)
      valid_mask <- matrix(1, nrow = n_rows, ncol = n_cols)
      test_indices <- indices[indices$fold == i, ]
      valid_mask[W==0] <- 0
      valid_mask[as.matrix(test_indices[, c("row", "col")])] <- 0
      # Store the mask
      folds[[i]] <- valid_mask
   }
   return(folds)
}


prepare_fold_data <- function(Y_train, Y_valid, W_fold, W, X, n1n2_optimized, theta_estimator) {
  n1 = dim(Y_train)[1]
  n2 = dim(Y_train)[2]
  m  = dim(X)[2]
  
  # The following two lines are as shown in (c) and (d)
  X.X = t(X) %*% X
  P_X = X %*% solve(X.X) %*% t(X)
  P_bar_X = diag(1,n1) - P_X 

  theta_hat = theta_estimator(W=W_fold, X=X)
  
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
   W_theta_Y = Y_train * theta_hat
   X.W.theta.Y = t(X) %*% W_theta_Y
   svdd = svd(P_bar_X %*% W_theta_Y)
   if(n1n2_optimized == TRUE){
      # this one is for equation 11, the product is also replace with the Eigen value of the SVD
      n1n2 = svdd$d[1]
   }else{
      n1n2 = n1 * n2
   }
   
   return(list(Y_valid=Y_valid, W_fold=W_fold, W=W, X=X, X.X=X.X, n1n2Im=n1n2Im, n1n2=n1n2,
               X.W.theta.Y = X.W.theta.Y, svdd=svdd))
}

Mao.fit_optimized <- function(data, lambda.1, lambda.2, alpha){
   
   beta_hat = solve( data$X.X + data$n1n2Im * lambda.1) %*% data$X.W.theta.Y
   T_c_D = data$svdd$u %*% (pmax(data$svdd$d - alpha*data$n1n2*lambda.2, 0) * t(data$svdd$v))
   # B hat as in (11)
   B_hat = (1 + 2 * (1-alpha) * data$n1n2 * lambda.2 )^(-1) * T_c_D
   # Estimate the matrix as given in the model at the top
   A_hat = data$X %*% beta_hat + B_hat
   
   return(A_hat[data$W_fold == 0 & data$W == 1])
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
   #indices = sample(cut(seq(1, nrow(A)), breaks=n_folds, labels=FALSE))
   best_score = Inf
   best_params = list(alpha = NA, lambda.1 = NA, lambda.2 = NA)
   
   folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, W)
   
   fold_data = lapply(1:n_folds, function(i) {
      #train_indices = which(indices != i, arr.ind = TRUE)
      W_fold = folds[[i]] #W[train_indices,]
      #---------------------------------------------------------------
      # EDIT: I implemented this above in k_fold_cells() no longer needed
      #W_fold[W==0] = 1 # This to avoid having the missing data as test set. 
      # Note that we don't have their original values so if they're passed to the validation step, 
      # their original will be equal to 0. We hope we have enough W_fold = 0 while W = 1.
      #---------------------------------------------------
      Y_train = Y * W_fold
      Y_valid = Y[W_fold==0 & W==1]
      prepare_fold_data(Y_train, Y_valid, W_fold, W, X, n1n2_optimized, theta_estimator)
   })
   
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
            # -- EDIT: Using Mao's formula in page 205 to compute the test error
            score = score + test_error(A_hat_test, data$Y_valid)
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
            # -- EDIT: Using Mao's formula in page 205 to compute the test error
            score = score + test_error(A_hat_test, data$Y_valid)
         }
         score = score / n_folds
         
         if(score < best_score){
            best_score = score
            best_params$alpha = alpha
         }
      }
   
   }else{ # Run on multiple cores
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
            score = score + test_error(A_hat_test, data$Y_valid)
         }
         score = score / n_folds
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
      best_score <- best_result[3]
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
         # -- EDIT: Using Mao's formula in page 205 to compute the test error
         score = score + test_error(A_hat_test, data$Y_valid)
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
