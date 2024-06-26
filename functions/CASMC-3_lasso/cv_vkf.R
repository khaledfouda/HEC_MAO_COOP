#------------------------------------------------------------------------------------------
CASMC3_kfold_M <-
  function(Y,
           # y_train is expected to be Incomplete
           X,
           obs_mask,
           n_folds = 3,
           # y: a final full-fit if provided. Expected to be Incomplete
           lambda.beta = .Machine$double.eps,
           learning.rate = 0.001,
           beta.iter.max = 20,
           # provide this if you need L2 regularization.
           error_function = error_metric$rmse,
           # tuning parameters for lambda
           lambda.factor = 1 / 4,
           lambda.init = NULL,
           n.lambda = 20,
           # tuning parameters for J
           rank.init = 2,
           rank.limit = 30,
           rank.step = 2,
           pct = 0.98,
           # laplacian parameters
           lambda.a = 0,
           S.a = NULL,
           lambda.b = 0,
           S.b = NULL,
           # stopping criteria
           early.stopping = 1,
           thresh = 1e-6,
           maxit = 100,
           # trace parameters
           trace = FALSE,
           print.best = TRUE,
           quiet = FALSE,
           # initial values.
           warm = NULL,
           # seed
           seed = NULL) {
    if (!is.null(seed))
      set.seed(seed)
    
    # prepare the sequence of lambda (nuclear regularization hyperparameter)
    if (is.null(lambda.init))
      lambda.init <-
        lambda0.cov_splr(Y, reduced_hat_decomp.H(X)) * lambda.factor
    lamseq <- seq(from = lambda.init,
                  to = .Machine$double.eps,
                  length = n.lambda)
    #-----------------------------------------------------------------------
    # prepare the folds
    folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, obs_mask)
    fold_data <- lapply(1:n_folds, function(i) {
      valid_mask = folds[[i]]
      valid_ind = valid_mask == 0 & obs_mask == 1
      
      Y_train = Y * valid_mask
      Y_train[valid_mask == 0] = NA
      Y_train = as(Y_train, "Incomplete")
      
      Y_valid = Y[valid_ind]
      
      valid_mask[valid_ind] = 1
      valid_mask[!valid_ind] = NA
      valid_mask <- as(valid_mask, "Incomplete")
      virow = valid_mask@i
      vpcol = valid_mask@p
      valid_mask <- NULL
      
      list(
        Y_train = Y_train,
        Y_valid = Y_valid,
        virow = virow,
        vpcol = vpcol
      )
    })
    #------------------------------------------------
    #-----------------------------------------------------------------------
    rank.max <- rank.init
    best_fit <- list(error = Inf)
    counter <- 0
    #---------------------------------------------------------------------
    fit_out <- list()
    fold_err <- rep(0, n_folds)
    for (i in seq(along = lamseq)) {
      err <- rank <- niter <- 0
      for (fold in 1:n_folds) {
        data = fold_data[[fold]]
        fit_out[[fold]] <-
          CASMC3_fit(
            y = data$Y_train,
            X = X,
            J = rank.max,
            lambda.M = lamseq[i],
            learning.rate = learning.rate,
            beta.iter.max = beta.iter.max,
            lambda.beta = lambda.beta,
            lambda.a = lambda.a,
            S.a = S.a,
            lambda.b = lambda.b,
            S.b = S.b,
            warm.start = warm,
            trace.it = F,
            thresh = thresh,
            maxit = maxit
          )
        # warm <- fit_out[[fold]]
        #--------------------------------------------------------------
        # predicting validation set and xbetas for next fit:
        XbetaValid = suvC(X, t(fit_out[[fold]]$beta), data$virow, data$vpcol)
        MValid = suvC(fit_out[[fold]]$u,
                      t(fit_out[[fold]]$d * t(fit_out[[fold]]$v)),
                      data$virow, data$vpcol)
        #--------------------------------------------
        fold_err[fold] <-
          error_function(MValid + XbetaValid, data$Y_valid)
        err = err + fold_err[fold]
        rank <- max(rank, sum(round(fit_out[[fold]]$d, 4) > 0))
        #rank <- max(rank, length(fit_out[[fold]]$d))
        # newly added, to be removed later
        # var_explained = fit_out[[fold]]$d ^ 2 / sum(fit_out[[fold]]$d ^ 2)
        # cum_var = cumsum(var_explained)
        # rank  <- rank + which(cum_var >= pct)[1]
        niter <- niter + fit_out[[fold]]$n_iter
      }
      err <- err / n_folds
      # rank <- as.integer(rank / n_folds)
      niter <- as.integer(niter / n_folds)
      #---------------------------------------------------------------------
      for (fold in 1:n_folds) {
        Jf = length(fit_out[[fold]]$d)
        if (Jf < rank) {
          to_fill = (Jf + 1):rank
          fit_out[[fold]]$d[to_fill] <- 0
          fit_out[[fold]]$u <- cbind(fit_out[[fold]]$u,
                                     matrix(0, nrow(fit_out[[fold]]$u), length(to_fill)))
          fit_out[[fold]]$v <- cbind(fit_out[[fold]]$v,
                                     matrix(0, nrow(fit_out[[fold]]$v), length(to_fill)))
        } else if (Jf > rank) {
          fit_out[[fold]]$d <- fit_out[[fold]]$d[1:rank]
          fit_out[[fold]]$u <- fit_out[[fold]]$u[, 1:rank]
          fit_out[[fold]]$v <- fit_out[[fold]]$v[, 1:rank]
        }
      }
      
      
      weights = 1 / fold_err #
      weights = weights / sum(weights)
      warm <- list(
        beta = weights[1] * fit_out[[1]]$beta,
        u = weights[1] * fit_out[[1]]$u,
        d = weights[1] * fit_out[[1]]$d,
        v = weights[1] * fit_out[[1]]$v
      )
      
      #M_avg <- weights[1] * unsvd(fit_out[[1]])
      if (n_folds > 1) {
        for (fold in 2:n_folds) {
          #M_avg <- M_avg + weights[fold] * unsvd(fit_out[[fold]])
          warm$beta <-
            warm$beta + weights[fold] * fit_out[[fold]]$beta
          warm$beta[fit_out[[fold]]$beta == 0] <- 0
          warm$u <-
            warm$u + weights[fold] * fit_out[[fold]]$u
          warm$d <-
            warm$d + weights[fold] * fit_out[[fold]]$d
          warm$v <-
            warm$v + weights[fold] * fit_out[[fold]]$v
          
        }
        #M_avg <- M_avg / n_folds
        #warm$beta <- warm$beta / n_folds
      }
      # warm$u <- svd(warm$u)$u
      # warm$v <- svd(warm$v)$v
      var_explained = warm$d ^ 2 / sum(warm$d ^ 2)
      cum_var = cumsum(var_explained)
      rank  <- rank + which(cum_var >= pct)[1]
      
      # Msvd <- fit_out[[n_folds]]#svd(M_avg)
      # warm$u <- Msvd$u
      # warm$d <- Msvd$d
      # warm$v <- Msvd$v
      warm$n_iter <- niter
      warm$J <- rank
      warm$lambda.M <- lamseq[i]
      #--------------------------------------------------------------------
      if (trace == TRUE)
        print(sprintf(
          paste0(
            "%2d lambda.M = %.3f, rank.max = %.0f  ==>",
            " rank = %.0f, error = %.5f, niter/fit = %.0f [Beta(lambda=%.3f)]"
          ),
          i,
          lamseq[i],
          rank.max,
          rank,
          err,
          niter,
          lambda.beta
        ))
      #-------------------------
      # register best fir
      if (err < best_fit$error) {
        best_fit$error = err
        best_fit$rank_M = rank
        best_fit$lambda.M = lamseq[i]
        best_fit$rank.max = rank.max
        best_fit$fit = warm
        best_fit$iter = i
        counter = 0
      } else
        counter = counter + 1
      if (counter >= early.stopping) {
        if (trace)
          print(
            sprintf(
              "Early stopping. Reached Peak point. Performance didn't improve for the last %.0f iterations.",
              counter
            )
          )
        break
      }
      # compute rank.max for next iteration
      rank.max <- min(rank + rank.step, rank.limit)
    }
    # fit one last time full model, if the train/valid is provided
    Y[Y == 0] = NA
    Y <- as(Y, "Incomplete")
    best_fit$fit <-
      CASMC3_fit(
        y = Y,
        X = X,
        J = best_fit$rank.max,
        lambda.M = best_fit$lambda.M,
        learning.rate = learning.rate,
        beta.iter.max = beta.iter.max,
        lambda.beta = lambda.beta,
        lambda.a = lambda.a,
        S.a = S.a,
        lambda.b = lambda.b,
        S.b = S.b,
        warm.start = best_fit$fit,
        trace.it = F,
        thresh = thresh,
        maxit = maxit
      )
    
    
    best_fit$lambda.beta = lambda.beta
    best_fit$lambda.a = lambda.a
    best_fit$lambda.b = lambda.b
    best_fit$learning.rate = learning.rate
    return(best_fit)
  }
#---------------------------------------------------------------------------------
CASMC3_kfold_M_2 <-
  function(Y,
           # y_train is expected to be Incomplete
           X,
           obs_mask,
           n_folds = 3,
           # y: a final full-fit if provided. Expected to be Incomplete
           lambda.beta = .Machine$double.eps,
           learning.rate = 0.001,
           beta.iter.max = 20,
           # provide this if you need L2 regularization.
           error_function = error_metric$rmse,
           # tuning parameters for lambda
           lambda.factor = 1 / 4,
           lambda.init = NULL,
           n.lambda = 20,
           # tuning parameters for J
           rank.init = 2,
           rank.limit = 30,
           rank.step = 2,
           pct = 0.98,
           # laplacian parameters
           lambda.a = 0,
           S.a = NULL,
           lambda.b = 0,
           S.b = NULL,
           # stopping criteria
           early.stopping = 1,
           thresh = 1e-6,
           maxit = 100,
           # trace parameters
           trace = FALSE,
           print.best = TRUE,
           quiet = FALSE,
           # initial values.
           warm = NULL,
           # seed
           seed = NULL) {
    if (!is.null(seed))
      set.seed(seed)
    
    # prepare the sequence of lambda (nuclear regularization hyperparameter)
    if (is.null(lambda.init))
      lambda.init <-
        lambda0.cov_splr(Y, reduced_hat_decomp.H(X)) * lambda.factor
    lamseq <- seq(from = lambda.init,
                  to = .Machine$double.eps,
                  length = n.lambda)
    #-----------------------------------------------------------------------
    # prepare the folds
    folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, obs_mask)
    fold_data <- lapply(1:n_folds, function(i) {
      valid_mask = folds[[i]]
      valid_ind = valid_mask == 0 & obs_mask == 1
      
      Y_train = Y * valid_mask
      Y_train[valid_mask == 0] = NA
      Y_train = as(Y_train, "Incomplete")
      
      Y_valid = Y[valid_ind]
      
      valid_mask[valid_ind] = 1
      valid_mask[!valid_ind] = NA
      valid_mask <- as(valid_mask, "Incomplete")
      virow = valid_mask@i
      vpcol = valid_mask@p
      valid_mask <- NULL
      
      list(
        Y_train = Y_train,
        Y_valid = Y_valid,
        virow = virow,
        vpcol = vpcol
      )
    })
    #------------------------------------------------
    #-----------------------------------------------------------------------
    #---------------------------------------------------------------------
    fold_fit <- vector("list", n_folds)
    for (fold in 1:n_folds) {
      warm <- NULL
      rank.max <- rank.init
      best_fit <- list(error = Inf)
      counter <- 0
      err <- rank <- niter <- 0
      data = fold_data[[fold]]
      for (i in seq(along = lamseq)) {
        fiti <-
          CASMC3_fit(
            y = data$Y_train,
            X = X,
            J = rank.max,
            lambda.M = lamseq[i],
            learning.rate = learning.rate,
            beta.iter.max = beta.iter.max,
            lambda.beta = lambda.beta,
            lambda.a = lambda.a,
            S.a = S.a,
            lambda.b = lambda.b,
            S.b = S.b,
            warm.start = warm,
            trace.it = F,
            thresh = thresh,
            maxit = maxit
          )
        warm = fiti
        # warm <- fit_out[[fold]]
        #--------------------------------------------------------------
        # predicting validation set and xbetas for next fit:
        XbetaValid = suvC(X, t(fiti$beta), data$virow, data$vpcol)
        MValid = suvC(fiti$u,
                      t(fiti$d * t(fiti$v)),
                      data$virow, data$vpcol)
        #--------------------------------------------
        err <- error_function(MValid + XbetaValid, data$Y_valid)
        #rank <- sum(round(fiti$d, 4) > 0)
        # newly added, to be removed later
        var_explained = fiti$d ^ 2 / sum(fiti$d ^ 2)
        cum_var = cumsum(var_explained)
        rank  <- rank + which(cum_var >= pct)[1]
        niter <- fiti$n_iter
        #---------------------------------------------------------------------
        #--------------------------------------------------------------------
        if (trace == TRUE)
          print(sprintf(
            paste0(
              "%2d lambda.M = %.3f, rank.max = %.0f  ==>",
              " rank = %.0f, error = %.5f, niter/fit = %.0f [Beta(lambda=%.3f)]"
            ),
            i,
            lamseq[i],
            rank.max,
            rank,
            err,
            niter,
            lambda.beta
          ))
        #-------------------------
        # register best fir
        if (err < best_fit$error) {
          best_fit$error = err
          best_fit$rank_M = rank
          best_fit$lambda.M = lamseq[i]
          best_fit$rank.max = rank.max
          best_fit$fit = fiti
          best_fit$iter = i
          counter = 0
        } else
          counter = counter + 1
        if (counter >= early.stopping) {
          if (trace)
            print(
              sprintf(
                "Early stopping. Reached Peak point. Performance didn't improve for the last %.0f iterations.",
                counter
              )
            )
          break
        }
        # compute rank.max for next iteration
        rank.max <- min(rank + rank.step, rank.limit)
      }
      #---
      # move to the next fold ->
      # what's shared?
      fold_fit[[fold]] <- best_fit
      
    }
    #---------------------------------------------------------------------
    # find the fit with the lowest validation error.
    best_fold = which.min(sapply(fold_fit, function(x)
      x$error))
    
    best_fit <- list(
      error = median(sapply(fold_fit, function(x)
        x$error)),
      rank_M = median(sapply(fold_fit, function(x)
        sum(x$fit$d > 0))) |> round(),
      lambda.M = median(sapply(fold_fit, function(x)
        x$lambda.M)),
      rank.max = median(sapply(fold_fit, function(x)
        x$rank.max)) |> round()
    )
    beta_3d = lapply(fold_fit, function(x)
      x$fit$beta) |> simplify2array()
    
    # alt beta method:
    # epsilon <- 1e-5
    # 
    # # Compute the mean of the elements, considering epsilon
    # beta <- apply(beta_3d, c(1, 2), function(x) {
    #   non_zeros <- x[x != 0]
    #   if (length(non_zeros) == 0) {
    #     return(0)
    #   } else {
    #     return(mean(non_zeros) * (length(non_zeros) / length(x)) + 
    #              epsilon * (length(x) - length(non_zeros)) / length(x))
    #   }
    # })
    # 
    # threshold <- 1e-4
    # beta[abs(beta) < threshold] <- 0
    #---
    zero_mask <- apply(beta_3d == 0, c(1, 2), any)
    beta <- apply(beta_3d, c(1, 2), function(x)
      mean(x[x != 0]))
    beta[zero_mask] <- 0
    
    # pad them:
    A_comp  <- B_comp <- M_comp <- 0
    print(best_fit$rank_M)
    
    for (fit in fold_fit) {
      u <- fit$fit$u
      d <- fit$fit$d
      v <- fit$fit$v
      # Jf = length(d)
      # if (Jf < best_fit$rank_M) {
      #   to_fill = (Jf + 1):best_fit$rank_M
      #   d[to_fill] <- 0
      #   u <- cbind(u, matrix(0, nrow(u), length(to_fill)))
      #   v <- cbind(v, matrix(0, nrow(v), length(to_fill)))
      # } else if (Jf > best_fit$rank_M) {
      #   d <- d[1:best_fit$rank_M]
      #   u <- u[, 1:best_fit$rank_M]
      #   v <- v[, 1:best_fit$rank_M]
      # }
      
      M_comp <- M_comp + u %*% (d * t(v))
      # A_comp <- A_comp + u %*% diag(sqrt(d))
      # B_comp <- B_comp + v %*% diag(sqrt(d))
    }
    # A_comp <- A_comp / n_folds
    # B_comp <- B_comp / n_folds
    # Msvd <- svd(A_comp %*% t(B_comp))
    # Msvd <- fold_fit[[best_fold]]$fit
    Msvd <- svd_trunc_simple(M_comp / n_folds, best_fit$rank_M)
    
    #print(length(fit$fit$d))
    #U_comp <- U_comp + fit$fit$u %*% diag(fit$fit$u)
    
    
    
    #beta = Reduce("+", lapply(fold_fit, function(x) x$fit$beta)) / n_folds
    #Msvd = svd(Reduce("+", lapply(fold_fit, function(x)
    #unsvd(x$fit))) / n_folds)
    best_fit$fit <- list(
      beta = beta,
      u = Msvd$u,
      d = Msvd$d,
      v = Msvd$v,
      J = best_fit$rank_M,
      lambda.M = best_fit$lambda.M,
      n_iter = median(sapply(fold_fit, function(x)
        x$fit$n_iter) |> round())
    )
    
    #best_fit$beta = fold_fit[[1]]$fit$beta
    # for(fold in 2:n_folds){
    #   best_fit$error <- best_fit$error + fold_fit[[fold]]$error
    #   best_fit$rank_M <- best_fit$rank_M + fold_fit[[fold]]$rank_M
    #   best_fit$lambda.M <- max(best_fit$lambda.M, fold_fit[[fold]]$lambda.M)
    #   best_fit$rank.max <- min(best_fit$rank.max, fold_fit[[fold]]$rank.max)
    #   #best_fit$fit$beta <- best_fit$fit$beta + fold_fit[[fold]]$fit$beta
    # }
    # best_fit$error <- best_fit$error / n_folds
    # best_fit$rank_M <- round(best_fit$rank_M / n_folds)
    # best_fit$lambda.M <- best_fit$lambda.M / n_folds
    # best_fit$rank.max <- round(best_fit$rank.max / n_folds)
    #best_fit$fit$beta = best_fit$fit$beta / n_folds
    #---------------------------------------------------------------------
    # fit one last time full model, if the train/valid is provided
    # Y[Y == 0] = NA
    # Y <- as(Y, "Incomplete")
    # best_fit$fit <-
    #   CASMC3_fit(
    #     y = Y,
    #     X = X,
    #     J = best_fit$rank.max,
    #     lambda.M = best_fit$lambda.M,
    #     learning.rate = learning.rate,
    #     beta.iter.max = beta.iter.max,
    #     lambda.beta = lambda.beta,
    #     lambda.a = lambda.a,
    #     S.a = S.a,
    #     lambda.b = lambda.b,
    #     S.b = S.b,
    #     warm.start = best_fit$fit,
    #     trace.it = F,
    #     thresh = thresh,
    #     maxit = maxit
    #   )
    
    
    best_fit$lambda.beta = lambda.beta
    best_fit$lambda.a = lambda.a
    best_fit$lambda.b = lambda.b
    best_fit$learning.rate = learning.rate
    #best_fit$rank_M = best_fit$fit$J
    return(best_fit)
  }

#---------------------------------------------------------------------------------



CASMC3_kfold <-
  function(Y,
           # Y is expected to be Incomplete
           obs_mask,
           X,
           n_folds = 3,
           # y: a final full-fit if provided. Expected to be Incomplete
           error_function = error_metric$rmse,
           # tuning parameters for lambda
           lambda.factor = 1 / 4,
           lambda.init = NULL,
           n.lambda = 20,
           # tuning parameters for J
           rank.M.init = 2,
           rank.M.limit = 30,
           rank.step = 2,
           pct = 0.98,
           # laplacian parameters
           lambda.a = 0,
           S.a = NULL,
           lambda.b = 0,
           S.b = NULL,
           # stopping criteria
           early.stopping = 1,
           thresh = 1e-6,
           maxit = 300,
           # trace parameters
           trace = 0,
           print.best = TRUE,
           quiet = FALSE,
           # initial values.
           warm = NULL,
           # L1 parameters
           lambda.beta.grid = "default1",
           learning.rate = .001,
           beta.iter.max = 20,
           track = FALSE,
           max_cores = 8,
           # seed
           seed = NULL) {
    if (identical(lambda.beta.grid, "default1"))
      lambda.beta.grid = sqrt((ncol(Y) * ncol(X)) / (nrow(Y))) *
        seq(10, .Machine$double.eps, length.out = 20)
    if (identical(lambda.beta.grid, "default2"))
      lambda.beta.grid = seq(propack.svd(naive_fit(Y, X, TRUE), 1)$d / 4,
                             .Machine$double.eps,
                             length.out = 20)
    #--------------------------------------------------------------------------
    
    num_cores = length(lambda.beta.grid)
    if (length(lambda.beta.grid) > max_cores)
      num_cores <-
        min(max_cores, ceiling(length(lambda.beta.grid) / 2))
    print(paste("Running on", num_cores, "cores."))
    
    best_fit <- list(error = Inf)
    results <- mclapply(lambda.beta.grid, function(lambda.beta) {
      fiti = tryCatch({
        CASMC3_kfold_M_2(
          Y = Y,
          X = X,
          n_folds = n_folds,
          obs_mask = obs_mask,
          learning.rate = learning.rate,
          lambda.beta = lambda.beta,
          beta.iter.max = beta.iter.max,
          error_function = error_function,
          lambda.factor = lambda.factor,
          lambda.init = lambda.init,
          n.lambda = n.lambda,
          trace = ifelse(trace == 2, TRUE, FALSE),
          print.best = print.best,
          thresh = thresh,
          maxit = maxit,
          rank.init = rank.M.init,
          rank.limit = rank.M.limit,
          rank.step = rank.step,
          pct = pct,
          warm = NULL,
          lambda.a = lambda.a,
          S.a = S.a,
          lambda.b = lambda.b,
          S.b = S.b,
          quiet = quiet,
          seed = seed
        )
      }, error = function(e)
        list(
          error_message = e,
          error = 999999,
          lambda.beta = lambda.beta
        ))
      fiti
    }, mc.cores = num_cores, mc.cleanup = TRUE)
    
    
    # showing errors, if any,
    sapply(results, function(x) {
      if (!is.null(x$error_message))
        print(
          paste(
            "Error encountered at lambda.beta = ",
            round(x$lambda.beta, 3),
            "with the following error message: ",
            x$error_message
          )
        )
    })
    
    # extract best fit
    best_fit <-
      results[[which.min(sapply(results, function(x)
        x$error))]]
    
    # print all fit output:
    if (trace > 0) {
      sapply(results, function(x)
        print(
          sprintf(
            paste0(
              "<< lambda.beta = %.3f, error = %.5f, niter/fit = %.0f, M = [%.0f,%.3f] >> "
            ),
            x$lambda.beta,
            x$error,
            x$fit$n_iter,
            x$fit$J,
            x$fit$lambda.M
          )
        ))
    }
    
    # print best if
    if (print.best)
      print(
        sprintf(
          paste0(
            "<< Best fit >> lambda.beta = %.3f, error = %.5f, niter/fit = %.0f, M = [%.0f,%.3f] > "
          ),
          best_fit$lambda.beta,
          best_fit$error,
          best_fit$fit$n_iter,
          best_fit$fit$J,
          best_fit$fit$lambda.M
        )
      )
    #-------------------------
    
    best_fit$hparams = data.frame(
      lambda.beta = best_fit$lambda.beta,
      lambda.M = best_fit$fit$lambda.M,
      learning.rate = learning.rate,
      rank.M = best_fit$fit$J
    )
    
    return(best_fit)
    
    
  }
