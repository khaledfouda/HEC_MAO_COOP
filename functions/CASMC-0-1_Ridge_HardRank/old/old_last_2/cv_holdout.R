CASMC_cv_holdout_with_reg <-
   function(y_train,
            X_r,
            y_valid,
            W_valid,
            y = NULL,
            lambda.beta.grid = "default",
            error_function = error_metric$rmse,
            lambda.factor = 1 / 4,
            lambda.init = NULL,
            n.lambda = 20,
            trace = FALSE,
            print.best = TRUE,
            early.stopping = 1,
            thresh = 1e-6,
            maxit = 100,
            rank.init = 2,
            rank.limit = 30,
            rank.step = 2,
            lambda.a = 0,
            S.a = NULL,
            lambda.b = 0,
            S.b = NULL,
            warm = NULL,
            track_beta = FALSE,
            max_cores = 10,
            pct = 0.98,
            quiet = FALSE,
            seed = NULL) {
      
      if(identical(lambda.beta.grid, "default"))
         lambda.beta.grid = sqrt((ncol(y_train)*ncol(X_r$X))/(nrow(y_train))) *  seq(0, 10, length.out = 20)
      
      num_cores = length(lambda.beta.grid)
      if(length(lambda.beta.grid) > max_cores)
         num_cores <- min(max_cores, ceiling(length(lambda.beta.grid)/2))
      print(paste("Running on", num_cores, "cores."))
      
      
      best_score = Inf
      best_fit = NULL
      results <- mclapply(lambda.beta.grid, function(lambda.beta) {
         
         Xterms = GetXterms(X_r$X, lambda.beta)
         fiti = CASMC_cv_holdout(
            y_train = y_train,
            X_r = X_r,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            Xterms = Xterms,
            r = NULL,
            error_function = error_function,
            lambda.factor = lambda.factor,
            lambda.init = lambda.init,
            n.lambda = n.lambda,
            trace = trace,
            print.best = print.best,
            early.stopping = early.stopping,
            thresh = thresh,
            maxit = maxit,
            rank.init = rank.init,
            rank.limit = rank.limit,
            rank.step = rank.step,
            pct = pct,
            warm = NULL,
            quiet = quiet,
            seed = seed,
            lambda.a = lambda.a,
            S.a = S.a,
            lambda.b = lambda.b,
            S.b = S.b
         )
         fiti$lambda.beta = lambda.beta
         fiti
      }, mc.cores = num_cores)
      
      
      best_fit <-
         results[[which.min(sapply(results, function(x)
            x$error))]]
      
      if (track_beta) {
         sapply(results, function(x)
            print(paste(x$lambda.beta, "-", x$error)))
      }
      if(print.best)
         print(paste("Best fit: lambda_beta = ",best_fit$lambda.beta, " - Validation Error: ", best_fit$error))
      
      #---------------
      # # beta rank reduction:
      # warm = best_fit$fit
      # warm$Beta = svd((as.matrix(warm$beta)))
      # rank_fit <- CASMC_cv_holdout_with_r(
      #    y_train = y_train,
      #    X_r = X_r,
      #    y_valid = y_valid,
      #    W_valid = W_valid,
      #    lambda.beta = best_fit$lambda.beta,
      #    lambda.M = best_fit$lambda,
      #    J = best_fit$rank.max,
      #    y = y,
      #    r_min = 0,
      #    r_max = qr(best_fit$fit$beta)$rank,
      #    error_function = error_function,
      #    trace = trace,
      #    print.best = print.best,
      #    thresh = thresh,
      #    maxit = maxit,
      #    lambda.a = lambda.a,
      #    S.a = S.a,
      #    lambda.b = lambda.b,
      #    S.b = S.b,
      #    warm = warm,
      #    track_r = track_beta,
      #    pct = pct,
      #    quiet = quiet,
      #    seed = seed
      # )
      # 
      # 
      # #rank_fit$lambda.beta = best_fit$lambda.beta
      # best_fit$rank_fit = rank_fit
      return(best_fit)
      
      
   }
#--------------------------------------------------------------------------------------
CASMC_cv_holdout_with_reg2 <-
   function(y_train,
            X_r,
            y_valid,
            W_valid,
            y = NULL,
            lambda.beta.grid = seq(0,2,length.out=10),
            error_function = error_metric$rmse,
            lambda.factor = 1 / 4,
            lambda.init = NULL,
            n.lambda = 20,
            trace = FALSE,
            print.best = TRUE,
            early.stopping = 1,
            thresh = 1e-6,
            maxit = 100,
            rank.init = 2,
            rank.limit = 30,
            rank.step = 2,
            lambda.a = 0,
            S.a = NULL,
            lambda.b = 0,
            S.b = NULL,
            warm = NULL,
            track_beta = FALSE,
            max_cores = 10,
            pct = 0.98,
            quiet = FALSE,
            seed = NULL) {
      num_cores = length(lambda.beta.grid)
      if(length(lambda.beta.grid) > max_cores)
         num_cores <- min(max_cores, ceiling(length(lambda.beta.grid)/2))
      print(paste("Running on", num_cores, "cores."))
      
      
      best_score = Inf
      best_fit = NULL
      # Initialize a list to store the results
      results <- list()
      fiti <- list(fit=NULL)
      # Loop through each lambda.beta value
      for (i in seq_along(lambda.beta.grid)) {
         lambda.beta <- lambda.beta.grid[i]
         
         # Perform the operations that were inside the mclapply function
         Xterms <- GetXterms(X_r$X, lambda.beta)
         fiti <- CASMC_cv_holdout(
            y_train = y_train,
            X_r = X_r,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            Xterms = Xterms,
            r = NULL,
            error_function = error_function,
            lambda.factor = lambda.factor,
            lambda.init = lambda.init,
            n.lambda = n.lambda,
            trace = trace,
            print.best = print.best,
            early.stopping = early.stopping,
            thresh = thresh,
            maxit = maxit,
            rank.init = rank.init,
            rank.limit = rank.limit,
            rank.step = rank.step,
            pct = pct,
            warm = fiti$fit,
            quiet = quiet,
            seed = seed,
            lambda.a = lambda.a,
            S.a = S.a,
            lambda.b = lambda.b,
            S.b = S.b
         )
         
         # Add lambda.beta to the result for traceability
         fiti$lambda.beta <- lambda.beta
         
         # Store the fit object in the results list
         results[[i]] <- fiti
      }
      
      # Identify the best fit based on minimum error
      best_fit <- results[[which.min(sapply(results, function(x) x$error))]]
      
      # If tracking beta, print lambda.beta and error
      if (track_beta) {
         sapply(results, function(x) print(paste(x$lambda.beta, "-", x$error)))
      }
      
      
      return(best_fit)
      
      
   }
#------------------------------------------------------------------------------------------
CASMC_cv_holdout_rank <-
   function(y_train,
            X_r,
            y_valid,
            W_valid,
            lambda.beta,
            y = NULL,
            r_min = 0,
            r_max = X_r$rank,
            error_function = error_metric$rmse,
            lambda.factor = 1 / 4,
            lambda.init = NULL,
            n.lambda = 20,
            trace = FALSE,
            print.best = TRUE,
            early.stopping = 1,
            thresh = 1e-6,
            maxit = 100,
            rank.init = 2,
            rank.limit = 30,
            rank.step = 2,
            lambda.a = 0,
            S.a = NULL,
            lambda.b = 0,
            S.b = NULL,
            warm = NULL,
            track_r = FALSE,
            max_cores = 12,
            pct = 0.98,
            quiet = FALSE,
            seed = NULL) {
      r_seq <- (max(r_min, 0)):(min(X_r$rank, r_max))
      Xterms = GetXterms(X_r$X, lambda.beta)
      best_score = Inf
      best_fit = NULL
      num_cores = min(max_cores, length(r_seq))
      print(paste("Running on", num_cores, "cores."))
      results <- mclapply(r_seq, function(r) {
         CASMC_cv_holdout(
            y_train = y_train,
            X_r = X_r,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            Xterms = Xterms,
            r = r,
            error_function = error_function,
            lambda.factor = lambda.factor,
            lambda.init = lambda.init,
            n.lambda = n.lambda,
            trace = trace,
            print.best = print.best,
            early.stopping = early.stopping,
            thresh = thresh,
            maxit = maxit,
            rank.init = rank.init,
            rank.limit = rank.limit,
            rank.step = rank.step,
            pct = pct,
            warm = NULL,
            quiet = quiet,
            seed = seed,
            lambda.a = lambda.a,
            S.a = S.a,
            lambda.b = lambda.b,
            S.b = S.b
         )
      }, mc.cores = num_cores)
      
      
      best_fit <-
         results[[which.min(sapply(results, function(x)
            x$error))]]
      
      if (track_r) {
         sapply(results, function(x)
            print(paste(x$r, "-", x$error)))
      }
      
      return(best_fit)
      
      
   }


#------------------------------------------------------------------------------------------
CASMC_cv_holdout_with_r_old <-
   function(y_train,
            X_r,
            y_valid,
            W_valid,
            y = NULL,
            r_min = 0,
            r_max = X_r$rank,
            error_function = error_metric$rmse,
            lambda.factor = 1 / 4,
            lambda.init = NULL,
            n.lambda = 20,
            Xterms=NULL,
            trace = FALSE,
            print.best = TRUE,
            early.stopping = 1,
            thresh = 1e-6,
            maxit = 100,
            rank.init = 2,
            rank.limit = 30,
            rank.step = 2,
            lambda.a = 0,
            S.a = NULL,
            lambda.b = 0,
            S.b = NULL,
            warm = NULL,
            track_r = FALSE,
            max_cores = 12,
            pct = 0.98,
            quiet = FALSE,
            seed = NULL) {
      r_seq <- (max(r_min, 0)):(min(X_r$rank, r_max))
      if(is.null(Xterms)) Xterms = GetXterms(X_r$X)
      best_score = Inf
      best_fit = NULL
      num_cores = min(max_cores, length(r_seq))
      print(paste("Running on", num_cores, "cores."))
      results <- mclapply(r_seq, function(r) {
         CASMC_cv_holdout(
            y_train = y_train,
            X_r = X_r,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            Xterms = Xterms,
            r = r,
            error_function = error_function,
            lambda.factor = lambda.factor,
            lambda.init = lambda.init,
            n.lambda = n.lambda,
            trace = trace,
            print.best = print.best,
            early.stopping = early.stopping,
            thresh = thresh,
            maxit = maxit,
            rank.init = rank.init,
            rank.limit = rank.limit,
            rank.step = rank.step,
            pct = pct,
            warm = NULL,
            type="rank",
            quiet = quiet,
            seed = seed,
            lambda.a = lambda.a,
            S.a = S.a,
            lambda.b = lambda.b,
            S.b = S.b
         )
      }, mc.cores = num_cores)
      
      
      best_fit <-
         results[[which.min(sapply(results, function(x)
            x$error))]]
      
      if (track_r) {
         sapply(results, function(x)
            print(paste(x$r, "-", x$error)))
      }
      
      return(best_fit)
      
      
   }





#---------------------------------------------------------------------------------------------
#' Covariate-Adjusted-Sparse-Matrix-completion
#' Cross-Validation function with holdout method (training / validation)
#'
#' @param y_train A sparse matrix of class incomplete (training)
#' @param X_r A list containing SvdH and X (see fit function for more details)
#' @param y_valid A vector of the validation set.
#' @param W_valid An indicator matrix where 0 = validation and 1 = train/test
#' @param y (optional) A sparse matrix of class Incomplete of the train and validation combined. A final fit step will be
#'                            applied if y is provided.
#' @return A list of u,d,v of M, Beta, and a vector of the observed Xbeta
#' @examples
#'  CASMC_fit(y,X,J=5)
#' @export
#'
CASMC_cv_holdout <-
   function(y_train,
            X_r,
            y_valid,
            W_valid,
            r = NULL,
            y = NULL,
            Xterms = NULL,
            error_function = error_metric$rmse,
            lambda.factor = 1 / 4,
            lambda.init = NULL,
            n.lambda = 20,
            trace = FALSE,
            print.best = TRUE,
            early.stopping = 1,
            thresh = 1e-6,
            maxit = 100,
            rank.init = 2,
            rank.limit = 30,
            rank.step = 2,
            lambda.a = 0,
            S.a = NULL,
            lambda.b = 0,
            S.b = NULL,
            warm = NULL,
            type = "reg",
            pct = 0.98,
            quiet = FALSE,
            seed = NULL) {
      if (!is.null(seed))
         set.seed(seed)
      #stopifnot(type %in% c("rank", "reg"))
      if(type=="rank"){
         fit_func = CASMC_fit_rank
      }else fit_func = CASMC_fit
      # prepare the sequence of lambda (nuclear regularization hyperparameter)
      if (is.null(lambda.init))
         lambda.init <-
            lambda0.cov_splr(y_train, X_r$svdH) * lambda.factor
      lamseq <- seq(from = lambda.init,
                    to = 0,
                    length = n.lambda)
      #----------------------------------------------------
      stopifnot(inherits(y_train, "dgCMatrix"))
      # we only need the indices for validation from W_valid
      W_valid[W_valid == 1] = NA
      W_valid[W_valid == 0] =  1
      W_valid <- as(W_valid, "Incomplete")
      virow = W_valid@i
      vpcol = W_valid@p
      W_valid = NULL
      #------------------------------------------------
      if (is.null(Xterms))
         Xterms = GetXterms(X_r$X)
      #-----------------------------------------------------------------------
      rank.max <- rank.init
      best_fit <- list(error = Inf, r = r)
      counter <- 0
      #---------------------------------------------------------------------
      for (i in seq(along = lamseq)) {
         fiti <-
            fit_func(
               y = y_train,
               X = X_r$X,
               #svdH = X_r$svdH,
               Xterms = Xterms,
               r = r,
               J = rank.max,
               lambda = lamseq[i],
               warm.start = warm,
               trace.it = F,
               thresh = thresh,
               lambda.a = lambda.a,
               S.a = S.a,
               lambda.b = lambda.b,
               S.b = S.b,
               init = "naive",
               final.svd = T,
               maxit = maxit
            )
         
         #--------------------------------------------------------------
         # predicting validation set and xbetas for next fit:
         #Beta = fiti$Beta
         XbetaValid = suvC(X_r$X , t(fiti$beta), virow, vpcol)
         #XbetaValid = suvC(X_r$X %*% Beta$v, t(Beta$d * t(Beta$u)), virow, vpcol)
         MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), virow, vpcol)
         #--------------------------------------------
         err = error_function(MValid + XbetaValid, y_valid)
         rank <- sum(round(fiti$d, 4) > 0)
         # newly added, to be removed later
         var_explained = fiti$d ^ 2 / sum(fiti$d ^ 2)
         cum_var = cumsum(var_explained)
         rank <- rank2 <- which(cum_var >= pct)[1]
         #print( fiti$d)
         #rank <- rank2 <- min(rank, 2, na.rm=TRUE)
         warm <- fiti # warm start for next
         #print(paste(rank,"-",rank2))
         #---------------------------------------------------------------------
         if (trace == TRUE)
            print(sprintf(
               paste0(
                  "%2d lambda=%9.5g, rank.max = %d  ==>",
                  " rank = %d, error = %.5f, niter/fit = %d"
               ),
               i,
               lamseq[i],
               rank.max,
               rank,
               err,
               fiti$n_iter
            ))
         #-------------------------
         # register best fir
         if (err < best_fit$error) {
            best_fit$error = err
            best_fit$rank_M = rank
            best_fit$lambda = lamseq[i]
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
                     "Early stopping. Reached Peak point. Performance didn't improve for the last %d iterations.",
                     counter
                  )
               )
            break
         }
         # compute rank.max for next iteration
         rank.max <- min(rank2 + rank.step, rank.limit)
      }
      # fit one last time full model, if the train/valid is provided
      if (!is.null(y)) {
         stopifnot(inherits(y, "dgCMatrix"))
         best_fit$fit <-
            CASMC_fit(
               y = y,
               X = X_r$X,
               #svdH = X_r$svdH,
               Xterms = Xterms,
               r = r,
               J = best_fit$rank.max,
               lambda = best_fit$lambda,
               warm.start = best_fit$fit,
               thresh = thresh,
               maxit = maxit,
               trace.it = F,
               lambda.a = lambda.a,
               S.a = S.a,
               lambda.b = lambda.b,
               S.b = S.b,
               final.svd = T
            )
         
      }
      return(best_fit)
   }
#---
