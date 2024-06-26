Mao_Sim_Wrapper <-
  function(dat,
           lambda.1_grid = seq(0, 1, length = 20),
           lambda.2_grid = seq(0.9, 0.1, length = 20),
           alpha_grid = c(1),
           ncores = 1,
           # keep it > 1
           n_folds = 5,
           weight_function = Mao_weights$uniform,
           ...) {
    start_time = Sys.time()
    fiti <- Mao.cv(
      Y = dat$Y,
      X = dat$X,
      W = dat$W,
      n_folds = n_folds,
      lambda.1_grid = lambda.1_grid,
      lambda.2_grid = lambda.2_grid,
      alpha_grid = alpha_grid,
      seed = 2023,
      numCores = ncores,
      n1n2_optimized = TRUE,
      test_error = error_metric$rmse,
      theta_estimator = weight_function,
      sequential = FALSE
    )
    
    fit. <- fiti$fit
    results = list(model = "Mao")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$lambda.M = fiti$best_parameters$lambda.2
    results$lambda.beta = fiti$best_parameters$lambda.1
    results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W ==
                                                                        0])
    results$error.all = test_error(fit.$estimates, dat$O)
    results$error.M = test_error(fit.$M, dat$M)
    results$error.beta = test_error(fit.$beta, dat$beta)
    results$rank_M = qr(fit.$M)$rank
    results$rank_beta = qr(fit.$beta)$rank
    results$sparse_in_sparse = sum(dat$beta == 0 & fit.$beta == 0) /
      (sum(dat$beta == 0) +  1e-17)
    results$sparse_in_nonsparse = sum(dat$beta != 0 &
                                        fit.$beta == 0) /
      (sum(dat$beta != 0) +  1e-17)
    results
  }

SImpute_Sim_Wrapper <- function(dat, ...) {
  start_time = Sys.time()
  fit. <- simpute.cv(
    Y_train = as.matrix(dat$fit_data$train),
    y_valid = dat$fit_data$valid,
    W_valid = dat$fit_data$W_valid,
    y = dat$Y,
    n.lambda = 20,
    trace = FALSE,
    print.best = FALSE,
    tol = 5,
    thresh = 1e-6,
    rank.init = 2,
    rank.limit = 30,
    rank.step = 2,
    maxit = 600,
    seed = NULL
  )
  results = list(model = "SoftImpute")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$lambda.beta = NA
  results$lambda.M = fit.$lambda
  results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W ==
                                                                      0])
  results$error.all = test_error(fit.$estimates, dat$O)
  results$error.M = NA
  results$error.beta = NA
  results$rank_M = fit.$rank_M
  results$rank_beta = NA
  results$sparse_in_sparse = NA
  results$sparse_in_nonsparse = NA
  results
}


#----------------------------------------------------------------------------------
CASMC_1_Sim_Wrapper <-
  function(dat,
           max_cores = 20,
           maxit = 300,
           ...) {
    start_time = Sys.time()
    
    fiti <- CASMC1_cv(
      y_train = dat$fit_data$train,
      X = dat$X,
      y_valid = dat$fit_data$valid,
      W_valid = dat$fit_data$W_valid,
      y = dat$fit_data$Y,
      error_function = error_metric$rmse,
      lambda.factor = 1 / 4,
      lambda.init = NULL,
      n.lambda = 20,
      rank.init = 2,
      rank.limit = 30,
      rank.step = 2,
      pct = 0.98,
      lambda.a = 0,
      S.a = NULL,
      lambda.b = 0,
      S.b = NULL,
      early.stopping = 1,
      thresh = 1e-6,
      maxit = maxit,
      trace = F,
      print.best = TRUE,
      quiet = FALSE,
      warm = NULL,
      r_min = 0,
      track = F,
      max_cores = max_cores,
      seed = NULL
    )
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model="CASMC-1")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$lambda.M = fit.$lambda
    results$lambda.beta = NA
    results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W == 0])
    results$error.all = test_error(fit.$estimates, dat$O)
    results$error.M = test_error(fit.$M, dat$M)
    results$error.beta = test_error(fit.$beta, dat$beta)
    results$rank_M = sum(fit.$d > 0)
    results$rank_beta = qr(fit.$beta)$rank
    results$sparse_in_sparse = sum(dat$beta == 0 & fit.$beta == 0) /
      (sum(dat$beta == 0) +  1e-17)
    results$sparse_in_nonsparse = sum(dat$beta != 0 &
                                        fit.$beta == 0) /
      (sum(dat$beta != 0) +  1e-17)
    results
  }

#--------------------------------------------------------------------------------------
CASMC_0_Sim_Wrapper <-
  function(dat,
           max_cores = 20,
           maxit = 300,
           ...) {
    start_time = Sys.time()
    
    fiti <- CASMC0_cv(
      y_train = dat$fit_data$train,
      X = dat$X,
      y_valid = dat$fit_data$valid,
      W_valid = dat$fit_data$W_valid,
      y = dat$fit_data$Y,
      error_function = error_metric$rmse,
      lambda.factor = 1 / 4,
      lambda.init = NULL,
      n.lambda = 20,
      rank.init = 2,
      rank.limit = 30,
      rank.step = 2,
      pct = 0.98,
      lambda.a = 0,
      S.a = NULL,
      lambda.b = 0,
      S.b = NULL,
      early.stopping = 1,
      thresh = 1e-6,
      maxit = maxit,
      trace = FALSE,
      print.best = F,
      quiet = FALSE,
      warm = NULL,
      lambda.beta.grid = "default",
      track = F,
      max_cores = max_cores,
      seed = NULL
    )
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC-0")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$lambda.M = fit.$lambda
    results$lambda.beta = fiti$lambda.beta
    results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W == 0])
    results$error.all = test_error(fit.$estimates, dat$O)
    results$error.M = test_error(fit.$M, dat$M)
    results$error.beta = test_error(fit.$beta, dat$beta)
    results$rank_M = sum(fit.$d > 0)
    results$rank_beta = qr(fit.$beta)$rank
    results$sparse_in_sparse = sum(dat$beta == 0 & fit.$beta == 0) /
      (sum(dat$beta == 0) +  1e-17)
    results$sparse_in_nonsparse = sum(dat$beta != 0 &
                                        fit.$beta == 0) /
      (sum(dat$beta != 0) +  1e-17)
    results
  }
#-------


CASMC_2_Sim_Wrapper <-
  function(dat,
           max_cores = 20,
           maxit = 300,
           ...) {
    start_time = Sys.time()
    
    fiti <- CASMC2_cv(
      y_train = dat$fit_data$train,
      X = dat$X,
      y_valid = dat$fit_data$valid,
      W_valid = dat$fit_data$W_valid,
      y = dat$fit_data$Y,
      error_function = error_metric$rmse,
      warm = NULL,
      quiet = F,
      rank.beta.init = 1,
      lambda.beta.grid = "default1",
      max_cores = max_cores,
      seed = NULL,
    )
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$beta = fit.$ub %*% (fit.$db^2) %*% t(fit.$vb)
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC-2")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$lambda.M = fiti$hparams$lambda.M
    results$lambda.beta = fiti$hparams$lambda.beta
    results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W == 0])
    results$error.all = test_error(fit.$estimates, dat$O)
    results$error.M = test_error(fit.$M, dat$M)
    results$error.beta = test_error(fit.$beta, dat$beta)
    results$rank_M = sum(fit.$d > 0)
    results$rank_beta = qr(fit.$beta)$rank
    results$sparse_in_sparse = sum(dat$beta == 0 & fit.$beta == 0) /
      (sum(dat$beta == 0) +  1e-17)
    results$sparse_in_nonsparse = sum(dat$beta != 0 &
                                        fit.$beta == 0) /
      (sum(dat$beta != 0) +  1e-17)
    results
  }
#----------------------------------------------------
CASMC_3a_Sim_Wrapper <-
  function(dat,
           max_cores = 20,
           maxit = 300,
           ...) {
    start_time = Sys.time()
    learning_rate = 1 / sqrt(sum((t(dat$X) %*% dat$X)^2))
    fiti <- CASMC3_cv_beta(
      y_train = dat$fit_data$train,
      X = dat$X,
      y_valid = dat$fit_data$valid,
      W_valid = dat$fit_data$W_valid,
      y = dat$fit_data$Y,
      trace = 0,
      print.best = T,
      warm = NULL,
      quiet = F, learning.rate = learning_rate,
      early.stopping = 1,
      lambda.beta.grid = seq(0,10,length.out=20),
      max_cores = max_cores
    ) 
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC-3a")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$lambda.M = fiti$hparams$lambda.M
    results$lambda.beta = fiti$hparams$lambda.beta
    results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W == 0])
    results$error.all = test_error(fit.$estimates, dat$O)
    results$error.M = test_error(fit.$M, dat$M)
    results$error.beta = test_error(fit.$beta, dat$beta)
    results$rank_M = sum(fit.$d > 0)
    results$rank_beta = qr(fit.$beta)$rank
    results$sparse_in_sparse = sum(dat$beta == 0 & fit.$beta == 0) /
      (sum(dat$beta == 0) +  1e-17)
    results$sparse_in_nonsparse = sum(dat$beta != 0 &
                                        fit.$beta == 0) /
      (sum(dat$beta != 0) +  1e-17)
    results
  }
#------
CASMC_3b_Sim_Wrapper <-
  function(dat,
           max_cores = 20,
           maxit = 300,
           ...) {
    start_time = Sys.time()
    learning_rate = 1 / sqrt(sum((t(dat$X) %*% dat$X)^2))
    fiti <- CASMC3_kfold(
      Y = dat$Y,
      X = dat$X,
      obs_mask = dat$W, n_folds = 10,
      trace = 0,
      print.best = T,
      warm = NULL,
      quiet = F, learning.rate = learning_rate,
      early.stopping = 1,
      lambda.beta.grid = seq(0,10,length.out=20),
      max_cores = max_cores
    ) 
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC-3b")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$lambda.M = fiti$hparams$lambda.M
    results$lambda.beta = fiti$hparams$lambda.beta
    results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W == 0])
    results$error.all = test_error(fit.$estimates, dat$O)
    results$error.M = test_error(fit.$M, dat$M)
    results$error.beta = test_error(fit.$beta, dat$beta)
    results$rank_M = sum(fit.$d > 0)
    results$rank_beta = qr(fit.$beta)$rank
    results$sparse_in_sparse = sum(dat$beta == 0 & fit.$beta == 0) /
      (sum(dat$beta == 0) +  1e-17)
    results$sparse_in_nonsparse = sum(dat$beta != 0 &
                                        fit.$beta == 0) /
      (sum(dat$beta != 0) +  1e-17)
    results
  }
#------





#----------------------------------------------
Naive_Sim_Wrapper <- function(dat, ...) {
  start_time = Sys.time()
  fiti <- naive_fit(dat$Y, dat$X)
  results = list(model = "Naive")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$lambda.M = NA
  results$lambda.beta = NA
  results$error.test = test_error(fiti$estimates[dat$W == 0], dat$O[dat$W == 0])
  results$error.all = test_error(fiti$estimates, dat$O)
  results$error.M = test_error(fiti$M, dat$M)
  results$error.beta = test_error(fiti$beta, dat$beta)
  results$rank_M = qr(fiti$M)$rank
  results$rank_beta = qr(fiti$beta)$rank
  results$sparse_in_sparse = sum(dat$beta == 0 & fiti$beta == 0) /
    (sum(dat$beta == 0) +  1e-17)
  results$sparse_in_nonsparse = sum(dat$beta != 0 &
                                      fiti$beta == 0) /
    (sum(dat$beta != 0) +  1e-17)
  results
}
