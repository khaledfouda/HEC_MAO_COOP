

missingness = 0.9

compare_and_save <- function(missingness,coll=TRUE, 
                             lambda.1_grid = seq(0,3,length=20),
                             lambda.2_grid = seq(.9, 0, length=20),
                             alpha_grid = seq(0.992, 1, length=10), plot=FALSE, tofile=FALSE, graph_label="",ncores=2){
   
   setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/")
   source("./Mao/SMC_functions6.R")
   source("./code_files/Mao_import_lib.R")
   data_dir = "./Mao/saved_data/"
   stopifnot(missingness %in% c(0,0.8, 0.9))
   dim = seq(400,1000,200)
   results <- data.frame(Dim=paste0(dim,"x",dim), true_rank=NA, 
                         Mao.alpha=NA, Mao.lambda.1=NA, Mao.lambda.2=NA, Mao.error.test=NA,
                         Mao.error.all=NA, Mao.error.B=NA, Mao.error.beta=NA, Mao.rank=NA, Mao.time=NA,
                         simputeCov.time = NA,
                         simputeCov.lambda.1 = NA,
                         simputeCov.error.test = NA,
                         simputeCov.error.all = NA,
                         simputeCov.error.beta = NA,
                         simputeCov.error.B = NA,
                         simputeCov.rank = NA,
                         simpute.time = NA,
                         simpute.lambda.1 = NA,
                         simpute.error.test = NA,
                         simpute.error.all = NA,
                         simpute.error.beta = NA,
                         simpute.error.B = NA,
                         simpute.rank = NA,
                         simputeL2.alpha=NA, simputeL2.lambda.1=NA, simputeL2.lambda.2=NA, simputeL2.error.test=NA,
                         simputeL2.error.all=NA, simputeL2.error.B=NA, simputeL2.error.beta=NA, simputeL2.rank=NA, simputeL2.time=NA,
                         simputeKF.alpha=NA, simputeKF.lambda.1=NA, simputeKF.lambda.2=NA, simputeKF.error.test=NA,
                         simputeKF.error.all=NA, simputeKF.error.B=NA, simputeKF.error.beta=NA, simputeKF.rank=NA, simputeKF.time=NA)
   
   for(i in 1:length(dim)){
      
      if(missingness == 0){
         gen.dat <- generate_simulation_data_mao(n1=dim[i],n2=dim[i],m=5,r=10, seed=2023)
      }else
         gen.dat <- generate_simulation_data_ysf(2,dim[i],dim[i],10,10, missing_prob = missingness,coll=coll)
      results$true_rank[i] = gen.dat$rank
      set.seed(2023)
      # fit 1. Mao
      print(i)
      start_time = Sys.time()
      cv.out <- Mao.cv(gen.dat$A, gen.dat$X, gen.dat$Y, gen.dat$W,
                       n_folds=5, 
                       lambda.1_grid = lambda.1_grid,
                       lambda.2_grid = lambda.2_grid,
                       alpha_grid = alpha_grid,
                       numCores = ncores,n1n2_optimized = TRUE,theta_estimator = theta_default)
      mao.out <- Mao.fit(gen.dat$Y, gen.dat$X, gen.dat$W, cv.out$best_parameters$lambda.1, 
                         cv.out$best_parameters$lambda.2, cv.out$best_parameters$alpha, 
                         theta_estimator = theta_default)
      results$Mao.time[i] = round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
      # results out of fit 1
      results$Mao.alpha[i] = cv.out$best_parameters$alpha
      results$Mao.lambda.1[i] = cv.out$best_parameters$lambda.1
      results$Mao.lambda.2[i] = cv.out$best_parameters$lambda.2
      results$Mao.error.test[i] = test_error(mao.out$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
      results$Mao.error.all[i] = test_error(mao.out$A_hat, gen.dat$A)
      results$Mao.error.B[i] = test_error(mao.out$B_hat, gen.dat$B)
      results$Mao.error.beta[i] = test_error(mao.out$beta_hat, gen.dat$beta)
      results$Mao.rank[i] = mao.out$rank
      print(".")
      #----------------------------------------------------------
      # soft Impute model without covariates
      set.seed(2023)
      # validation set to be used for the next two models
      W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
      Y_train <- gen.dat$Y * W_valid
      Y_valid <- gen.dat$Y[W_valid==0]
      #---------------------------------------------------------
      start_time = Sys.time()
      sout <- simpute.orig(Y_train, W_valid, gen.dat$Y, trace=FALSE, rank.limit = 30,print.best=FALSE, rank.step = 4)
      results$simpute.time[i] =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
      results$simpute.lambda.2[i] = sout$lambda
      results$simpute.error.test[i] = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
      results$simpute.error.all[i] = test_error(sout$A_hat, gen.dat$A)
      results$simpute.rank[i] = sout$rank_A 
      print("..")
      #----------------------------------------------------------------------------
      # soft Impute model with covariates
      start_time = Sys.time()
      sout <- simpute.cov.cv(Y_train, gen.dat$X, W_valid, Y_valid, trace=FALSE, rank.limit = 30, 
                             print.best=FALSE, rank.step=4, type="als", lambda1=0, tol=2)
      results$simputeCov.time[i] =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
      results$simputeCov.lambda.1[i] = 0
      results$simputeCov.lambda.2[i] = sout$lambda
      results$simputeCov.error.test[i] = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
      results$simputeCov.error.all[i] = test_error(sout$A_hat, gen.dat$A)
      results$simputeCov.error.beta[i] = test_error(sout$beta_hat, gen.dat$beta)
      results$simputeCov.error.B[i] = test_error(sout$B_hat, gen.dat$B)
      results$simputeCov.rank[i] = sout$rank_A 
      print("...")
      #-------------------------------------------------------------------------------------
      # Soft Impute with Covariates and With L2 regularization on the covariates
      set.seed(2023)
      start_time = Sys.time()
      sout <- simpute.cov.cv(Y_train, gen.dat$X, W_valid, Y_valid, trace=FALSE, rank.limit = 30, 
                             print.best=FALSE, rank.step=4, type="als", lambda1=0, tol=2)
      sout <- simpute.cov.cv.lambda1(Y_train, gen.dat$X, W_valid, Y_valid, sout$lambda, sout$rank.max, print.best = FALSE,
                                    trace=FALSE, lambda1.grid = seq(0,20,length.out=20) ,n1n2 = 1, warm=NULL)
      
      results$simputeL2.time[i] =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
      results$simputeL2.lambda.1[i] = sout$lambda1
      results$simputeL2.lambda.2[i] = sout$lambda2
      results$simputeL2.error.test[i] = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
      results$simputeL2.error.all[i] = test_error(sout$A_hat, gen.dat$A)
      results$simputeL2.error.beta[i] = test_error(sout$beta_hat, gen.dat$beta)
      results$simputeL2.error.B[i] = test_error(sout$B_hat, gen.dat$B)
      results$simputeL2.rank[i] = sout$rank_A
      print("....")
      #-------------------------------------------------------------------------------------
      # Soft Impute with Covariates and With L2 regularization on the covariates and K-fold cross-validation
      set.seed(2023)
      start_time = Sys.time()
      
      sout <- simpute.cov.kfold(gen.dat$Y, gen.dat$X, gen.dat$W, n_folds = 3, print.best = FALSE,
                                 trace=FALSE, rank.limit = 30, lambda1=0,n1n2 = 1, warm=NULL,tol = 2)
      sout <- simpute.cov.kfold.lambda1(gen.dat$Y, gen.dat$X, gen.dat$W, sout$lambda2, n_folds = 3, print.best = FALSE, 
                                         trace=FALSE,lambda1.grid = seq(0,20,length.out=20) ,n1n2 = 1, warm=NULL,
                                         J=c(sout$J))
      
      results$simputeKF.time[i] =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
      results$simputeKF.lambda.1[i] = sout$lambda1
      results$simputeKF.lambda.2[i] = sout$lambda2
      results$simputeKF.error.test[i] = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
      results$simputeKF.error.all[i] = test_error(sout$A_hat, gen.dat$A)
      results$simputeKF.error.beta[i] = test_error(sout$beta_hat, gen.dat$beta)
      results$simputeKF.error.B[i] = test_error(sout$B_hat, gen.dat$B)
      results$simputeKF.rank[i] = sout$rank_A
      print(".....")
      #--------------------------------------------------------------------------------
      # saving plots to disk
      if(plot==TRUE){
         filename = paste0(graph_label,"_theta" ,missingness, c("_A","_beta", "_B"), "_dim",dim[i],"_coll",coll)
         plot_actual.vs.estimated.v2(gen.dat$A[gen.dat$W==0], mao.out$A_hat[gen.dat$W==0], fit_mao_x$Ahat[gen.dat$W==0],
                                     "New", "Old", "A", expression(hat("A")), "",sample.size = 10000, tofile=tofile, filename=filename[1])
         plot_actual.vs.estimated.v2(gen.dat$beta, mao.out$beta_hat, fit_mao_x$betahat,
                                     "New", "Old", "Beta", expression(hat("Beta")), "",sample.size = 10000, tofile=tofile, filename=filename[2])
         plot_actual.vs.estimated.v2(gen.dat$B, mao.out$B_hat, fit_mao_x$Bhat,
                                     "New", "Old", "B", expression(hat("B")), "",sample.size = 10000, tofile=tofile, filename=filename[3])
      }
      print(results[i,])
      
   }
   if(missingness == 0){
      filename = "compare_Mao_Youssef_Implementation_MaoSim.csv"
   }else
      filename = paste0("compare_Mao_Youssef_Implementation_YsfSim_",round(missingness*100),
                        "_coll_",coll, ".csv")
   
   write.csv(results, file = paste0(data_dir, filename), row.names = FALSE)
    
   #----------------------------
}

alpha_grid = c(1)
ncores = 1
 
compare_and_save(0.8, TRUE, plot = FALSE, tofile = TRUE,
                 lambda.1_grid = seq(0,2,length=20),lambda.2_grid = seq(.9, 0, length=20),
                 alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0.9, TRUE, plot = FALSE, tofile = TRUE,
                 lambda.1_grid = seq(0,2,length=20),lambda.2_grid = seq(.9, 0, length=20),
                 alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0.8, FALSE, plot = FALSE, tofile = TRUE,
                  lambda.1_grid = seq(0,2,length=20),lambda.2_grid = seq(.9, 0, length=20),
                  alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0.9, FALSE, plot = FALSE, tofile = TRUE,
                 lambda.1_grid = seq(0,2,length=20),lambda.2_grid = seq(.9, 0, length=20),
                 alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0,  FALSE, plot = FALSE, tofile = TRUE,
                 lambda.1_grid = seq(0,2,length=20), lambda.2_grid = seq(.9, 0, length=20), 
                 alpha_grid = alpha_grid, ncores=ncores)

#------------------------------------------
image_files <- function(miss,coll, dim, var) 
   paste0("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/Mao/saved_data/rds_plots/_theta",
          miss,"_",var,"_dim",dim,"_coll",coll,".rds")
 