
source("./code_files/Mao_import_lib.R")
library(softImpute)
dim = c(600)
missingness = 0.9
i=1
coll=TRUE

if(missingness == 0){
   gen.dat <- generate_simulation_data_mao(n1=dim[i],n2=dim[i],m=5,r=10, seed=2023)
}else
   gen.dat <- generate_simulation_data_ysf(2,dim[i],dim[i],10,10, missing_prob = missingness,coll=coll)

W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_train = (gen.dat$Y * W_valid) 
Y_train[Y_train==0] = NA

start_time <- Sys.time()
sout1 <- simpute.als.cov(Y_train, gen.dat$X, 3, 1e-3, 30,trace.it = FALSE)
sout1 <- simpute.als.cov(Y_train, gen.dat$X, 3, 1e-3, 30,trace.it = FALSE)
sout1 <- simpute.als.cov(Y_train, gen.dat$X, 3, 1e-3, 30,trace.it = FALSE)
round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))


sout2 <- simpute.svd.cov(Y_train, gen.dat$X, 3, 1e-3, 30,trace.it = TRUE)



lambda1.grid <- seq(0,100,length.out=4)
model_output <- list()
validation_errors <- rep(NA, length(lambda1.grid))

for(i in 1:length(lambda1.grid)){
   model_output[[i]] <- simpute.cov.cv(gen.dat$Y*W_valid, gen.dat$X, W_valid, gen.dat$Y,
                                       trace=FALSE, rank.limit = 30, lambda1=lambda1.grid[i])
}
validation_error <- lapply(model_output, function(d) d$error) %>% unlist()
best_fit = model_output[[which.min(validation_error)]]

#----------

lambda1.grid <- seq(0, 20, length.out=10)

softImputeALS_L2 <- function(Y.train, Y.valid, W.valid, X, lambda1.grid=seq(0,20,length.out=10), n1n2=1, no_cores=NA, max_cores=20){
   
   if(is.na(no_cores)) no_cores = length(lambda1.grid)
   no_cores = min(max_cores, no_cores)
   
   model_output <- mclapply(lambda1.grid, function(lambda) {
      simpute.cov.cv(Y.train, X, W.valid, Y.valid,
                     trace=FALSE, rank.limit = 30, lambda1=lambda,n1n2 = n1n2)
   }, mc.cores = no_cores)
   
   valid_errors <- unlist(lapply(model_output, function(d) d$error))
   best_index = which.min(valid_errors)
   best_fit <- model_output[[best_index]]
   list(best_score=valid_errors[best_index], best_fit=best_fit, lambda1 = lambda1.grid[best_index])
}

results <- softImputeALS_L2(gen.dat$Y*W_valid, gen.dat$Y[W_valid==0], W_valid, gen.dat$X)

no_cores = min(20, length(lambda1.grid))

model_output <- mclapply(lambda1.grid, function(lambda) {
   simpute.cov.cv(gen.dat$Y*W_valid, gen.dat$X, W_valid, gen.dat$Y,
                  trace=FALSE, rank.limit = 30, lambda1=lambda,n1n2 = 1)
}, mc.cores = no_cores)
validation_error <- unlist(lapply(model_output, function(d) d$error))
best_fit <- model_output[[which.min(validation_error)]]
validation_error
test_errors <- unlist(lapply(model_output, function(d) test_error(d$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])))

# [1] 0.3043862 [0.2720237] 0.2746299 0.2816688 0.2904483 0.2963708 0.3022300 0.3068791 0.3126155 0.3168772 # 1 / 0-100
#     0.3021313 0.2896868 0.2810955 0.2764015 0.2736201 0.2730740 [0.2716675] 0.2723489 0.2729733 0.2734022     # 1 / 0-20
# [1] [0.3028395] 0.3511031 0.3544883 0.3553735 0.3507690 0.3494529 0.3524645 0.3510712 0.3543921 0.3490639 # 2 / 0-100
# [1] [0.3022256] 0.3284005 0.3440804 0.3456595 0.3452691 0.3488261 0.3493347 0.3461140 0.3515047 0.3471571 # 2 / 0-2
# [1] [0.3026262] 0.3499032 0.3527240 0.3493431 0.3512216 0.3526009 0.3505592 0.3525565 0.3537018 0.3490670 # 3 / 0-100
# [1] [0.3010404] 0.3501780 0.3525480 0.3524318 0.3509615 0.3490343 0.3492069 0.3505750 0.3488609 0.3485438   # 3 / 0-1
#      0.3028448 0.3522959 0.3509149 0.3502864 0.3504158 0.3486029 0.3508370 0.3544691 0.3503435 0.3500651
#      0.3024280 0.2873747 0.2782052 0.2737401 [0.2725309] 0.2728103 0.2733185 0.2741148 0.2751018 0.2767083 # 4 / 0-1
# [1]  0.3031338 0.2938073 0.2871283 0.2816428 0.2771472 0.2762637 0.2737229 [0.2721922] 0.2729892 0.2730445 # 4 / 0-0.5
X.X = t(gen.dat$X) %*% gen.dat$X
svd(gen.dat$X)$d[1]


test_error(sout1$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
test_error(sout1$beta_hat, gen.dat$beta)
test_error(sout1$B_hat, gen.dat$B)


gen.dat$Y*W_valid - sout2$B_hat

start_time <- Sys.time()
sout2 <- simpute.cov.cv.v2(gen.dat$Y*W_valid, gen.dat$X, W_valid, gen.dat$Y, trace=TRUE,
                           lambda1.optimize = TRUE,
                           lambda1.grid = seq(0,1000, length.out=20),
                           rank.limit = 30,type="als")
round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))

test_error(sout2$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
test_error(sout2$beta_hat, gen.dat$beta)
test_error(sout2$B_hat, gen.dat$B)
sout2$lambda1


sout3 <- simpute.cov.cv(gen.dat$Y*W_valid, gen.dat$X, W_valid, gen.dat$Y, trace=TRUE,
                           rank.limit = 30,type="als")
test_error(sout3$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
test_error(sout3$beta_hat, gen.dat$beta)
test_error(sout3$B_hat, gen.dat$B)
sout2$lambda1


sout <- simpute.orig(gen.dat$Y*W_valid, W_valid, gen.dat$Y, trace=TRUE, rank.limit = 30)
test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
test_error(sout$beta_hat, gen.dat$beta)
test_error(sout$B_hat, gen.dat$B)



lambda.1_grid = seq(0,3,length=20)
lambda.2_grid = seq(.9, 0, length=20)
alpha_grid = c(1)#seq(0.992, 1, length=10)

cv.out <- Mao.cv(gen.dat$A, gen.dat$X, gen.dat$Y, gen.dat$W,
                 n_folds=5, 
                 lambda.1_grid = lambda.1_grid,
                 lambda.2_grid = lambda.2_grid,
                 alpha_grid = alpha_grid,
                 numCores = 1,n1n2_optimized = TRUE,theta_estimator = theta_default)
mao.out <- Mao.fit(gen.dat$Y, gen.dat$X, gen.dat$W, cv.out$best_parameters$lambda.1, 
                   cv.out$best_parameters$lambda.2, cv.out$best_parameters$alpha, 
                   theta_estimator = theta_default)

test_error(mao.out$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
test_error(mao.out$beta_hat, gen.dat$beta)
test_error(mao.out$B_hat, gen.dat$B)
#-----------------------------------------------------------
# 1. Soft Impute with/without X
with_X = FALSE

if (with_X==TRUE){new_Y = gen.dat$Y - gen.dat$X %*% mao.out$beta_hat}else new_Y = gen.dat$Y
new_Y[gen.dat$W==0] = NA

#xs <- as(new_Y, "Incomplete")
lam0 <- lambda0(new_Y)
lamseq <- exp(seq(from=log(lam0), to=log(1), length=20))

fits <- as.list(lamseq)
ranks <- as.integer(lamseq)
rank.max <- 10
warm <- NULL


for(i in seq(along=lamseq)) {
   fiti <- softImpute(new_Y, lambda=lamseq[i], rank.max=rank.max, warm=warm)
   ranks[i] <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
   rank.max <- min(ranks[i]+2, 50)
   if (with_X==TRUE){ soft_estim  = complete(new_Y, fiti) +  gen.dat$X %*% mao.out$beta_hat}else soft_estim  =  complete(new_Y, fiti)
   
   err = test_error(soft_estim[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   warm <- fiti # warm start for next 
   fits[[i]] <- fiti
   cat(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
               i, lamseq[i], rank.max, ranks[i], err))
}
#---------------------------------------------
fitslo <- simpute.svd.cov(new_Y, X, W,lambda = 7, trace.it=FALSE, warm.start = fitslo)

v=as.matrix(fitslo$v)
vd=v*outer(rep(1,nrow(v)),fitslo$d)
soft_estim = fitslo$u %*% t(vd)  + X %*% fitslo$beta.estim
test_error(soft_estim[W==0], A[W==0])

fitslo$d
#-----------------------------------------------------------------------------------
# with covariates - loop


#-------------------------------



fits <- softImpute(new_Y, trace=FALSE)
full.Y <- complete(new_Y, fits)
soft_estim <- full.Y + X %*% mao.out$beta_hat
test_error(soft_estim[W==0], A[W==0])


dim(new_Y)
X <- X

qr(full.Y)$rank
fits$d
Y <- Y

#######################
#-----------------------------------------------------------------------------------
# with covariates - loop
new_Y = gen.dat$Y
new_Y[gen.dat$W==0] = NA

#xs <- as(new_Y, "Incomplete")
lam0 <- lambda0(new_Y)
lam0 <- lambda0.cov(new_Y, gen.dat$X, gen.dat$W)
#lam0 <- 40
lamseq <- seq(from=40, to=0, length=20)

fits <- as.list(lamseq)
ranks <- as.integer(lamseq)
rank.max <- 2
warm <- NULL
best_fit <- list(error=Inf, rank=NA, lambda=NA, rank.max=NA)

for(i in seq(along=lamseq)) {
   fiti <- simpute.svd.cov(new_Y, gen.dat$X, gen.dat$W,lambda = lamseq[i], J=rank.max, warm.start = warm)
   ranks[i] <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
   rank.max <- min(ranks[i]+2, 50)
   # get estimates
   v=as.matrix(fiti$v)
   vd=v*outer(rep(1,nrow(v)),fiti$d)
   soft_estim = fiti$u %*% t(vd)  + gen.dat$X %*% fiti$beta.estim
   #----------------------------
   err = test_error(soft_estim[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   warm <- fiti # warm start for next 
   fits[[i]] <- fiti
   cat(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
               i, lamseq[i], rank.max, ranks[i], err))
   #-------------------------
   # register best fir
   if(err < best_fit$error){
      
      best_fit$error = err
      best_fit$rank_B = ranks[i]
      best_fit$rank_A = qr(soft_estim)$rank
      best_fit$lambda = lamseq[i]
      best_fit$rank.max = rank.max
   } 
}
print(best_fit)

#-------------------------------
###################################################



optimize_lambda1 <- function(W, Y.minus.B, X, lambda1.grid, Y.valid, trace=FALSE, acc=4){
   #' W: Validation mask where Wij=0 if the cell belongs to the validation set and 1 if belongs to the training or test
   #' fiti: training object containing the SVD decomposition of the low-rank matrix (UDV) 
   #'  lambda1.grid: grid of possible values for lambda1
   #'  Y.valid: A vector of the true values of Y for the validation set.
   #'  ------------------
   #'  Returns: best.lambda1: optimal value
   #'           Xbeta: the term X^T beta corresponding to the optimal lambda
   X.X = t(X) %*% X
   XY <- t(X) %*% Y.minus.B
   validation_indices <- which(W == 0, arr.ind = TRUE)
   best_score = Inf
   best_lambda1 = NULL
   best_beta = NULL
   nr <- nrow(X.X)
   n1n2 = svd(X.X)$d[1]
   Y.train <- Y.minus.B
   
   for(lambda1 in lambda1.grid){
      beta_partial = solve(X.X + n1n2*lambda1* diag(1, nr)) %*% t(X)
      for(i in 1:10){
         #print(i)
         beta.estim = beta_partial %*% Y.train
         soft_estim = Y.minus.B  + X %*% beta.estim
         Y.train[W==0] = soft_estim[W==0]
         err = round(test_error(soft_estim[W==0], Y.valid),acc)
         #print(err)
         #print(test_error(beta.estim, gen.dat$beta))
      }
      #soft_estim = compute_test_estimates(validation_indices, Y.minus.B, X, beta.estim)
      if(err < best_score){
         best_score = err
         best_lambda1 = lambda1
         beta_partial = beta_partial
         best_beta = beta.estim
      }
         print(test_error(beta.estim, gen.dat$beta))
         print(paste(err, lambda1))
   }
   results = list(lambda1 = best_lambda1, score=best_score, beta_partial =beta_partial%*% t(X), beta=best_beta)
   results
}

W <- W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y.minus.B <- Y_train <- gen.dat$Y * W_valid
Y.valid <- gen.dat$Y[W_valid ==0]
lambda1.grid <- seq(0,2,length.out=10)
X <- gen.dat$X
acc <- 4
lambda1 <- 0


opt.out <- optimize_lambda1(W_valid, Y_train, gen.dat$X, seq(0,2,length.out=10), Y_valid)

test_error(sout1$beta_hat, gen.dat$beta)
