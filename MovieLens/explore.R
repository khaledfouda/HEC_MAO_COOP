setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
suppressMessages(suppressWarnings(source("./code_files/import_lib.R", local = FALSE)))
#https://paperswithcode.com/sota/collaborative-filtering-on-movielens-100k

# ratings <-
#  read_csv("MovieLens/ml-latest/ratings.csv", show_col_types = FALSE)
# # links <- read_csv("MovieLens/ml-latest/links.csv",show_col_types = FALSE)
# movies <-
#  read_csv("MovieLens/ml-latest/movies.csv", show_col_types = FALSE)
# genome_scores <-
#  read_csv("MovieLens/ml-latest/genome-scores.csv", show_col_types = FALSE)
# genome_tags <-
#  read_csv("MovieLens/ml-latest/genome-tags.csv", show_col_types = FALSE)
# tags <-
#  read_csv("MovieLens/ml-latest/tags.csv", show_col_types = FALSE)
#
#
# head(tags)
# head(genome_tags)
# head(genome_scores)
# head(movies)
# head(ratings)

# read the data, and clean it. >>>>
dtrain <- read_tsv(
  "MovieLens/ml-100k/ua.base",
  col_names = c("user_id", "movie_id", "rating", "timestamp"),
  show_col_types = FALSE
)
dtest <- read_tsv(
  "MovieLens/ml-100k/ua.test",
  col_names = c("user_id", "movie_id", "rating", "timestamp"),
  show_col_types = FALSE
)
dX <- read_delim(
  "MovieLens/ml-100k/u.user",
  delim = "|",
  col_names = c("user_id", "age", "sex", "occupation", "zipcode"),
  show_col_types = FALSE
)
dX %>%
  arrange(user_id) %>%
  transmute(Male = sex == "M",
            age = age) %>%  #(age - mean(age)) / sd(age)) %>%
  # mutate(Male = (Male - mean(Male)) / sd(Male)) %>%
  as.matrix ->
  dX
head(dtrain)
dim(dtrain)
head(dtest)
dim(dtest)
head(dX)
dim(dX)
#-------------------------- drop movies with no ratings
bad_movies <- (dtest %>% 
  filter(!movie_id %in% dtrain$movie_id))$movie_id
bad_movies
dtest %<>% filter(! movie_id %in% bad_movies)
dtrain %<>% filter(! movie_id %in% bad_movies)
# reset the index 
dtrain %<>%
  #mutate(movie_id = as.integer(factor(movie_id)))
  mutate(movie_id = ifelse(movie_id < bad_movies[1], movie_id, movie_id-1)) %>% 
  mutate(movie_id = ifelse(movie_id < bad_movies[2], movie_id, movie_id-1)) 
dtest %<>%
  #mutate(movie_id = as.integer(factor(movie_id)))
  mutate(movie_id = ifelse(movie_id < bad_movies[1], movie_id, movie_id-1)) %>% 
  mutate(movie_id = ifelse(movie_id < bad_movies[2], movie_id, movie_id-1)) 


#----------------------------------------
# prepare sparse matrix and X decomposition
dtrain.mat <- sparseMatrix(
  i = dtrain$user_id,
  j = dtrain$movie_id,
  x = dtrain$rating,
  dims = c(max(dtrain$user_id), max(dtrain$movie_id))
)

dtrain.mat <- as(dtrain.mat, "Incomplete")
dim(dtrain.mat)
X_r <- reduced_hat_decomp(dX, 1e-2)
#-----------------------------------------------------------


##----------------------------------------------------------------
#
# #
# fits <- simpute.als.fit_splr(y=dtrain.mat, svdH=X_r$svdH,  trace=T, J=10,
#                              thresh=1e-6, lambda=10, init = "naive",
#                              final.svd = T,maxit = 500, warm.start = NULL)
# length(fits$xbeta.obs)
# length(dtrain.mat@x)
#
# dtrain.mat
#
# naive.y = naive_MC(as.matrix(dtrain.mat))
# naive.y[1:5,1:5]
# propack.svd(naive.y, 1)
# svd(naive.y)$d[1]
#
# min(as.matrix(dtrain.mat))
# table(dtrain.mat@x,useNA = "ifany")
#
# gen.dat <-  generate_simulation_data_mao(400,300,10,10,2024,method="MAR")
#
#=-======================================================
## prepare validation data
biScale.out <-biScaleMatrix(as.matrix(dtrain.mat)) 
dtrain.Y = biScale.out$scaledMat
dtrain.mat <- as(dtrain.Y, "Incomplete")

dW <- !is.na(dtrain.Y)
sum(dW == 0) / length(dW)
dtrain.Y[is.na(dtrain.Y)] = 0
W_valid <- matrix.split.train.test(dW, testp = 0.2)
Y_train = (dtrain.Y * W_valid)
Y_valid = dtrain.Y[W_valid == 0]
test_error <- RMSE_error
#------------------------------------------------------------------
# fit splr model and get RMSE
start_time <- Sys.time()
best_fit = simpute.cov.cv_splr(
  Y_train,
  X_r,
  Y_valid,
  W_valid,
  dtrain.Y,
  trace = T,
  thresh = 1e-6,
  n.lambda = 30,
  rank.step = 4,
  rank.limit = 50
)
print(paste("Execution time is", round(as.numeric(
  difftime(Sys.time(), start_time, units = "secs")
), 2), "seconds"))
fit1 = best_fit$fit1
fit2 = best_fit$fit2
# get estimates and validate
beta =  fit2$u %*% (fit2$d * t(fit2$v))
M = fit1$u %*% (fit1$d * t(fit1$v))
A = M + X_r$X %*% t(beta)
A = revertBiScaledMatrix(A, biScale.out)

#qr(A)$rank
dtest$preds <- A[cbind(dtest$user_id, dtest$movie_id)]
RMSE_error(dtest$preds, dtest$rating)

mean(beta[, 1])
mean(beta[, 2])
sd(beta[, 1])
sd(beta[, 2])
#------------------------------------------------------------------
# fit the original soft Impute and get the results
start_time = Sys.time()

#Y_train[Y_train==0] <- NA
#xc=biScale(Y_train,col.scale=FALSE,row.scale=FALSE,trace=TRUE)
#biScale(dtrain.mat) -> scaleout
#scaleout[1:10,1:10]
#Y_train[1:10,1:10]


sout <- simpute.orig(
  Y_train,
  W_valid,
  dtrain.Y,
  trace = TRUE,
  rank.limit = 50,
  print.best = FALSE,
  rank.step = 7,
  biscale=F
)
print(paste("Execution time is", round(as.numeric(
  difftime(Sys.time(), start_time, units = "secs")
), 2), "seconds"))

sout$A_hat = revertBiScaledMatrix(sout$A_hat, biScale.out)
dtest$preds_orig <- sout$A_hat[cbind(dtest$user_id, dtest$movie_id)]

RMSE_error(dtest$preds_orig, dtest$rating)
sout$rank_A
#----------------------------------------------------------
# repeat with mao
lambda.1_grid = seq(0, 2, length = 20)
lambda.2_grid = seq(0.9, 0, length = 20)
alpha_grid = seq(0.992, 1, length = 20)
start_time = Sys.time()
cv.out <- suppressMessages(suppressWarnings(
  Mao.cv(
    dtrain.Y,
    X_r$X,
    dtrain.Y,
    dW,
    n_folds = 3,
    lambda.1_grid = lambda.1_grid,
    lambda.2_grid = lambda.2_grid,
    alpha_grid = alpha_grid,
    numCores = 1,
    n1n2_optimized = T,
    theta_estimator = MaoUniWeights
  )
))
mao.out <- suppressMessages(suppressWarnings(
  Mao.fit(
    dtrain.Y,
    X_r$X,
    dW,
    cv.out$best_parameters$lambda.1,
    cv.out$best_parameters$lambda.2,
    cv.out$best_parameters$alpha,
    theta_estimator = MaoUniWeights,
    n1n2_optimized = T
  )
))
print(paste("Execution time is", round(as.numeric(
  difftime(Sys.time(), start_time, units = "secs")
), 2), "seconds"))
mao.out$A_hat = revertBiScaledMatrix(mao.out$A_hat, biScale.out)
dtest$preds_mao <-
  mao.out$A_hat[cbind(dtest$user_id, dtest$movie_id)]
RMSE_error(dtest$preds_mao, dtest$rating)
mao.out$rank
cv.out$best_parameters
#---------------------------------------------------------------
# applying the K-fold method
#
#
# start_time = Sys.time()
# soutl2 <-
# simpute.cov.kfold(
#  dtrain.Y,
#  X_r$X,
#  dW,
#  n_folds = 3,
#  print.best =TRUE,
#  trace = TRUE,
#  rank.limit = 30,
#  lambda1 = 0,
#  n1n2 = 1,
#  warm = NULL,
#  tol = 2,
#  rank.step = 2,
#  n.lambda = 30
# )
# soutl2. <-
#  simpute.cov.kfold.lambda1(
#   dtrain.Y,
#   X_r$X,
#   dW,
#   soutl2$lambda2,
#   n_folds = 3,
#   print.best = TRUE,
#   trace = TRUE,
#   lambda1.grid = lambda.1_grid ,
#   n1n2 = 1,
#   warm = NULL,
#   J = c(soutl2$J)
#  )
# print(paste("Execution time is", round(as.numeric(
#   difftime(Sys.time(), start_time, units = "secs")
# ), 2), "seconds"))
# dtest$preds_l2 <-
#  soutl2.$A_hat[cbind(dtest$user_id, dtest$movie_id)]
# RMSE_error(dtest$preds_l2, dtest$rating)
# soutl2.$rank
#