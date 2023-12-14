library(knitr)
library(kableExtra)
library(tidyverse)
library(magrittr)
require(foreach)
require(doParallel)
library(ggplot2)
library(hexbin)
library(patchwork)
library(softImpute)


#setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/")
path_to_code = "/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/code_files/"
path_to_data = "/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/Mao/saved_data/"

source(paste0(path_to_code,"Mao_fit.R"))
source(paste0(path_to_code,"Mao_cv.R"))
source(paste0(path_to_code,"Mao_sim.R")) 
source(paste0(path_to_code,"Ysf_sim.R")) 
source(paste0(path_to_code,"graph_estim.R"))
source(paste0(path_to_code,"graph_estim_2methods.R"))
source(paste0(path_to_code,"matrix_split_train_valid.R"))
source(paste0(path_to_code,"SoftImpute_fit_covariates.R"))
source(paste0(path_to_code,"SoftImputeALS_fit_covariates.R"))
source(paste0(path_to_code,"SoftImpute_cv_covariates.R"))