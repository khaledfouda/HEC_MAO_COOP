setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/")

source("./code_files/Mao_import_lib.R")
data_dir = "./Mao/saved_data/"


results <- read.csv(paste0(data_dir,"sim_results_MAO_table1_B30_1.csv"))
#results
kable(results,format = "pipe")

 
# a

missingness = seq(.3,.9,.1)
names = c("theta_default", "theta_simple", "theta_random")



results <- rbind( read.csv(paste0(data_dir,"a_",names[1],"_.csv")) %>% mutate(label="Logistic model"),
                  read.csv(paste0(data_dir,"a_",names[2],"_.csv")) %>% mutate(label="Proportion within column"),
                  read.csv(paste0(data_dir,"a_",names[3],"_.csv")) %>% mutate(label="Proportion overall")) 


results %>%
   ggplot(aes(x = missingness, y=test_errors, group=label, color=label)) +
   geom_line() + 
   geom_point(size = 3, shape=17) +
   theme_bw() +  
   scale_color_manual(values = c("red","pink", "cyan")) +
   labs(x = "missing %", y = "RMSE", title = "Model 2; n1=n2=300; m=5;", 
        subtitle = "RMSE on the missing data.", color="Theta Method") +
   theme(legend.position = "bottom") 

# a2

missingness = seq(.3,.9,.1)
names = c("theta_default", "theta_simple", "theta_random")



results <- rbind( read.csv(paste0(data_dir,"a_",names[1],"_.csv")) %>% mutate(label="Logistic model"),
                  read.csv(paste0(data_dir,"a_",names[2],"_.csv")) %>% mutate(label="Proportion within column"),
                  read.csv(paste0(data_dir,"a_",names[3],"_.csv")) %>% mutate(label="Proportion overall")) 


results %>%
   ggplot(aes(x = missingness, y=test_errors, group=label, color=label)) +
   geom_line() + 
   geom_point(size = 3, shape=17) +
   theme_bw() +  
   scale_color_manual(values = c("red","pink", "cyan")) +
   labs(x = "missing %", y = "RMSE", title = "Model 2; n1=n2=300; m=5;", 
        subtitle = "RMSE on the missing data.", color="Theta Method") +
   theme(legend.position = "bottom") 

# b



graphs_A = list()
graphs_beta = list()
graphs_B = list()
part_b_data <- read.csv(paste0(data_dir,"sim_results_MAO_table1_B30_1.csv"))[c(1,3,4),]



for(i in 1:nrow(part_b_data)){
   gen.dat <- generate_simulation_data_mao(part_b_data$dimension[i],part_b_data$dimension[i])
   #gen.dat <- generate_simulation_data_ysf(2, results$dim[i],results$dim[i],5,5,0.5)
   mao.out <- Mao.fit(gen.dat$Y, gen.dat$X, gen.dat$W, part_b_data$lambda_1[i], 
                      part_b_data$lambda_2[i], part_b_data$alpha[i],n1n2_optimized = T,theta_estimator = theta_default)
   graphs_beta[[i]] = plot_actual.vs.estimated(gen.dat$beta, mao.out$beta_hat, "Beta", expression(hat("Beta")),title="")
   graphs_B[[i]] = plot_actual.vs.estimated(gen.dat$B[gen.dat$W==0], mao.out$B_hat[gen.dat$W==0], "B", expression(hat("B")),title="")
   graphs_A[[i]] = plot_actual.vs.estimated(gen.dat$A[gen.dat$W==0], mao.out$A_hat[gen.dat$W==0], "A", expression(hat("A")),title="")
}

combined_plot <- (graphs_A[[1]] | graphs_beta[[1]] | graphs_B[[1]]) / 
   (graphs_A[[2]] | graphs_beta[[2]] | graphs_B[[2]]) / 
   (graphs_A[[3]] | graphs_beta[[3]] | graphs_B[[3]]) +
   plot_layout(guides = "collect") +
   plot_annotation(title = "Comparing original and estimated values (Only data with W=0 is included)",
                   subtitle = "with n1 = n2 = 300, 500, and 700; m=5, and probability of missingness is 0.8",)

combined_plot & theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),legend.position = "none") 


# c


names = c("theta_default", "theta_simple", "theta_random")
results <- rbind( read.csv(paste0(data_dir,"a_",names[1],"_.csv")) %>% mutate(label="Logistic model"),
                  read.csv(paste0(data_dir,"a_",names[2],"_.csv")) %>% mutate(label="Proportion within column"),
                  read.csv(paste0(data_dir,"a_",names[3],"_.csv")) %>% mutate(label="Proportion overall")) 


results %>%
   ggplot(aes(x = missingness)) +
   geom_line(aes( y = alpha, color="Alpha")) + 
   geom_point(aes( y = alpha), size = 3) +
   geom_line(aes( y = lambda.1, color="Lambda 1")) + 
   geom_point(aes( y = lambda.1), size = 3) +
   geom_line(aes( y = lambda.2, color="Lambda 2")) + 
   geom_point(aes( y = lambda.2), size = 3) +
   theme_bw() +  
   scale_color_manual(values = c("Alpha" = "darkgreen", "Lambda 1" = "darkred",
                                 "Lambda 2"="darkblue")) +
   scale_shape_manual(values = c("Alpha" = 17, "Lambda 1" = 15, "Lambda 2" = 18)) +
   labs(x = "missing %", y = "Parameter Value", title = "Hyperparameters vs missing percentage", 
        subtitle = "Settings are the same as the previous plot", color="Hyperparameter") +
   theme(legend.position = "bottom") +
   facet_wrap(~label)

# d

RMSE_vals = readRDS(paste0(data_dir,"part_d_data.rds"))
lambda.1 = seq(0,2,length=10)
lambda.2 = seq(.9, 0.1, length=10)
alpha = seq(0.992, 1, length=10)

data.frame(RMSE=RMSE_vals, labels=rep(c("Lambda 1", "Lambda 2", "Alpha"),each=10),
           values=c(lambda.1, lambda.2, alpha)) %>%
   ggplot(aes(x = values, y = RMSE, group = labels, color = labels, shape = labels)) +
   geom_line() +  # Add lines
   geom_point(size = 3) +  # Add points with a specified size
   scale_color_manual(values = c("Alpha" = "red", "Lambda 1" = "green", "Lambda 2" = "blue")) +
   scale_shape_manual(values = c("Alpha" = 17, "Lambda 1" = 15, "Lambda 2" = 18)) +
   theme_minimal() +  # Use a minimal theme
   labs(x = "Hyperparameter value", y = "RMSE", title = "Hyperparameters change vs RMSE", 
        subtitle = "Model 2; n1=n2=700; m=5; RMSE on the missing data.") +
   theme(legend.position = "bottom")  
