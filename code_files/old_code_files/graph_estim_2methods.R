plot_actual.vs.estimated.v2 <- function(Actual, Predicted1, Predicted2, label1, label2,
                                        xlab, ylab, title1="Comparison of Actual and Predicted Values",
                                        tofile=FALSE, img_folder="default", sample.size=10000, filename=""){
   
   if(img_folder=="default" & tofile==TRUE)
      img_folder =  "/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/Mao/saved_data/rds_plots/"
   
   if(tofile==TRUE)
      stopifnot(filename != "")
   
   data <- data.frame(
      Actual = as.vector(Actual),
      Predicted1 = as.vector(Predicted1),
      Predicted2 = as.vector(Predicted2)
   )
   lm_fit <- lm(Predicted1 ~ Actual, data = data)
   coefficients1 <- coef(lm_fit)
   lm_fit <- lm(Predicted2 ~ Actual, data = data)
   coefficients2 <- coef(lm_fit)
   
   # middle point to align text
   axis_range = range(data)
   mid_pts = c(median(data$Actual), quantile(c(data$Predicted1,data$Predicted1),c(0.1,0.9)))
   
   # Create the equation text
   equation_text1 <- sprintf("y = %.2fx + %.3f", coefficients1[2], coefficients1[1])
   equation_text2 <- sprintf("y = %.2fx + %.3f", coefficients2[2], coefficients2[1])
   
   data %>% 
      sample_n(min(sample.size,nrow(data))) %>% 
   ggplot(aes(x = Actual, y=Actual)) + 
      
      geom_abline(intercept = 0, slope = 1, linetype = 1, linewidth=1, color = "black") + 
      geom_point(aes(y=Predicted1), color = "pink",alpha=0.3) +
      geom_point(aes(y=Predicted2), color = "#b5eff5", alpha=0.3) +
      #geom_hex(color = "grey80",bins=200) +  
      
      geom_smooth(aes(y=Predicted1, color = "LM Fit1"),data=data, method = "lm", formula= y~x, se = TRUE) +
      geom_smooth( aes(y=Predicted2, color = "LM Fit2"),data=data, method = "lm",formula= y~x, se = TRUE) +
      
      #geom_line(data = data.frame(x = range(data$Actual), y = range(data$Predicted)),
      #          aes(x = x, y = y, linetype = "Perfect Fit"), alpha=0, color = "white") +  
      
      annotate("text",  x = mid_pts[1], y =mid_pts[2], label = equation_text1,  size=5, color = "#de1425") + 
      annotate("text", x = mid_pts[1], y =mid_pts[3], label = equation_text2, size=5, color = "#217880") +
      
      scale_color_manual("", values = c("LM Fit2" = "#217880","LM Fit1" = "#de1425"),labels=c(label1,label2)) +
      
      #xlim(axis_range)+
      #ylim(axis_range)+
      #scale_linetype_manual("", values = c("Perfect Fit" = 1)) +
      theme_minimal() + 
      theme(
         legend.position = "top", 
         legend.title = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_blank(),
         panel.border = element_blank()) + 
      labs(title = title1,
           x = xlab,
           y = ylab) +
      guides(fill = FALSE) -> generated_graph
   
   if(tofile==TRUE){
      filepath = paste0(img_folder, str_split(filename,".rds")[[1]][1], ".rds")
      print(paste0("Saving generated graph to ", filepath))
      saveRDS(generated_graph, filepath)
      #return()
   }else
      return(generated_graph)
}

# example use:
#X = matrix(runif(100000,0,3),nrow=50000)
#P1 = X + matrix(rnorm(100000,0,3),nrow=50000)
#P2 = X + matrix(rnorm(100000,1,2),nrow=50000)
#plot_actual.vs.estimated.v2(X, P1, P2, "Method 1", "Method 2", "Actual", "Predicted",sample.size = 10000,tofile=TRUE, filename="test.rds")
#readRDS("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/Mao/saved_data/rds_plots/test.rds")
#---------------------------------------------------------------

plot_A_B_beta.v2 <- function(name_function, missingness, coll, tofile=FALSE)
{
   title = sprintf("Comparing original and estimated values using the two methods and with a probability of missing values of %.2f",
                   missingness)
   subtitle = paste0("with n1 = n2 = 400, 600, 800, and 1000 and with ",ifelse(coll==TRUE,"","no "),"collinearity.")
   dim = seq(400,1000,200)
   
   graphs_A = list()
   graphs_beta = list()
   graphs_B = list()
   
   
   for(i in 1:length(dim)){
      graphs_A[[i]] = readRDS(name_function(missingness, coll, dim[i], "A"))
      graphs_beta[[i]] = readRDS(name_function(missingness, coll, dim[i], "beta")) 
      graphs_B[[i]] = readRDS(name_function(missingness, coll, dim[i], "B"))
   }
   
   
   combined_plot <- (graphs_A[[1]] | graphs_beta[[1]] | graphs_B[[1]]) / 
      (graphs_A[[2]] | graphs_beta[[2]] | graphs_B[[2]]) / 
      (graphs_A[[3]] | graphs_beta[[3]] | graphs_B[[3]]) /
      (graphs_A[[4]] | graphs_beta[[4]] | graphs_B[[4]]) +
      plot_layout(guides = "collect") +
      plot_annotation(title = title,
                      subtitle = subtitle)
   
   combined_plot <- combined_plot &
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = "top") #,legend.position = "none"
   
   return(combined_plot)
   
}
image_files <- function(miss,coll, dim, var) 
   paste0("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/Mao/saved_data/rds_plots/_theta",
          miss,"_",var,"_dim",dim,"_coll",coll,".rds")
#plot_A_B_beta.v2(image_files, 0, FALSE)
#plot_A_B_beta.v2(image_files, 0.8, FALSE)
#plot_A_B_beta.v2(image_files, 0.9, FALSE)
#plot_A_B_beta.v2(image_files, 0.8, TRUE)
#plot_A_B_beta.v2(image_files, 0.9, TRUE)


