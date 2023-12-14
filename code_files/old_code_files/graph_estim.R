plot_actual.vs.estimated <- function(Actual, Predicted, xlab, ylab,title1="Comparison of Actual and Predicted Values"){
   
   data <- data.frame(
      Actual = as.vector(Actual),
      Predicted = as.vector(Predicted)
   )
   lm_fit <- lm(Predicted ~ Actual, data = data)
   coefficients <- coef(lm_fit)
   
   # Create the equation text
   equation_text <- sprintf("y = %.2fx + %.3f", coefficients[2], coefficients[1])
   
   ggplot(data, aes(x = Actual, y = Predicted)) + 
      geom_hex(color = "grey80",bins=200) +  
      geom_smooth(method = "lm", aes(color = "LM Fit"), se = TRUE) +
      #geom_line(data = data.frame(x = range(data$Actual), y = range(data$Predicted)),
      #          aes(x = x, y = y, linetype = "Perfect Fit"), alpha=0, color = "white") +  
      annotate("text", x = Inf, y = Inf, label = equation_text, hjust = 2.1, size=5, vjust = 3.5, color = "red", size = 2.5) + 
      scale_color_manual("", values = c("LM Fit" = "red")) + 
      geom_abline(intercept = 0, slope = 1, linetype = 1, linewidth=1, color = "black") + 
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
      guides(fill = FALSE) 
}
