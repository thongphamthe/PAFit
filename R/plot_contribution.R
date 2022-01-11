
plot_contribution <- function(simulated_object,
                              original_result,
                              which_plot = "PA",
                              y_label = ifelse("PA" == which_plot,"Contribution of the rich-get-richer effect",
                                              "Contribution of the fit-get-richer effect"),
                              legend_pos_x = 0.75,
                              legend_pos_y = 0.9) {
  x <- y <- ymin <- ymax <- y_obs <- NULL
  if (class(simulated_object) != "Simulated_Data_From_Fitted_Model") {
     stop("simulated_object needs to be an object from generate_simulated_data_from_estimated_model")
  }
  if (class(original_result) != "Full_PAFit_result") {
    stop("original_result needs to be an object from joint_estimate")
  }
  result_list <- simulated_object$result_list
  M <- length(result_list)
  T <- length(result_list[[1]]$contribution$PA_contribution) + 1
PA_contrib_replicate  <- matrix(0,nrow = T - 1, ncol = M)
fit_contrib_replicate <- matrix(0,nrow = T - 1, ncol = M)

for (mm in 1:M) {
  PA_contrib_replicate[,mm]  <- result_list[[mm]]$contribution$PA_contribution
  fit_contrib_replicate[,mm] <- result_list[[mm]]$contribution$fit_contribution
}

simulated_mean_PA_contrib  <- rowMeans(PA_contrib_replicate)
simulated_sd_PA_contrib    <- apply(PA_contrib_replicate,MARGIN = 1, sd)
simulated_mean_fit_contrib <- rowMeans(fit_contrib_replicate)
simulated_sd_fit_contrib   <- apply(fit_contrib_replicate,MARGIN = 1, sd)
sd_log_PA                  <- simulated_sd_PA_contrib/simulated_mean_PA_contrib
sd_log_fit                 <- simulated_sd_fit_contrib/simulated_mean_fit_contrib
sd_log_PA[is.nan(sd_log_PA)] <- 0
sd_log_fit[is.nan(sd_log_fit)] <- 0

if ("PA" == which_plot) {
  y_dat      <- simulated_mean_PA_contrib
  y_upper    <- simulated_mean_PA_contrib + 2 * simulated_sd_PA_contrib 
  y_lower    <- simulated_mean_PA_contrib - 2 * simulated_sd_PA_contrib 
  y_observed <- original_result$contribution$PA_contribution 
} else {
  y_dat      <- simulated_mean_fit_contrib
  y_upper    <- simulated_mean_fit_contrib + 2 * simulated_sd_fit_contrib 
  y_lower    <- simulated_mean_fit_contrib - 2 * simulated_sd_fit_contrib 
  y_observed <- original_result$contribution$fit_contribution  
}

y_lower <- pmax(y_lower,0)

dat <- data.frame(x = 1:(T-1),y = y_dat, ymin = y_lower, ymax = y_upper,
                  y_obs = y_observed)

colors <- c("Calculated from simulations" = "#000000","Calculated from observed data" = "#e31a1c")
p<- ggplot(dat, aes(x=x, y=y)) + 
  scale_color_manual(values = colors) + 
  geom_point(alpha = 0.25,aes(colour = "Calculated from simulations")) + 
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.1, alpha = 0.25) +
  geom_point(aes(x = x, y = y_obs,colour = "Calculated from observed data"), alpha = 0.5,
             size = 1) +
  ylab(y_label) + xlab("Time-step") + ylim(0,max(y_upper,y_observed)) + 
  theme_bw() +
  theme(
    plot.margin =  margin(0.5, 0.25, 0.25, 0.25, "lines"), aspect.ratio=1,
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.title.y = element_text(angle=90,vjust =2),
    axis.title.x = element_text(vjust = -0.2),
    axis.text = element_text(), 
    axis.line = element_line(colour="black"),
    axis.ticks = element_line()) +
  theme(legend.position = c(legend_pos_x, legend_pos_y),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"))

plot(p)
}
