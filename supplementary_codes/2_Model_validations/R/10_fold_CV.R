# Clear workspace
rm(list = ls())

#---------Set working directory--------------
setwd('/users/research/mupton/Splines/NI_GAM_RSL_Sep22/0.NIGAM_east_coast_fixed_regional_update_rates_update_model_validations')
library(tidyverse)
library(gtools) # ordering the list of files properly

# Test data 
test_data_files <- list.files(path = "data/test_sets",
                             pattern = "*.csv", 
                             full.names = T) #%>% map_df(~read_csv(.))
test_data_list <- lapply(mixedsort(test_data_files),read_csv)
#test<-do.call("cbind",test_data_list)
# Training data
training_data_files <- list.files(path = "data/training_sets",
                             pattern = "*.csv", 
                             full.names = T)# %>% map_df(~read_csv(.))
training_data_list <- lapply(mixedsort(training_data_files),read_csv)
# 10 Model runs
model_runs_files <- list.files(
  #path = "output/no_noise_no_local",
  path = "output/noise_full_model",
                              pattern = "*.rds", 
                              full.names = T) 
model_runs_list <- lapply(mixedsort(model_runs_files), readRDS)


# RMSE, MSE --------------------------------------
# residual = observed - estimated
model_residual_MSE <- NA# MSE
model_residual_RMSE <- NA # RMSE
for(i in 1:10){
  model_residual_MSE[i] <- mean(model_runs_list[[i]]$BUGSoutput$sims.list$residuals)^2
  model_residual_RMSE[i] <- sqrt(mean(model_runs_list[[i]]$BUGSoutput$sims.list$residuals)^2)
}

# Plotting overall True vs predicted for non local plot----------------------------- 
pred_true_df <- list()
for(i in 1:10){
  pred_true_df[[i]] <- data.frame(
    true_RSL = test_data_list[[i]]$RSL,
    Age = test_data_list[[i]]$Age, 
    SiteName = test_data_list[[i]]$SiteName,
    pred = model_runs_list[[i]]$BUGSoutput$mean$mu_pred,
    pred_lwr = apply(model_runs_list[[i]]$BUGSoutput$sims.list$y_pred,2, 'quantile', 0.025),
    pred_upper = apply(model_runs_list[[i]]$BUGSoutput$sims.list$y_pred, 2, 'quantile',0.975),
    run_id = as.character(i))
  
}
pred_true_full_df<-bind_rows(pred_true_df,.id = "column_label")#do.call("bind",pred_true_df)

true_pred_plot<-ggplot()+
  geom_errorbar(data = pred_true_full_df, 
                aes(x = true_RSL,ymin = pred_lwr,ymax = pred_upper),width=0)+
  #colour=run_id)
  geom_point (data = pred_true_full_df, aes(x = true_RSL, y = pred),fill="white",shape = 21)+
                                            #colour = run_id))+

  #ggtitle("True vs Predicted for \n the 10 fold Model without the Local Component")+
  ggtitle("True vs Predicted for \n the 10 fold full Model with Noise")+
  theme_bw()+
  theme(plot.title = element_text(size=22),
        axis.title=element_text(size=14,face="bold"),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12))+
  labs(x = "True RSL",y = "Predicted RSL")+
  geom_abline()#+
#facet_wrap(~SiteName)
true_pred_plot
#ggsave(true_pred_plot,file = "fig/validations/true_pred_full_model_site_10_fold.pdf", width = 10, height = 6)
ggsave(true_pred_plot,file = "fig/validations/true_pred_full_model_10_fold.pdf", width = 10, height = 6)
#ggsave(true_pred_plot,file = "fig/validations/true_pred_no_local_site_10_fold.pdf", width = 10, height = 6)
#ggsave(true_pred_plot,file = "fig/validations/true_pred_no_local_10_fold.pdf", width = 10, height = 6)

# Calculating the overall coverage
pred_true_full_df$obs_in_PI <- NA
for(i in 1:nrow(pred_true_full_df)){
  pred_true_full_df$obs_in_PI[i] <- ifelse(between(pred_true_full_df$true_RSL[i],
                                                   pred_true_full_df$pred_lwr[i],
                                                   pred_true_full_df$pred_upper[i]),TRUE,FALSE)}

# Total coverage is trues/ number of rows
total_coverage <- length(which(pred_true_full_df$obs_in_PI == "TRUE"))/nrow(pred_true_full_df)
total_coverage
site_specific_coverage <- pred_true_full_df %>% 
  group_by(SiteName) %>% 
  summarise(coverage = length(which(obs_in_PI == "TRUE"))/n())
site_specific_coverage
run_specific_coverage <- pred_true_full_df %>% 
  group_by(run_id) %>% 
  summarise(coverage = length(which(obs_in_PI == "TRUE"))/n())


# Filter
pred_true_full_4sites_df <- pred_true_full_df %>% 
  filter(SiteName %in% c("East River Marsh,\n Connecticut",
                         "Little Manatee River,\n Florida",
                         "Big River Marsh,\n Newfoundland",
                         "Cedar Island,\n North Carolina")) %>%
  mutate(SiteName = as.factor(SiteName))
true_pred_4_sites_plot<-ggplot()+
  geom_errorbar(data = pred_true_full_4sites_df, 
                aes(x = true_RSL,ymin = pred_lwr,
                    ymax = pred_upper),width=0)+
                    #colour=run_id)
  geom_point (data = pred_true_full_4sites_df, aes(x = true_RSL,
                                                   y = pred),fill="white",shape = 21)+
  #colour = run_id))+
  #ggtitle("True vs Predicted for \n the 10 fold Model without the Local Component")+
  ggtitle("True vs Predicted for \n the 10 fold full Model with Noise")+
  theme_bw()+
  theme(plot.title = element_text(size=22),
        axis.title=element_text(size=14,face="bold"),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12))+
  labs(x = "True RSL (m)",y = "Predicted RSL (m)")+
  geom_abline(linetype = "dashed")+
facet_wrap(~SiteName)
true_pred_4_sites_plot
ggsave(true_pred_4_sites_plot,file = "fig/validations/true_pred_4_sites_full_model_site_10_fold.pdf", width = 10, height = 6)
#ggsave(true_pred_4_sites_plot,file = "fig/validations/true_pred_4_sites no_local_site_10_fold.pdf", width = 10, height = 6)

# Plotting RSL data with predictions over time
RSL_pred_plot<-ggplot()+
  #colour=run_id)
  geom_point (data = pred_true_full_4sites_df, aes(y = pred,
                                                   x = Age*1000),fill="white",shape = 0)+
  geom_point (data = pred_true_full_4sites_df, aes(y = true_RSL,
                                                   x = Age*1000))+
  geom_errorbar(data = pred_true_full_4sites_df, 
                aes(x = Age*1000,ymin = pred_lwr,
                    ymax = pred_upper),width=0)+
  #colour = run_id))+
  #ggtitle("True vs Predicted for \n the 10 fold Model without the Local Component")+
  ggtitle("True vs Predicted for \n the 10 fold full Model with Noise")+
  theme_bw()+
  theme(plot.title = element_text(size=22),
        axis.title=element_text(size=14,face="bold"),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12))+
  labs(x = "Year (CE)",y = "RSL (m)")+
  facet_wrap(~SiteName)
RSL_pred_plot
ggsave(RSL_pred_plot,file = "fig/validations/cv_full_model_4site.pdf", width = 10, height = 6)
#ggsave(RSL_pred_plot,file = "fig/validations/cv_no_local_4site.pdf", width = 10, height = 6)

