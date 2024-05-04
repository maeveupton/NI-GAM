# Clear workspace
rm(list = ls())

#---------Set working directory--------------
setwd('/users/research/mupton/1. NIGAM_paper_2022/2.Model_validation_separate_more_sites_proxy_only')
library(tidyverse)
library(gtools) # ordering the list of files properly


model_validation_tests <- function(model_type_path, save_name,plot_title){
  # Test data from Noisy part?
  test_data_files <- list.files(
    #path = "data/test_sets",
    path = "data/test_sets_noise",
                             pattern = "*.csv", 
                             full.names = T)
  test_data_list <- lapply(mixedsort(test_data_files),read_csv)
 
  # # Training data
  # training_data_files <- list.files(path = "data/training_sets",
  #                            pattern = "*.csv",
  #                            full.names = T)
  # training_data_list <- lapply(mixedsort(training_data_files),read_csv)
  
  # 10 Model runs
  model_runs_files <- list.files(
    path = model_type_path,
    #path = "output/no_noise_no_local",
    #path = "output/noise_full_model",
                              pattern = "*.rds", 
                              full.names = T) 
  model_runs_list <- lapply(mixedsort(model_runs_files), readRDS)


  # RMSE, MSE overall--------------------------------------
  # residual = observed - estimated
  model_residual_ME <- NA# Mean Error
  model_residual_MAE <- NA# Mean Absolute Error
  model_residual_RMSE <- NA # RMSE
  model_residual_site_run <- list() 
  sum_mod_res_site_run <- list()# By Site
  for(i in 1:10){
    model_residual_ME[i] <- mean(model_runs_list[[i]]$BUGSoutput$sims.list$residuals)
    model_residual_MAE[i] <- mean(abs(model_runs_list[[i]]$BUGSoutput$sims.list$residuals))
    model_residual_RMSE[i] <- sqrt(mean(model_runs_list[[i]]$BUGSoutput$sims.list$residuals)^2)
    # All sites
    model_residual_all <-as.data.frame(model_runs_list[[i]]$BUGSoutput$sims.list$residuals_pred)
    rownames(model_residual_all) <- paste0("sample_",1:nrow(model_residual_all))
    model_residual_all_long_df<- pivot_longer(
      data = as.data.frame(cbind(Age = test_data_list[[i]]$Age,
                                 SiteName = test_data_list[[i]]$SiteName,
                                 t(model_residual_all))),
                                            cols = starts_with("sample_"))
    model_residual_site_run[[i]]<- model_residual_all_long_df%>% 
      mutate(residual_pred = as.numeric(value)) %>% 
      as.data.frame()
    sum_mod_res_site_run[[i]] <- model_residual_site_run[[i]] %>% 
      group_by(SiteName) %>% 
      summarise(model_resid_ME_site = mean(residual_pred),
                model_resid_MAE_site = mean(abs(residual_pred),
                                            model_resid_RMSE = sqrt(mean(residual_pred)^2)))
  }

  # Plotting overall True vs predicted ----------------------------- 
  pred_true_df <- list()
  for(i in 1:10){
    pred_true_df[[i]] <- data.frame(
      true_RSL = test_data_list[[i]]$RSL,
      Age = test_data_list[[i]]$Age, 
      SiteName = test_data_list[[i]]$SiteName,
      Longitude = test_data_list[[i]]$Longitude,
      Latitude = test_data_list[[i]]$Latitude,
      pred = model_runs_list[[i]]$BUGSoutput$mean$mu_pred,
      pred_lwr = apply(model_runs_list[[i]]$BUGSoutput$sims.list$y_pred,2, 'quantile', 0.025),
      pred_upper = apply(model_runs_list[[i]]$BUGSoutput$sims.list$y_pred, 2, 'quantile',0.975),
      pred_lwr_50 = apply(model_runs_list[[i]]$BUGSoutput$sims.list$y_pred,2, 'quantile', 0.25),
      pred_upper_50 = apply(model_runs_list[[i]]$BUGSoutput$sims.list$y_pred, 2, 'quantile',0.75),
      run_id = as.character(i))
  
  }
  pred_true_full_df<-bind_rows(pred_true_df,.id = "column_label")
  # #----True vs Pred Plot for TG----
  # true_pred_plotTG<-ggplot()+
  #   geom_errorbar(data =pred_true_full_df, 
  #                 aes(x = true_RSL,ymin = pred_lwr,ymax = pred_upper),colour = "red3",
  #                 width=0,alpha = 0.5)+
  #   geom_point(data = pred_true_full_df, 
  #              aes(x = true_RSL, y = pred,colour="95% Prediction Interval"),
  #              fill="white",size = 0.5,shape = 0)+
  #   geom_abline(data = pred_true_full_df,
  #               aes(intercept = 0, slope = 1, colour ="True = Predicted"))+
  #   theme_bw()+
  #   theme(plot.title = element_text(size=22),
  #         axis.title=element_text(size=7,face="bold"),
  #         axis.text=element_text(size=7),
  #         strip.background =element_rect(fill=c("white")),
  #         strip.text = element_text(size = 4),
  #         legend.text=element_text(size=7),
  #         legend.title= element_blank(),
  #         legend.justification = c(1, 0),
  #         axis.text.x = element_text(size = 3),
  #         axis.text.y = element_text(size = 3),
  #         legend.position = c(1, 0))+
  #   labs(x = "True RSL",y = "Predicted RSL")+
  #   scale_colour_manual("", 
  #                       values = c("95% Prediction Interval"="red3",
  #                                  "True = Predicted" = "black"))+
  #   guides(colour = guide_legend(
  #     override.aes = list(pch = c(0,NA)),
  #     size = c(2,1),
  #     linetype = c(0,1))) +
  #                       # guide = guide_legend(
  #                       #   override.aes = list(pch = c(0,NA)),
  #                       #   size = c(100,1),
  #                       #   linetype = c(0,1))) +
  #   facet_wrap(~SiteName,scales = "free")
  # true_pred_plotTG
  # ggsave(true_pred_plotTG,file = paste0("fig/validations/true_pred",save_name,"TG.pdf"), 
  #        width = 12, height = 7)
  # 
  #----True vs Pred Plot for proxy----
  true_pred_plot_proxy<-ggplot()+
    geom_errorbar(data = pred_true_full_df, 
                  aes(x = true_RSL,ymin = pred_lwr,ymax = pred_upper),colour = "red3",
                  width=0,alpha = 0.5)+
    geom_point(data = pred_true_full_df, 
               aes(x = true_RSL, y = pred,colour="95% Prediction Interval"),
               fill="white",size = 0.5,shape = 0)+
    geom_abline(data = pred_true_full_df,
                aes(intercept = 0, slope = 1, colour ="True = Predicted"),alpha = 0.2)+
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          strip.background =element_rect(fill=c("white")),
          legend.text=element_text(size=12),
          legend.title= element_blank(),
          legend.justification = c(1, 0),
          legend.position = c(1, 0))+
    labs(x = "True RSL",y = "Predicted RSL")+
    scale_colour_manual("", 
                        values = c("95% Prediction Interval"="red3",
                                   "True = Predicted" = "black"))+
    guides(colour = guide_legend(
      override.aes = list(pch = c(0,NA)),
      size = c(2,1),
      linetype = c(0,1))) +
    # guide = guide_legend(
    #   override.aes = list(pch = c(0,NA)),
    #   size = c(100,1),
    #   linetype = c(0,1))) +
    facet_wrap(~SiteName,scales = "free")
  true_pred_plot_proxy
  ggsave(true_pred_plot_proxy,file = paste0("fig/validations/true_pred",save_name,"proxy.pdf"), width = 10, height = 6)
  
  # #----Plotting RSL data with predictions over time using TG----
  # RSL_pred_plotTG<-ggplot()+
  #   geom_errorbar(data = pred_true_full_df, 
  #                 aes(x = Age*1000,ymin = pred_lwr,
  #                     ymax = pred_upper, colour = "95% Prediction Interval"),
  #                 #colour = "red3",
  #                 width=0,alpha = 0.5)+
  #   geom_point(data = pred_true_full_df, aes(y = pred,x = Age*1000,
  #                                            colour = "Predicted Mean"),
  #              fill="white",shape = 0,alpha = 0.6,
  #              size = 1)+
  #   geom_point(data = pred_true_full_df, aes(y = true_RSL,
  #                                                    x = Age*1000,colour="True Mean"),
  #              shape = 20,
  #              size = 0.8)+
  #   theme_bw()+
  #   theme(plot.title = element_text(size=22),
  #         axis.title=element_text(size=14,face="bold"),
  #         axis.text=element_text(size=12),
  #         strip.background =element_rect(fill=c("white")),
  #         strip.text = element_text(size = 4),
  #         legend.text=element_text(size=12),
  #         axis.text.x = element_text(size = 3),
  #         axis.text.y = element_text(size = 3))+
  #   labs(x = "Year (CE)",y = "RSL (m)")+
  #   theme(legend.title= element_blank(),
  #         legend.justification = c(1, 0),
  #         legend.margin=margin(c(1,5,5,5)),
  #         legend.position = c(1, 0.1))+
  #   scale_colour_manual("", 
  #                       values = c("95% Prediction Interval"="red3",
  #                                  "Predicted Mean"="red3",
  #                                  "True Mean" = "gray13"),
  #                       guide = guide_legend(
  #                         override.aes = list(pch = c(NA, 0,20),size = c(0.5,2,2),
  #                                             linetype = c(1, 0, 0)))) +
  #   facet_wrap(~SiteName)
  # RSL_pred_plotTG
  # ggsave(RSL_pred_plotTG,file = paste0("fig/validations/CV_all_sites",save_name,"TG.pdf"),
  #        width = 12, height = 7)
  
  #----Plotting RSL data with predictions over time for proxy----
  RSL_pred_plot_proxy<-ggplot()+
    geom_errorbar(data = pred_true_full_df, 
                  aes(x = Age*1000,ymin = pred_lwr,
                      ymax = pred_upper, colour = "95% Prediction Interval"),
                  #colour = "red3",
                  width=0,alpha = 0.5)+
    geom_point(data =pred_true_full_df, aes(y = pred,x = Age*1000,
                                             colour = "Predicted Mean"),
               fill="white",shape = 0,alpha = 0.6,
               size = 1)+
    geom_point(data = pred_true_full_df, aes(y = true_RSL,
                                             x = Age*1000,colour="True Mean"),
               shape = 20,
               size = 0.8)+
    
    #ggtitle(paste0("Cross Validation for \n the 10 fold",plot_title))+
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          strip.background =element_rect(fill=c("white")),
          legend.text=element_text(size=12))+
    labs(x = "Year (CE)",y = "RSL (m)")+
    theme(legend.title= element_blank(),
          legend.justification = c(1, 0),
          legend.margin=margin(c(1,5,5,5)),
          legend.position = c(1, 0))+
    scale_colour_manual("", 
                        values = c("95% Prediction Interval"="red3",
                                   "Predicted Mean"="red3",
                                   "True Mean" = "gray13"),
                        guide = guide_legend(
                          override.aes = list(pch = c(NA, 0,20),size = c(0.5,2,2),
                                              linetype = c(1, 0, 0)))) +
    facet_wrap(~SiteName)
  RSL_pred_plot_proxy
  ggsave(RSL_pred_plot_proxy,file = paste0("fig/validations/CV_all_sites",save_name,"proxy.pdf"),
         width = 10, height = 6)
  
  
  # Calculating the overall coverage for 95% PI
  pred_true_full_df$obs_in_PI <- NA
  for(i in 1:nrow(pred_true_full_df)){
    pred_true_full_df$obs_in_PI[i] <- ifelse(between(pred_true_full_df$true_RSL[i],
                                                   pred_true_full_df$pred_lwr[i],
                                                   pred_true_full_df$pred_upper[i]),TRUE,FALSE)}

  # Calculating the overall coverage for 50% PI
  pred_true_full_df$obs_in_PI_50 <- NA
  for(i in 1:nrow(pred_true_full_df)){
    pred_true_full_df$obs_in_PI_50[i] <- ifelse(between(pred_true_full_df$true_RSL[i],
                                                     pred_true_full_df$pred_lwr_50[i],
                                                     pred_true_full_df$pred_upper_50[i]),TRUE,FALSE)}
  
  # Total coverage is trues/ number of rows with 95% PI
  total_coverage <- 
    length(which(pred_true_full_df$obs_in_PI == "TRUE"))/nrow(pred_true_full_df)
  total_coverage
  # Total coverage with 50% PI
  total_coverage_50 <- 
    length(which(pred_true_full_df$obs_in_PI_50 == "TRUE"))/nrow(pred_true_full_df)
  total_coverage_50
  # Coverage by site with 95% PI
  site_specific_coverage <- pred_true_full_df %>% 
    group_by(SiteName) %>% 
    summarise(coverage = length(which(obs_in_PI == "TRUE"))/n(),
              PI_width = mean(pred_upper-pred_lwr),
              coverage_50 = length(which(obs_in_PI_50 == "TRUE"))/n(),
              PI_width_50 = mean(pred_upper_50-pred_lwr_50 ))
  site_specific_coverage
  # Median of prediction and see if it above of below. It should be 50:50
  median_total <-pred_true_full_df %>% 
    mutate(median_pred  = median(pred_true_full_df$pred)) %>% 
    mutate(median_pred_condition = ifelse(true_RSL>=median_pred,1,0))
  median_total_test_upr <- median_total %>%
    summarise(above_median = sum(median_pred_condition)/n(),
              below_median = 1 - above_median)
   median_site <- pred_true_full_df %>% 
    group_by(SiteName) %>% 
    mutate(median_pred  = median(pred_true_full_df$pred)) %>% 
    mutate(median_pred_condition = ifelse(true_RSL>=median_pred,1,0)) %>% 
    summarise(above_median = sum(median_pred_condition)/n(),
              below_median = 1 - above_median)
  
  # 95% PI
  interval_size_df <- pred_true_full_df %>% 
    group_by(SiteName) %>% 
    summarise(PI_width = mean(pred_upper-pred_lwr ))
  # 50% PI
  interval_size_df_50 <- pred_true_full_df %>% 
    group_by(SiteName) %>% 
    summarise(PI_width_50 = mean(pred_upper_50-pred_lwr_50 ))
  
  # Saving all the MSE, RMSE & coverage values
  model_test_res <- save(total_coverage,
                         total_coverage_50,
                         median_site,
                         median_total,
                         median_total_test_upr,
                         site_specific_coverage,
                         interval_size_df,
                         interval_size_df_50,
                         model_residual_ME,
                         model_residual_site_run,
                         model_residual_MAE,
                         model_residual_RMSE,
                         sum_mod_res_site_run,
                         #model_residual_site,
                         pred_true_full_df,
                         file = paste0("data/model_val/model_validation_res",save_name,".RData"))
  
  # Plots for 4 sites only
  pred_true_full_4sites_df <- pred_true_full_df %>%
    filter(SiteName %in% c("East River Marsh,\n Connecticut",
                         "Swan Key,\n Florida",
                         "Placentia,\n Newfoundland",
                         "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  #------Ordering the sites-------
  all_data_sites<-factor(paste0(pred_true_full_4sites_df$Longitude, 
                                 pred_true_full_4sites_df$Longitude), labels = 1:4)#For ordering sites
  pred_true_full_4sites_df <- cbind(pred_true_full_4sites_df,all_data_sites=all_data_sites)
   order_sites <- pred_true_full_4sites_df %>%  group_by(SiteName,all_data_sites) %>%
     dplyr::summarise(n = n()) %>%
     arrange(all_data_sites)
  pred_true_full_4sites_df$SiteName <- factor(pred_true_full_4sites_df$SiteName,
                                              levels = unique(order_sites$SiteName))
  
  true_pred_4_sites_plot<-ggplot()+
    geom_errorbar(data = pred_true_full_4sites_df, 
                  aes(x = true_RSL,ymin = pred_lwr,ymax = pred_upper),colour = "red3",
                  width=0,alpha = 0.5)+
      geom_point(data = pred_true_full_4sites_df, 
               aes(x = true_RSL, y = pred,colour="95% Prediction Interval"),
               fill="white",size = 0.5,shape = 0)+
    geom_abline(data = pred_true_full_4sites_df,
                aes(intercept = 0, slope = 1, colour ="True = Predicted"),alpha = 0.2)+
    #ggtitle(paste0("True vs Predicted for \n the 10 fold",plot_title))+
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          strip.background =element_rect(fill=c("white")),
          legend.text=element_text(size=9),
          legend.title= element_blank(),
          legend.margin=margin(c(1,5,5,5)),
          legend.justification = c(1, 0),
          legend.position = c(1, 0.002))+
    scale_colour_manual("", 
                        values = c("95% Prediction Interval"="red3",
                                   "True = Predicted" = "black"),
                        guide = guide_legend(
                          override.aes = list(pch = c(0,NA)),
                          linetype = c("blank", "dashed"))) +
    labs(x = "True RSL",y = "Predicted RSL")+
    facet_wrap(~SiteName)
  true_pred_4_sites_plot
  ggsave(true_pred_4_sites_plot,file = paste0("fig/validations/true_pred_4_sites",save_name,".pdf"), width = 10, height = 6)
  #ggsave(true_pred_4_sites_plot,file = paste0("fig/validations/true_pred_all_sites",save_name,".pdf"), width = 10, height = 6)


  # Plotting RSL data with predictions over time
  RSL_pred_plot<-ggplot()+
    geom_errorbar(data = pred_true_full_4sites_df, 
                  aes(x = Age*1000,ymin = pred_lwr,
                      ymax = pred_upper, colour = "95% Prediction Interval"),
                  #colour = "red3",
                  width=0,alpha = 0.5)+
    geom_point(data = pred_true_full_4sites_df, aes(y = pred,x = Age*1000,
                                             colour = "Predicted Mean"),
               fill="white",shape = 0,alpha = 0.6,
               size = 1)+
    geom_point(data = pred_true_full_4sites_df, aes(y = true_RSL,
                                             x = Age*1000,colour="True Mean"),
               shape = 20,
               size = 1)+
    
    #ggtitle(paste0("Cross Validation for \n the 10 fold",plot_title))+
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          strip.background =element_rect(fill=c("white")),
          legend.text=element_text(size=12))+
    labs(x = "Year (CE)",y = "RSL (m)")+
    theme(legend.title= element_blank(),
          legend.margin=margin(c(1,5,5,5)),
          legend.justification = c(1, 0),
          legend.position = c(0.95, 0.005))+
    scale_colour_manual("", 
                        values = c("95% Prediction Interval"="red3",
                                   "Predicted Mean"="red3",
                                   "True Mean" = "gray13"),
                        guide = guide_legend(
                          override.aes = list(pch = c(NA, 0,20),size = c(0.5,3,3),
                                              linetype = c(1, 0, 0)))) +
    facet_wrap(~SiteName)
  RSL_pred_plot
  ggsave(RSL_pred_plot,file =
           paste0("fig/validations/CV_4_sites",save_name,".pdf"),
         width = 10, height = 6)
  # ggsave(RSL_pred_plot,file = 
  #          paste0("fig/validations/CV_all_sites",save_name,".pdf"), 
  #        width = 10, height = 6)
}

full_model_validation <-
  model_validation_tests(
    #model_type_path="output/no_noise_no_local",
    model_type_path="output/noise_full_model",
                         save_name = "_full_model",
                         plot_title = " Full Model with Noise")

load("data/model_val/model_validation_res_full_model.RData")

#-- Overall & Site Specific MAE & RMSE----
model_residual_site_run_df<-
  bind_rows(model_residual_site_run,.id = "column_label")# bind runs for means

model_residual_overall_summary <-model_residual_site_run_df %>% 
  summarise(model_resid_ME = mean(residual_pred),
            model_resid_MAE = mean(abs(residual_pred)),
            model_resid_RMSE_overall = sqrt(mean(residual_pred^2)))

model_residual_site_summary <-model_residual_site_run_df %>% group_by(SiteName) %>% 
  summarise(model_resid_ME_site = mean(residual_pred),
            model_resid_MAE_site = mean(abs(residual_pred)),
            model_resid_RMSE = sqrt(mean(residual_pred^2)))
#----Total coverage----
total_coverage
total_coverage_50

#--Writing results in table for Latex
library(xtable)
print(xtable(model_residual_site_summary,digits=6,type = "latex"),
      file = "latex_output/model_resid_bias_tests_proxy.tex",include.rownames=FALSE)

print(xtable(site_specific_coverage,type = "latex"),
      file = "latex_output/site_sepecific_coverage_proxy.tex",include.rownames=FALSE)


model_residual_site_summary_4 <-
  model_residual_site_summary %>%
  filter(SiteName %in% c("Placentia,\n Newfoundland","East River Marsh,\n Connecticut",
                         "Cedar Island,\n North Carolina","Swan Key,\n Florida"))
print(xtable(model_residual_site_summary_4,digits=6,type = "latex"),
      file = "latex_output/model_resid_bias_tests_4_sites.tex",include.rownames=FALSE)


site_specific_coverage_4_sites <-
  site_specific_coverage %>% filter(SiteName %in% c("Placentia,\n Newfoundland",
                                                    "East River Marsh,\n Connecticut",
                                                    "Cedar Island,\n North Carolina",
                                                    "Swan Key,\n Florida"))
print(xtable(site_specific_coverage_4_sites,type = "latex"),
      file = "latex_output/site_sepecific_coverage_4_sites.tex",
      include.rownames=FALSE)

