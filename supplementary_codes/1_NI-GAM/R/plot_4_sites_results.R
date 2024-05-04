plot_4_sites_results <- function(SL_df,model_run,basis_fun_list,
                                 save_option,set_up_option){
  
  dir.create(paste0("fig/",save_option),showWarnings = F)
  
  
  #-----------Plotting Overall Output----------
  mu_post <- model_run$BUGSoutput$sims.list$mu
  #---Get estimates and uncertainty bounds--
  RSL_mod<-apply(mu_post,2,mean)
  RSL_mod_upr<-apply(mu_post,2,quantile,probs=0.025)
  RSL_mod_lwr<-apply(mu_post,2,quantile,probs=0.975)
  lwr_50<- apply(mu_post,2,quantile,probs=0.25)
  upr_50<- apply(mu_post,2,quantile,probs=0.75)
  
  #-----Create data frame for plotting for total-----
  mod_output_df<-data.frame(RSL_mod,RSL_mod_upr,RSL_mod_lwr,lwr_50,upr_50,SL_df$Age,
                            SL_df$SiteName,SL_df$Longitude,SL_df$Latitude,
                            ID = "Total Posterior Model")
  names(mod_output_df)<-c("RSL","upr","lwr","upr_50","lwr_50", "Age",
                          "SiteName","Longitude","Latitude",
                          "ID")
  # Filter
  mod_output_df <- mod_output_df %>% 
    filter(SiteName %in% c(#"ARGENTIA",
      "East River Marsh,\n Connecticut",
                           "Swan Key,\n Florida",
                           "Placentia,\n Newfoundland",
                           "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  #------Ordering the sites total component-------
  all_data_sites<-factor(paste0(mod_output_df$Longitude,
                                mod_output_df$Longitude), labels = 1:4)#For ordering sites
  mod_output_df <- cbind(mod_output_df,all_data_sites=all_data_sites)
  order_sites <- mod_output_df %>%  group_by(SiteName,all_data_sites) %>%
    dplyr::summarise(n = n()) %>%
    arrange(all_data_sites)
  mod_output_df$SiteName <- factor(mod_output_df$SiteName,
                                          levels = unique(order_sites$SiteName))
  
  
  #------------Regional Component: Spline in Time--------------
  time_component <- model_run$BUGSoutput$sims.list$regional
  RSL_mod_NI<-apply(time_component,2,mean)
  RSL_mod_lwr_NI<-apply(time_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(time_component,2,quantile,probs=0.975)
  lwr_50<-apply(time_component,2,quantile,probs=0.25)
  upr_50<-apply(time_component,2,quantile,probs=0.75)
  
  time_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
                                lwr_50,upr_50,
                                Age = SL_df$Age,
                                SiteName = SL_df$SiteName,
                                ID = "Regional Component")
  names(time_component_df)<-c("RSL","upr","lwr","lwr_50","upr_50", "Age","SiteName","ID")
  
  time_component_df <- time_component_df %>%
    filter(SiteName %in% c(#"ARGENTIA",
      "East River Marsh,\n Connecticut",
                           "Swan Key,\n Florida",
                           "Placentia,\n Newfoundland",
                           "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  
  #---GIA---
  # GIA_component <- model_run$BUGSoutput$sims.list$GIA_correction
  # RSL_mod_NI<-apply(GIA_component,2,mean)
  # RSL_mod_lwr_NI<-apply(GIA_component,2,quantile,probs=0.025)
  # RSL_mod_upr_NI<-apply(GIA_component,2,quantile,probs=0.975)
  # GIA_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
  #                              Age = SL_df$Age,
  #                              SiteName = SL_df$SiteName,
  #                              ID = "Linear Local Component")
  # names(GIA_component_df)<-c("RSL","upr","lwr","Age","SiteName","ID")
  # 
  # GIA_component_df <- GIA_component_df %>%  
  #   filter(SiteName %in% c(#"ARGENTIA",
  #     "East River Marsh,\n Connecticut",
  #                          "Swan Key,\n Florida",
  #                          "Placentia,\n Newfoundland",
  #                          "Cedar Island,\n North Carolina")) %>%
  #   mutate(SiteName = as.factor(SiteName))
  
  #---GIA + yoffset---
  GIA_yoffset_component <- model_run$BUGSoutput$sims.list$GIA_yoffset
  RSL_mod_NI<-apply(GIA_yoffset_component,2,mean)
  RSL_mod_lwr_NI<-apply(GIA_yoffset_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(GIA_yoffset_component,2,quantile,probs=0.975)
  GIA_yoffset_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
                               Age = SL_df$Age,
                               SiteName = SL_df$SiteName,
                               ID = "Linear Local Component \n and site-specific vertical offset")
  names(GIA_yoffset_component_df)<-c("RSL","upr","lwr","Age","SiteName","ID")

  GIA_yoffset_component_df <- GIA_yoffset_component_df %>%
    filter(SiteName %in% c(#"ARGENTIA",
      "East River Marsh,\n Connecticut",
                           "Swan Key,\n Florida",
                           "Placentia,\n Newfoundland",
                           "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  
  #------Non-Linear Local Component: Spline in Space Time------
  space_time_component <- model_run$BUGSoutput$sims.list$local
  RSL_mod_NI<-apply(space_time_component,2,mean)
  RSL_mod_lwr_NI<-apply(space_time_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(space_time_component,2,quantile,probs=0.975)
  lwr_50<-apply(space_time_component,2,quantile,probs=0.25)
  upr_50<-apply(space_time_component,2,quantile,probs=0.75)
  space_time_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
                                      lwr_50,upr_50,
                                      Age = SL_df$Age,
                                      Longitude = SL_df$Longitude,
                                      Latitude = SL_df$Latitude,
                                      SiteName = SL_df$SiteName,
                                      ID = "Non-Linear Local Component")
  names(space_time_component_df)<-c("RSL","upr", "lwr","lwr_50","upr_50", "Age",
                                    "Longitude","Latitude","SiteName","ID")
  
  space_time_component_df <- space_time_component_df %>%
    filter(SiteName %in% c(#"ARGENTIA",
      "East River Marsh,\n Connecticut",
                           "Swan Key,\n Florida",
                           "Placentia,\n Newfoundland",
                           "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  
  
  # #--Filtering for 4 sites---
  SL_df <-SL_df %>%   filter(SiteName %in% c(#"ARGENTIA",
    "East River Marsh,\n Connecticut",
                                             "Swan Key,\n Florida",
                                             "Placentia,\n Newfoundland",
                                             "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  #------Ordering the sites original Data set-------
  all_data_sites<-factor(paste0(SL_df$Longitude, SL_df$Longitude), labels = 1:4)#For ordering sites
  SL_df <- cbind(SL_df,all_data_sites=all_data_sites)
  order_sites <- SL_df %>%  group_by(SiteName,all_data_sites) %>%
    dplyr::summarise(n = n()) %>%
    arrange(all_data_sites)
  SL_df$SiteName <- factor(SL_df$SiteName,levels = unique(order_sites$SiteName))

  
  #------Ordering the sites Space time component-------
  # all_data_sites<-factor(paste0(space_time_component_df$Longitude,
  #                               space_time_component_df$Longitude), labels = 1:4)#For ordering sites
  # space_time_component_df <- cbind(space_time_component_df,all_data_sites=all_data_sites)
  # order_sites <- space_time_component_df %>%  group_by(SiteName,all_data_sites) %>%
  #   dplyr::summarise(n = n()) %>%
  #   arrange(all_data_sites)
  space_time_component_df$SiteName <- factor(SL_df$SiteName,levels = unique(order_sites$SiteName))
  
  #---Total plot----
  plot_result <- ggplot() +
    facet_wrap(~SiteName)+
    geom_point(data = SL_df, aes(y = RSL, x =  Age*1000), size = 0.5) +
    geom_rect(data = SL_df, aes(
      xmin = Age*1000 - Age_er_average*1000, xmax = Age*1000 + Age_er_average*1000,
      ymin = RSL - RSL_er_average, ymax = RSL + RSL_er_average
    ), alpha = 0.4) +
    geom_line(data = mod_output_df, aes(x = Age*1000, y = RSL), colour = "purple3") +
    geom_ribbon(data=mod_output_df,
                aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000),fill = "purple3",alpha=0.2)+
    geom_ribbon(data=mod_output_df,
                aes(y = RSL,ymin=lwr_50,ymax=upr_50,x=Age*1000),fill = "purple3",alpha=0.3)+
        xlab("Year (CE)") +
    ylab("RSL (m)") +
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    theme(strip.text.x = element_text(size = 14),
          strip.background =element_rect(fill=c("white")))
    #ggtitle("Total Model Fit for East Coast of North America")
    #ggtitle(paste0("Total Model Fit for East Coast of North America"," ",set_up_option, sep = " "))
  plot_result
  ggsave(paste0("fig/",save_option,"/",
                unique(mod_output_df$ID), ".pdf", sep = ""),
         plot_result, width = 10, height = 6)
  #ggsave(plot_result, file = "fig/total_model_fit.pdf", width = 10, height = 6)


  #---Regional Plot----
  regional_plot<-ggplot()+
    geom_line(data=time_component_df,aes(x=Age*1000,y=RSL),colour="#3b47ad")+
    geom_ribbon(data=time_component_df,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#3b47ad",alpha=0.2)+
    geom_ribbon(data=time_component_df,aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="#3b47ad",alpha=0.3)+
    #ggtitle("Regional Component")+
    #ggtitle(paste0("Regional Component",set_up_option, sep = " "))+
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    ylab('Sea Level (m)')+
    xlab('Year(CE)')
  regional_plot
  ggsave(paste0("fig/",save_option,'/' ,
                unique(time_component_df$ID), ".pdf", sep = ""),
         regional_plot, width = 10, height = 6)
  #ggsave(regional_plot,file = "fig/regional_component.pdf", width = 10, height = 6)
  
  cat("Regional Component plotted")
  
  # #------Linear Local Component: Random Effect------
  # GIA_plot<-ggplot()+
  #   geom_line(data=GIA_component_df,aes(x=Age*1000,y=RSL),colour="#5bac06")+
  #   geom_ribbon(data=GIA_component_df,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#5bac06",alpha=0.3)+
  #   theme_bw()+
  #   #ggtitle("Linear Local Component")+
  #   #ggtitle(paste0("Linear Local Component",set_up_option, sep = " "))+
  #   theme(plot.title = element_text(size=22),
  #         axis.title=element_text(size=14,face="bold"),
  #         axis.text=element_text(size=12),
  #         legend.text=element_text(size=12))+
  #   theme(strip.text.x = element_text(size = 14),
  #         strip.background =element_rect(fill=c("white")))+
  #   ylab('Sea Level (m)')+
  #   facet_wrap(~SiteName) +
  #   xlab('Year(CE)')
  # GIA_plot
  # ggsave(paste0("fig/",save_option,"/",
  #               unique(GIA_component_df$ID), ".pdf", sep = ""),
  #        GIA_plot, width = 10, height = 6)
  # #ggsave(GIA_plot,file = "fig/GIA_component.pdf", width = 10, height = 6)
  # 
  # cat("Linear Local Component plotted")
  
  #------Linear Local Component + y offset: Random Effect------
  # GIA_yoffset_plot<-ggplot()+
  #   geom_line(data=GIA_yoffset_component_df,aes(x=Age*1000,y=RSL),colour="#5bac06")+
  #   geom_ribbon(data=GIA_yoffset_component_df,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#5bac06",alpha=0.3)+
  #   theme_bw()+
  #     theme(plot.title = element_text(size=22),
  #         axis.title=element_text(size=14,face="bold"),
  #         axis.text=element_text(size=12),
  #         legend.text=element_text(size=12))+
  #   theme(strip.text.x = element_text(size = 14),
  #         strip.background =element_rect(fill=c("white")))+
  #   ylab('Sea Level (m)')+
  #   facet_wrap(~SiteName) +
  #   xlab('Year(CE)')
  # GIA_yoffset_plot
  # ggsave(paste0("fig/",save_option,"/",
  #               unique(GIA_yoffset_component_df$ID), ".pdf", sep = ""),
  #        GIA_yoffset_plot, width = 10, height = 6)
  # #ggsave(GIA_plot,file = "fig/GIA_component.pdf", width = 10, height = 6)
  
  cat("GIA + y offset Component plotted")
  
  #----Local Plot---  
  local_plot<-ggplot()+
    geom_line(data=space_time_component_df,aes(x=Age*1000,y=RSL),colour="#ad4c14")+
    geom_ribbon(data=space_time_component_df,
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="#ad4c14",alpha=0.2)+
    geom_ribbon(data=space_time_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#ad4c14",alpha=0.3)+
    geom_hline(yintercept = 0)+
    theme_bw()+
    ylab('Sea Level (m)')+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    facet_wrap(~SiteName) +
    theme(strip.text.x = element_text(size = 14),
          strip.background =element_rect(fill=c("white")))+
    xlab('Year(CE)')
  local_plot
  #ggsave(local_plot,file = "fig/local_component.pdf", width = 10, height = 6)
  ggsave(paste0("fig/",save_option,"/",
                unique(space_time_component_df$ID), ".pdf", sep = ""),
         local_plot, width = 10, height = 6)
  cat("Non-Linear Local Component plotted")
  
  #-------y offset------
  # y_offset_post <-model_run$BUGSoutput$sims.list$intercept
  # RSL_mod_NI<-apply(y_offset_post,2,mean)
  # RSL_mod_lwr_NI<-apply(y_offset_post,2,quantile,probs=0.025)
  # RSL_mod_upr_NI<-apply(y_offset_post,2,quantile,probs=0.975)
  # y_offset_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
  #                                   Age = SL_df$Age,
  #                                   SiteName = SL_df$SiteName,
  #                                   ID = "y offset Component")
  # names(y_offset_component_df)<-c("RSL","upr","lwr","Age","SiteName","ID")
  
  
  #----Separate Components on one plot with CI----
  all_components_CI_plot <- ggplot()+
    # # all components & The data
    # geom_rect(data = SL_df, aes(
    #   xmin = Age*1000 - Age_er_average*1000, xmax = Age*1000 + Age_er_average*1000,
    #   ymin = RSL - RSL_er_average, ymax = RSL + RSL_er_average#,fill = "Observed Uncertainty"
    # ), alpha = 0.2) +
    # geom_point(data = SL_df, aes(y = RSL, x = Age*1000), size = 0.5) +
    # #guides(color = guide_legend(override.aes = list(size=2)))+
    
    # Local
    geom_line(data=space_time_component_df,
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=space_time_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
    
    # Linear Local Component + y offset
    geom_line(data=GIA_yoffset_component_df,
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=GIA_yoffset_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill= ID),alpha=0.3)+
    # Global Component
    geom_line(data=time_component_df,
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=time_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
    #Total Model
    geom_line(data=mod_output_df,
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=mod_output_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
    #scale_fill_manual('',values = c("#5bac06","dodgerblue3","darkorange4","purple3"))+
    scale_fill_manual(name = '',values = c("#5bac06","#ad4c14","#3b47ad","purple3"),
                      guide = guide_legend(override.aes = list(alpha = 0.1)))+
    #guides(fill = guide_legend(override.aes = list(shape = c(19,19,19,19,NA,NA),
    #                                               linetype = c(0,0,0,0, 0,1))))+
    scale_colour_manual(name = '',values = c("#5bac06","#ad4c14","#3b47ad","purple3"))+
    # scale_fill_manual('',
    #                    values = 'grey',  
    #                    guide = guide_legend(override.aes = list(alpha = 0.5))) +
    theme_bw()+
    theme(plot.title = element_text(size=18),
          axis.title=element_text(size=12,face="bold"),
          axis.text=element_text(size=8),
          legend.text=element_text(size=10))+
    theme(strip.text.x = element_text(size = 14),
          strip.background =element_rect(fill=c("white")))+
    theme(legend.position="bottom", legend.box = "horizontal")+
    guides(color = guide_legend(override.aes = list(size = 0.5)))+
    ylab('RSL (m)')+
    facet_wrap(~SiteName) +
    xlab(' Year(CE)')
  
  all_components_CI_plot
  #ggsave(all_components_CI_plot,file = "fig/all_components_one_plot.pdf", width = 10, height = 6)
  ggsave(paste0("fig/",save_option,"/","all_component_plot.pdf", sep = ""),
         all_components_CI_plot, width = 10, height = 6)
  
  
  # #------------Rate of change-----------------
  # SL_df_separate_allTG <- read_csv("data/result_csv/regional_rate.csv")
  # SL_df_separate_allTG <- SL_df_separate_allTG %>%
  #   filter(SiteName %in% c("East River Marsh,\n Connecticut",
  #                          "Swan Key,\n Florida",
  #                          "Placentia,\n Newfoundland",
  #                          "Cedar Island,\n North Carolina")) %>%
  #   mutate(SiteName = as.factor(SiteName))
  # 
  # #---Regional Rate Plot----
  # regional_rate_plot<-ggplot()+
  #   geom_line(data=SL_df_separate_allTG,aes(x=Age*1000,y=RSL),colour="#3b47ad")+
  #   geom_ribbon(data=SL_df_separate_allTG,
  #               aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#3b47ad",alpha=0.2)+
  #   geom_ribbon(data=SL_df_separate_allTG,
  #               aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="#3b47ad",alpha=0.3)+
  # #  ggtitle("Regional Rate Plot comparison")+
  #   theme_bw()+
  #   theme(plot.title = element_text(size=22),
  #         axis.title=element_text(size=14,face="bold"),
  #         axis.text=element_text(size=12),
  #         legend.text=element_text(size=12))+
  #   ylab('Rate of Change (mm/yr)')+
  #   #facet_wrap(~SiteName)+
  #   xlab('Year(CE)')
  # regional_rate_plot
  # ggsave(paste0("fig/",save_option,"/","regional_rate_plot.pdf", sep = ""),
  #        regional_rate_plot, width = 10, height = 6)
  
  
  
  # #-----------Total Rate of Change--------------
  # SL_df_separate_allTG <- read_csv("data/result_csv/total_rate.csv")
  # SL_df_separate_allTG <- SL_df_separate_allTG %>%
  #   filter(SiteName %in% c("East River Marsh,\n Connecticut",
  #                          "Swan Key,\n Florida",
  #                          "Placentia,\n Newfoundland",
  #                          "Cedar Island,\n North Carolina")) %>%
  #   mutate(SiteName = as.factor(SiteName))
  # SL_df_separate_allTG$SiteName <- factor(SL_df$SiteName,levels = unique(order_sites$SiteName))
  # #---Total Rate Plot----
  # total_rate_plot<-ggplot()+
  #   geom_line(data=SL_df_separate_allTG,aes(x=Age*1000,y=RSL),colour="purple3")+
  #   geom_ribbon(data=SL_df_separate_allTG,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="purple3",alpha=0.2)+
  #   geom_ribbon(data=SL_df_separate_allTG,aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="purple3",alpha=0.3)+
  #   #ggtitle("Total Rate Plot comparison")+
  #   theme_bw()+
  #   theme(plot.title = element_text(size=22),
  #         axis.title=element_text(size=14,face="bold"),
  #         axis.text=element_text(size=12),
  #         legend.text=element_text(size=12))+
  #   ylab('Rate of Change (mm/yr)')+
  #   facet_wrap(~SiteName)+
  #   xlab('Year (CE)')+
  #   theme(strip.text.x = element_text(size = 14),
  #         strip.background =element_rect(fill=c("white")))
  # total_rate_plot
  # ggsave(paste0("fig/",save_option,"/","total_rate_plot.pdf", sep = ""),
  #        total_rate_plot, width = 10, height = 6)
  # 
  
  cat("All Components plotted")
}