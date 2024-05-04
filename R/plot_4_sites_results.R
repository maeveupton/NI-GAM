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
  mod_output_df<-data.frame(RSL_mod,RSL_mod_upr,
                            RSL_mod_lwr,
                            lwr_50,upr_50,
                            SL_df$Age,
                            SL_df$SiteName,SL_df$Longitude,SL_df$Latitude,
                            ID = "Total Posterior Model")
  names(mod_output_df)<-c("RSL","upr","lwr","upr_50","lwr_50", "Age",
                          "SiteName","Longitude","Latitude",
                          "ID")
  # Filter for 4 sites 
  mod_output_df <- mod_output_df %>% 
    filter(SiteName %in% c(
      "East River Marsh,\n Connecticut",
      "Swan Key,\n Florida",
      "Placentia,\n Newfoundland",
      "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  
  #------Ordering the sites total component-------
  all_data_sites<-factor(paste0(mod_output_df$Longitude,
                                mod_output_df$Longitude), labels = 1:4)
  mod_output_df <- cbind(mod_output_df,all_data_sites=all_data_sites)
  order_sites <- mod_output_df %>%  group_by(SiteName,all_data_sites) %>%
    dplyr::summarise(n = n()) %>%
    arrange(all_data_sites)
  mod_output_df$SiteName <- factor(mod_output_df$SiteName,
                                   levels = unique(order_sites$SiteName))
  
  
  #------------Regional Component: Spline in Time--------------
  time_component <- model_run$BUGSoutput$sims.list$r
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
  names(time_component_df)<-c("RSL","upr","lwr","lwr_50",
                              "upr_50", "Age","SiteName","ID")
  
  time_component_df <- time_component_df %>%
    filter(SiteName %in% c(
      "East River Marsh,\n Connecticut",
      "Swan Key,\n Florida",
      "Placentia,\n Newfoundland",
      "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  
  #---Linear Local + Site specific vertical offset---
  g_h_component <- model_run$BUGSoutput$sims.list$g_h_z_x
  RSL_mod_NI<-apply(g_h_component,2,mean)
  RSL_mod_lwr_NI<-apply(g_h_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(g_h_component,2,quantile,probs=0.975)
  g_h_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
                               Age = SL_df$Age,
                               SiteName = SL_df$SiteName,
                               ID = "Linear Local Component \n and site-specific vertical offset")
  names(g_h_component_df)<-c("RSL","upr","lwr","Age","SiteName","ID")
  
  g_h_component_df <- g_h_component_df %>%
    filter(SiteName %in% c(
      "East River Marsh,\n Connecticut",
      "Swan Key,\n Florida",
      "Placentia,\n Newfoundland",
      "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  
  #------Non-Linear Local Component: Spline in Space Time------
  space_time_component <- model_run$BUGSoutput$sims.list$l
  RSL_mod_NI<-apply(space_time_component,2,mean)
  RSL_mod_lwr_NI<-apply(space_time_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(space_time_component,2,quantile,probs=0.975)
  lwr_50<-apply(space_time_component,2,quantile,probs=0.25)
  upr_50<-apply(space_time_component,2,quantile,probs=0.75)
  space_time_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,
                                      RSL_mod_lwr_NI,
                                      lwr_50,upr_50,
                                      Age = SL_df$Age,
                                      Longitude = SL_df$Longitude,
                                      Latitude = SL_df$Latitude,
                                      SiteName = SL_df$SiteName,
                                      ID = "Non-Linear Local Component")
  names(space_time_component_df)<-c("RSL","upr", "lwr","lwr_50","upr_50", "Age",
                                    "Longitude","Latitude","SiteName","ID")
  
  space_time_component_df <- space_time_component_df %>%
    filter(SiteName %in% c(
      "East River Marsh,\n Connecticut",
      "Swan Key,\n Florida",
      "Placentia,\n Newfoundland",
      "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  
  
  # --Filtering for 4 sites---
  SL_df <-SL_df %>%   filter(SiteName %in% c("East River Marsh,\n Connecticut",
                                             "Swan Key,\n Florida",
                                             "Placentia,\n Newfoundland",
                                             "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  
  #------Ordering the sites original Data set-------
  all_data_sites<-factor(paste0(SL_df$Longitude, SL_df$Longitude), labels = 1:4)
  SL_df <- cbind(SL_df,all_data_sites=all_data_sites)
  order_sites <- SL_df %>%  group_by(SiteName,all_data_sites) %>%
    dplyr::summarise(n = n()) %>%
    arrange(all_data_sites)
  SL_df$SiteName <- factor(SL_df$SiteName,levels = unique(order_sites$SiteName))
  space_time_component_df$SiteName <- factor(SL_df$SiteName,levels = unique(order_sites$SiteName))
  
  #---Total plot----
  plot_result <- ggplot() +
    facet_wrap(~SiteName)+
    geom_rect(data = SL_df, aes(
      xmin = Age*1000 - Age_er_average*1000, xmax = Age*1000 + Age_er_average*1000,
      ymin = RSL - RSL_er_average, ymax = RSL + RSL_er_average,fill = "Uncertainty"), alpha = 0.5) +
    geom_point(data = SL_df, aes(y = RSL, x =  Age*1000, colour = "black"), size = 0.6) +
    geom_line(data = mod_output_df, aes(x = Age*1000, y = RSL,colour="mean")) +
    geom_ribbon(data=mod_output_df,
                aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill="CI95"),alpha=0.3)+
    geom_ribbon(data=mod_output_df,
                aes(y = RSL,ymin=lwr_50,ymax=upr_50,x=Age*1000,fill="CI50"),alpha=0.4)+
    xlab("Year (CE)") +
    ylab("Relative Sea Level (m)") +
    theme_bw()+
    theme(axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=10),
          strip.text.x = element_text(size = 14),
          strip.background =element_rect(fill=c("white")))+
    ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
    ggplot2::scale_fill_manual("",
                               values = c(
                                 "Uncertainty" = ggplot2::alpha("grey", 0.6),
                                 "CI95" = ggplot2::alpha("purple3", 0.3),
                                 "CI50" = ggplot2::alpha("purple3", 0.4)
                               ),
                               labels = c(
                                 CI95 = paste0("95% Credible Interval"),
                                 CI50 = paste0("50% Credible Interval"),
                                 Uncertainty = expression(paste("1",sigma," " ,"RSL and age uncertainty"))
                               )
    ) +
    ggplot2::scale_colour_manual("",
                                 values = c(
                                   "black" = "black",
                                   "mean" = "purple3"
                                 ),
                                 # "mean" = "purple3"),
                                 labels = c("Proxy Records", "Posterior Mean Fit")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(
        alpha = c(0.6, 0.3, 0.4),
        size = 0.75
      )),
      colour = ggplot2::guide_legend(override.aes = list(
        linetype = c(0, 1),
        shape = c(16, NA),
        size = 1
      ))
    ) 
  plot_result
  ggsave(paste0("fig/",save_option,"/",
                unique(mod_output_df$ID), ".pdf", sep = ""),
         plot_result, width = 10, height = 6)
  
  #---Regional Plot----
  regional_plot<-ggplot()+
    geom_line(data=time_component_df,aes(x=Age*1000,y=RSL),colour="dodgerblue1")+#"#3b47ad")+
    geom_ribbon(data=time_component_df,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="dodgerblue1",#"#3b47ad",
                alpha=0.3)+
    geom_ribbon(data=time_component_df,aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="dodgerblue1",#"#3b47ad",
                alpha=0.4)+
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
  
  cat("Regional Component plotted")
  
  #------Linear Local Component + site specific vertical offset: Random Effect------
  g_h_plot<-ggplot()+
    geom_line(data=g_h_component_df,aes(x=Age*1000,y=RSL),colour="#5bac06")+
    geom_ribbon(data=g_h_component_df,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#5bac06",alpha=0.3)+
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    theme(strip.text.x = element_text(size = 14),
          strip.background =element_rect(fill=c("white")))+
    ylab('Sea Level (m)')+
    facet_wrap(~SiteName) +
    xlab('Year(CE)')
  g_h_plot
  ggsave(paste0("fig/",save_option,"/",
                unique(g_h_component_df$ID), ".pdf", sep = ""),
         g_h_plot, width = 10, height = 6)
  
  cat("Linear Local Component + site specific vertical offset plotted")
  
  #----Local Plot---  
  local_plot<-ggplot()+
    geom_line(data=space_time_component_df,aes(x=Age*1000,y=RSL,colour="mean"))+
    geom_ribbon(data=space_time_component_df,
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000,fill="CI50"),alpha=0.4)+
    geom_ribbon(data=space_time_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill="CI95"),alpha=0.3)+
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
    xlab('Year (CE)')+
    ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
    ggplot2::scale_fill_manual("",
                               values = c(
                                 "CI95" = ggplot2::alpha("#ad4c14", 0.3),
                                 "CI50" = ggplot2::alpha("#ad4c14", 0.4)
                               ),
                               labels = c(
                                 CI95 = paste0("95% Credible Interval"),
                                 CI50 = paste0("50% Credible Interval")
                               )
    ) +
    ggplot2::scale_colour_manual("",
                                 values = c(
                                   "mean" = "#ad4c14"
                                 ),
                                 labels = c("Posterior Mean Fit")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(
        alpha = c(0.4, 0.3),
        size = 0.75
      )),
      colour = ggplot2::guide_legend(override.aes = list(
        linetype = c( 1),
        shape = c(NA),
        size = 1
      ))
    ) 
  
  local_plot
  #ggsave(local_plot,file = "fig/local_component.pdf", width = 10, height = 6)
  ggsave(paste0("fig/",save_option,"/",
                unique(space_time_component_df$ID), ".pdf", sep = ""),
         local_plot, width = 10, height = 6)
  cat("Non-Linear Local Component plotted")
  
  
  #----Separate Components on one plot with CI----
  all_components_CI_plot <- ggplot()+
    # Local
    geom_line(data=space_time_component_df,
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=space_time_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
    
    # Linear Local Component + y offset
    geom_line(data=g_h_component_df,
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=g_h_component_df,
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
    scale_fill_manual(name = '',values = c("#5bac06","#ad4c14","dodgerblue1",#"#3b47ad",
                                           "purple3"),
                      guide = guide_legend(override.aes = list(alpha = 0.4)))+
    scale_colour_manual(name = '',values = c("#5bac06","#ad4c14","dodgerblue1",#"#3b47ad",
                                             "purple3"))+
    theme_bw()+
    theme(
          axis.title=element_text(size=12,face="bold"),
          axis.text=element_text(size=10),
          legend.text=element_text(size=10),
          strip.text.x = element_text(size = 14),
          strip.background =element_rect(fill=c("white")))+
    theme(legend.position="bottom", legend.box = "horizontal")+
    guides(color = guide_legend(override.aes = list(size = 0.5)))+
    ylab('Sea Level (m)')+
    facet_wrap(~SiteName) +
    xlab(' Year (CE)')
  
  all_components_CI_plot
  ggsave(paste0("fig/",save_option,"/","all_component_plot.pdf", sep = ""),
         all_components_CI_plot, width = 10, height = 6)
  
  cat("All Components plotted")
}