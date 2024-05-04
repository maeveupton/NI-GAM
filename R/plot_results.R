plot_results <- function(SL_df,model_run,basis_fun_list,
                         save_option,set_up_option,n_sites){
  
  dir.create(paste0("fig/",save_option),showWarnings = F)
  
  #-----------Plotting Overall Output----------
  mu_post <- model_run$BUGSoutput$sims.list$mu
  #---Get estimates and uncertainty bounds--
  RSL_mod<-apply(mu_post,2,mean)
  RSL_mod_upr<-apply(mu_post,2,quantile,probs=0.025)
  RSL_mod_lwr<-apply(mu_post,2,quantile,probs=0.975)
  upr_50<-apply(mu_post,2,quantile,probs=0.25)
  lwr_50<-apply(mu_post,2,quantile,probs=0.75)
  
  
  #-----Create data frame for plotting-----
  mod_output_df<-data.frame(RSL_mod,RSL_mod_upr,RSL_mod_lwr,
                            upr_50,lwr_50,
                            SL_df$Age,
                            SL_df$SiteName,
                            ID = "Total Posterior Model",
                            data_type_id = SL_df$data_type_id)
  names(mod_output_df)<-c("RSL","upr","lwr","upr_50","lwr_50",
                          "Age",
                          "SiteName",
                          "ID","data_type_id")
  
  #---Total plot----
  plot_result <- ggplot() +
    facet_wrap(~SiteName)+
    geom_rect(data = subset(SL_df,data_type_id =="ProxyData"), aes(
      xmin = Age*1000 - Age_er_average*1000, xmax = Age*1000 + Age_er_average*1000,
      ymin = RSL - RSL_er_average, ymax = RSL + RSL_er_average,fill = "Uncertainty"), alpha = 0.9) +
    geom_point(data = subset(SL_df,data_type_id =="ProxyData"),
               aes(y = RSL, x =  Age*1000, colour = "black"), size = 0.3) +
    geom_line(data = subset(mod_output_df,data_type_id =="ProxyData"), 
              aes(x = Age*1000, y = RSL,colour="mean")) +
    geom_ribbon(data=subset(mod_output_df,data_type_id =="ProxyData"),
                aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill="CI95"),alpha=0.3)+
    geom_ribbon(data=subset(mod_output_df,data_type_id =="ProxyData"),
                aes(y = RSL,ymin=lwr_50,ymax=upr_50,x=Age*1000,fill="CI50"),alpha=0.4)+
    xlab("Year (CE)") +
    ylab("Relative Sea Level (m)") +
    theme_bw()+
    theme(axis.title=element_text(size=10,face="bold"),
          axis.text=element_text(size=8),
          legend.text=element_text(size=9),
          strip.text.x = element_text(size = 8),
          strip.background =element_rect(fill=c("white")))+
    #ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
    ggplot2::scale_fill_manual("",
                               values = c(
                                 "Uncertainty" = ggplot2::alpha("grey", 0.9),
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
                                 labels = c("Proxy Records", "Posterior Mean Fit")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(
        alpha = c(0.9, 0.3, 0.4),
        size = 0.75
      )),
      colour = ggplot2::guide_legend(override.aes = list(
        linetype = c(0, 1),
        shape = c(16, NA),
        size = 1
      ))
    ) +
    theme(
      #legend.position=c(0.90,-0.05),
      legend.position=c(1,0),
          legend.direction = "horizontal",
          legend.justification = c(1,0),
          legend.spacing.y = unit(0.001, 'cm'),
          legend.title=element_blank(),
          legend.margin=margin(c(1,1,1,1)))
  plot_result
  ggsave(paste0("fig/",save_option,"/",
                unique(mod_output_df$ID), "proxy.pdf", sep = ""),
         plot_result, width = 10, height = 6)

  #---Total plot for TG----
  SL_df_TG <- SL_df %>%
    filter(data_type_id == "TideGaugeData") %>% 
    group_by(SiteName) %>% 
    mutate(Age_range = max(Age*1000) - min(Age*1000))
  SL_df_TG_names <- SL_df_TG %>% 
    filter(Age_range>=50) %>% 
    dplyr::select(SiteName) %>% unique()
  mod_output_df_TG <- mod_output_df %>% 
    group_by(SiteName) %>% 
    filter(data_type_id == "TideGaugeData") %>% 
    mutate(Age_range = max(Age*1000) - min(Age*1000)) #%>% 
    #filter(SiteName %in% c(SL_df_TG_names))
  
  plot_result_TG <- ggplot() +
    facet_wrap(~SiteName)+
    geom_rect(data = subset(SL_df_TG,Age_range >= 50), aes(
      xmin = Age*1000 - Age_er_average*1000, xmax = Age*1000 + Age_er_average*1000,
      ymin = RSL - RSL_er_average, ymax = RSL + RSL_er_average,fill = "Uncertainty"), alpha = 0.9) +
    geom_point(data = subset(SL_df_TG,Age_range >= 50),
               aes(y = RSL, x =  Age*1000, colour = "black"), size = 0.3) +
    geom_line(data =  subset(mod_output_df_TG,Age_range >= 50), 
              aes(x = Age*1000, y = RSL,colour="mean")) +
    geom_ribbon(data= subset(mod_output_df_TG,Age_range >= 50),
                aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill="CI95"),alpha=0.3)+
    geom_ribbon(data= subset(mod_output_df_TG,Age_range >= 50),
                aes(y = RSL,ymin=lwr_50,ymax=upr_50,x=Age*1000,fill="CI50"),alpha=0.4)+
    xlab("Year (CE)") +
    ylab("Relative Sea Level (m)") +
    theme_bw()+
    theme(axis.title=element_text(size=10,face="bold"),
          axis.text=element_text(size=8),
          legend.text=element_text(size=9),
          strip.text.x = element_text(size = 8),
          strip.background =element_rect(fill=c("white")))+
    #ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
    ggplot2::scale_fill_manual("",
                               values = c(
                                 "Uncertainty" = ggplot2::alpha("grey", 0.9),
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
                                 labels = c("Tide Gauge Data", "Posterior Mean Fit")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(
        alpha = c(0.9, 0.3, 0.4),
        size = 0.75
      )),
      colour = ggplot2::guide_legend(override.aes = list(
        linetype = c(0, 1),
        shape = c(16, NA),
        size = 1
      ))
    ) +
    theme(
      #legend.position=c(0.90,-0.05),
      legend.position=c(1,-0.06),
      #legend.direction = "horizontal",
      legend.justification = c(1,0),
      legend.spacing.y = unit(0.001, 'cm'),
      legend.title=element_blank(),
      legend.margin=margin(c(1,1,1,1)))
  plot_result_TG
  ggsave(paste0("fig/",save_option,"/",
                unique(mod_output_df_TG$ID), "TG.pdf", sep = ""),
         plot_result_TG, width = 10, height = 6)
  
    
  #------------Regional Component: Spline in Time--------------
  time_component <- model_run$BUGSoutput$sims.list$r
  RSL_mod_NI<-apply(time_component,2,mean)
  RSL_mod_lwr_NI<-apply(time_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(time_component,2,quantile,probs=0.975)
  lwr_50<-apply(time_component,2,quantile,probs=0.25)
  upr_50<-apply(time_component,2,quantile,probs=0.75)
  time_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,lwr_50,upr_50,
                                Age = SL_df$Age,
                                SiteName = SL_df$SiteName,
                                ID = "Regional Component",
                                data_type_id = SL_df$data_type_id)
  names(time_component_df)<-c("RSL","upr","lwr","lwr_50","upr_50",
                              "Age","SiteName","ID","data_type_id")
  
  
  # regional_plot<-ggplot()+
  #   geom_line(data=time_component_df,aes(x=Age*1000,y=RSL),colour="dodgerblue1")+#"#3b47ad")+
  #   geom_ribbon(data=time_component_df,
  #               aes(ymin=lwr,ymax=upr,x=Age*1000),fill="dodgerblue1",#"#3b47ad",
  #               alpha=0.2)+
  #   geom_ribbon(data=time_component_df,
  #               aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="dodgerblue1",#"#3b47ad",
  #               alpha=0.3)+
  #   theme_bw()+
  #   theme(plot.title = element_text(size=15),
  #         axis.title=element_text(size=12,face="bold"),
  #         axis.text=element_text(size=8),
  #         legend.text=element_text(size=12))+
  #   scale_x_continuous(breaks = round(seq(min(basis_fun_list$Age_grid*1000)-7,
  #                                         max(basis_fun_list$Age_grid*1000)+7, by = 100),1))+
  #   ylab('Sea Level (m)')+
  #   xlab('Year (CE)')
  # regional_plot
  # ggsave(paste0("fig/",save_option,'/' ,
  #               unique(time_component_df$ID), ".pdf", sep = ""),
  #        regional_plot, width = 10, height = 6)
  # 
  # cat("Regional Component plotted")

  
  #------------Regional Grid Component--------------
  time_grid_component <- model_run$BUGSoutput$sims.list$r_grid
  RSL_mod_NI<-apply(time_grid_component,2,mean)
  RSL_mod_lwr_NI<-apply(time_grid_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(time_grid_component,2,quantile,probs=0.975)
  lwr_50<-apply(time_grid_component,2,quantile,probs=0.25)
  upr_50<-apply(time_grid_component,2,quantile,probs=0.75)
  time_grid_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
                                     lwr_50,upr_50, 
                                     Age = basis_fun_list$Age_grid,
                                     ID = "Regional Grid Component")
  names(time_grid_component_df)<-c("RSL","upr","lwr","lwr_50","upr_50",
                                   "Age","ID")
  
  regional_grid_plot<-ggplot()+
    geom_line(data=time_grid_component_df,aes(x=Age*1000,y=RSL),colour="dodgerblue1")+#"#3b47ad")+
    geom_ribbon(data=time_grid_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000),fill="dodgerblue1",#"#3b47ad",
                alpha=0.2)+
    geom_ribbon(data=time_grid_component_df,
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="dodgerblue1",#"#3b47ad",
                alpha=0.3)+
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=8),
          legend.text=element_text(size=12))+
    scale_x_continuous(breaks = round(seq(min(basis_fun_list$Age_grid*1000)-7,
                                          max(basis_fun_list$Age_grid*1000)+7, by = 100),1))+
    ylab('Sea Level (m)')+
    xlab('Year (CE)')
  regional_grid_plot
  ggsave(paste0("fig/",save_option,'/' ,
                unique(time_grid_component_df$ID), ".pdf", sep = ""),
         regional_grid_plot, width = 10, height = 6)
  cat("Regional Grid Component plotted")
  
  
  #----- Plotting all posterior for regional component-----
  time_component_grid <- model_run$BUGSoutput$sims.list$r_grid
  #post_regional <- as.data.frame(time_component_grid[1:5,])
  #post_regional <- as.data.frame(time_component_grid[sample(nrow(time_component_grid),5),])
  post_regional <- as.data.frame(time_component_grid[sample(nrow(time_component_grid),10),])
  rownames(post_regional) <- paste0("sample_",1:nrow(post_regional))
  post_regional_long_df <- pivot_longer(data = as.data.frame(cbind(Age = basis_fun_list$Age_grid,t(post_regional))),
                                        cols = starts_with("sample_"))
  post_regional_long_df <- post_regional_long_df %>% mutate(post_run = as.factor(name))
  post_regional_long_df <- post_regional_long_df %>% dplyr::select(-name)%>% 
    group_by(post_run) %>% arrange(desc(Age))
  regional_post_plot <- ggplot(data  = post_regional_long_df,
                               aes(Age,value,colour = post_run))+
    #aes(Age,value,colour = name))+
    geom_line()
  ggsave(paste0("fig/",save_option,'/',"regional_post_plot.pdf", sep = ""),
         regional_post_plot, width = 10, height = 6)
  
  
  # regional_grid_plot<-ggplot()+
  #   geom_line(data=time_grid_component_df,aes(x=Age*1000,y=RSL),colour="#3b47ad")+
  #   geom_ribbon(data=time_grid_component_df,
  #               aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#3b47ad",alpha=0.2)+
  #   geom_ribbon(data=time_grid_component_df,
  #               aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="#3b47ad",alpha=0.3)+
  #   theme_bw()+
  #   theme(legend.position="none")+
  #   theme(plot.title = element_text(size=22),
  #         axis.title=element_text(size=14,face="bold"),
  #         axis.text=element_text(size=12),
  #         legend.text=element_text(size=12))+
  #   ylab('Sea Level (m)')+
  #   xlab('Year (CE)')
  # regional_grid_plot
  # ggsave(paste0("fig/",save_option,'/' ,
  #               unique(time_grid_component_df$ID), ".pdf", sep = ""),
  #        regional_grid_plot, width = 10, height = 6)
  
  
  regional_grid_plot_post<-ggplot()+
    geom_line(data=time_grid_component_df,aes(x=Age*1000,y=RSL),colour="dodgerblue1")+
    geom_ribbon(data=time_grid_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000),fill="dodgerblue1",alpha=0.2)+
    geom_ribbon(data=time_grid_component_df,
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="dodgerblue1",alpha=0.3)+
    geom_line(data  = post_regional_long_df,aes(Age*1000,value,
                                                colour=post_run),alpha = 0.5)+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=12,face="bold"),
          axis.text=element_text(size=8),
          legend.text=element_text(size=12))+
    scale_colour_manual(values=c("gray48", "gray48","gray48","gray48","gray48",
                                 "gray48","gray48","gray48","gray48","gray48"))+
    annotate(geom="text", x=1500, y=-0.19, label="10 Posterior Samples",
             color="gray48")+
    ylab('Sea Level (m)')+
    scale_x_continuous(breaks = round(seq(min(basis_fun_list$Age_grid*1000)-7,
                                          max(basis_fun_list$Age_grid*1000)+7, by = 100),1))+
    xlab('Year (CE)')
  regional_grid_plot_post
  ggsave(paste0("fig/",save_option,'/' ,
                unique(time_grid_component_df$ID), "post.pdf", sep = ""),
         regional_grid_plot_post, width = 10, height = 6)
  
  #------Linear Local Component: Random Effect------
  lin_loc_component <- model_run$BUGSoutput$sims.list$g_z_x
  RSL_mod_NI<-apply(lin_loc_component,2,mean)
  RSL_mod_lwr_NI<-apply(lin_loc_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(lin_loc_component,2,quantile,probs=0.975)
  lin_loc_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
                                   Age = SL_df$Age,
                                   SiteName = SL_df$SiteName,
                                   data_type_id = SL_df$data_type_id,
                                   ID = "Linear Local Component")
  names(lin_loc_component_df)<-c("RSL","upr","lwr","Age","SiteName","ID","data_type_id")
  
  #---Linear Local + Site specific vertical offset---
  g_h_component <- model_run$BUGSoutput$sims.list$g_h_z_x
  RSL_mod_NI<-apply(g_h_component,2,mean)
  RSL_mod_lwr_NI<-apply(g_h_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(g_h_component,2,quantile,probs=0.975)
  g_h_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
                               Age = SL_df$Age,
                               SiteName = SL_df$SiteName,
                               ID = "Linear Local Component and site-specific vertical offset",
                               data_type_id = SL_df$data_type_id)
  names(g_h_component_df)<-c("RSL","upr","lwr","Age","SiteName","ID", "data_type_id")
  
  g_h_component_df_plot<-ggplot()+
    geom_line(data=subset(g_h_component_df,data_type_id =="ProxyData"),
              aes(x=Age*1000,y=RSL),colour="#5bac06")+
    geom_ribbon(data=subset(g_h_component_df,data_type_id =="ProxyData"),
                aes(ymin=lwr,ymax=upr,x=Age*1000),
                fill="#5bac06",alpha=0.3)+
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    theme(strip.text.x = element_text(size = 14),
          strip.background =element_rect(fill=c("white")))+
    ylab('Sea Level (m)')+
    facet_wrap(~SiteName) +
    xlab('Year (CE)')
  g_h_component_df_plot
  ggsave(paste0("fig/",save_option,"/",
                unique(g_h_component_df_plot$ID), ".pdf", sep = ""),
         g_h_component_df_plot, width = 10, height = 6)
  
  #------Non-Linear Local Component: Spline in Space Time------
  space_time_component <- model_run$BUGSoutput$sims.list$l
  RSL_mod_NI<-apply(space_time_component,2,mean)
  RSL_mod_lwr_NI<-apply(space_time_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(space_time_component,2,quantile,probs=0.975)
  lwr_50<-apply(space_time_component,2,quantile,probs=0.25)
  upr_50<-apply(space_time_component,2,quantile,probs=0.75)
  space_time_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,lwr_50,upr_50,
                                      Age = SL_df$Age,
                                      SiteName = SL_df$SiteName,
                                      ID = "Non-Linear Local Component",
                                      data_type_id = SL_df$data_type_id)
  names(space_time_component_df)<-c("RSL","upr","lwr","lwr_50","upr_50",
                                    "Age","SiteName","ID","data_type_id")
  
  local_plot1<-ggplot()+
    geom_line(data=subset(space_time_component_df,data_type_id =="ProxyData"),
              aes(x=Age*1000,y=RSL,colour="mean"))+
    geom_ribbon(data=subset(space_time_component_df,data_type_id =="ProxyData"),
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000,fill="CI50"),alpha=0.4)+
    geom_ribbon(data=subset(space_time_component_df,data_type_id =="ProxyData"),
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill="CI95"),alpha=0.3)+
    geom_hline(yintercept = 0)+
    theme_bw()+
    ylab('Sea Level (m)')+
    theme(
          axis.title=element_text(size=10,face="bold"),
          axis.text=element_text(size=9),
          legend.text=element_text(size=10))+
    facet_wrap(~SiteName) +
    theme(strip.text.x = element_text(size = 8),
          strip.background =element_rect(fill=c("white")))+
    xlab('Year (CE)')+
    #ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
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
        alpha = c(0.5, 0.3),
        size = 0.75
      )),
      colour = ggplot2::guide_legend(override.aes = list(
        linetype = c( 1),
        shape = c(NA),
        size = 1
      ))
    ) +
    theme(legend.position=c(0.95,-0.05),
          legend.justification = c(1,0),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title=element_blank(),
          legend.margin=margin(c(1,1,1,1)))
  
  local_plot1
  
  ggsave(paste0("fig/",save_option,"/",
                unique(space_time_component_df$ID), "proxy.pdf", sep = ""),
         local_plot1, width = 10, height = 6)
  
  #----Separate Components on one plot with CI proxy----
  all_components_CI_plot <- ggplot()+
    # Local
    geom_line(data=subset(space_time_component_df,data_type_id =="ProxyData"),
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=subset(space_time_component_df,data_type_id =="ProxyData"),
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
    # Linear Local Component + site specific vertical offset
    geom_line(data=subset(g_h_component_df,data_type_id == "ProxyData"),
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=subset(g_h_component_df,data_type_id == "ProxyData"),
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
    
    # Regional Component
    geom_line(data=subset(time_component_df,data_type_id == "ProxyData"),
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=subset(time_component_df,data_type_id == "ProxyData"),
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
    # Total Model
    geom_line(data=subset(mod_output_df,data_type_id == "ProxyData"),
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=subset(mod_output_df,data_type_id == "ProxyData"),
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
    
    theme_bw()+
    theme(strip.text.x = element_text(size = 9),
          strip.background =element_rect(fill=c("white")))+
    scale_fill_manual(name = '',values = c("#5bac06","#ad4c14","dodgerblue1",#3b47ad",
                                           "purple3"),
                      guide = guide_legend(override.aes = list(alpha = 0.4)))+
    scale_colour_manual(name = '',values = c("#5bac06","#ad4c14","dodgerblue1",#"#3b47ad",
                                             "purple3"))+
    
    ylab('Sea Level (m)')+
    facet_wrap(~SiteName) +
    theme(
      axis.title=element_text(size=10,face="bold"),
      axis.text=element_text(size=9))+
    xlab('Year (CE)')+
    theme(legend.position=c(0.95,-0.05),
          legend.justification = c(1,0),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title=element_blank(),
          legend.margin=margin(c(1,1,1,1)))
    #theme(legend.position="bottom", legend.box = "horizontal")+
    guides(color = guide_legend(override.aes = list(size = 0.5)))
  all_components_CI_plot
  ggsave(paste0("fig/",save_option,"/","all_component_plot.pdf", sep = ""),
         all_components_CI_plot, width = 10, height = 6)
  
  cat("All Components plotted")
}


