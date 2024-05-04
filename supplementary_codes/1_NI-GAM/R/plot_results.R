plot_results <- function(SL_df,model_run,basis_fun_list,
                         save_option,set_up_option,n_sites){
  
  dir.create(paste0("fig/",save_option),showWarnings = F)
  
  #-----------Plotting Overall Output----------
  #mu_post <- model_run$BUGSoutput$sims.list$mu_overall
  mu_post <- model_run$BUGSoutput$sims.list$mu
  #---Get estimates and uncertainty bounds--
  RSL_mod<-apply(mu_post,2,mean)
  RSL_mod_upr<-apply(mu_post,2,quantile,probs=0.025)
  RSL_mod_lwr<-apply(mu_post,2,quantile,probs=0.975)
  upr_50<-apply(mu_post,2,quantile,probs=0.25)
  lwr_50<-apply(mu_post,2,quantile,probs=0.75)
  
  # #------Ordering the sites-------
  # all_data_sites<-factor(paste0(SL_df$Longitude, SL_df$Longitude), labels = 1:n_sites)#For ordering sites
  # SL_df <- cbind(SL_df,all_data_sites=all_data_sites)
  # order_sites <- SL_df %>%  group_by(SiteName,all_data_sites) %>%
  #   dplyr::summarise(n = n()) %>%
  #   arrange(all_data_sites)
  # SL_df$SiteName <- factor(SL_df$SiteName,levels = unique(order_sites$SiteName))
  
  
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
  write_csv(mod_output_df,"data/result_csv/total_model_fit.csv")
  plot_result <- ggplot() +
    geom_point(data = subset(SL_df,data_type_id =="ProxyData"), aes(y = RSL, x = Age*1000), size = 0.5) +
     geom_rect(data = subset(SL_df,data_type_id =="ProxyData"), aes(
      xmin = Age*1000 - Age_er_average*1000, xmax = Age*1000 + Age_er_average*1000,
      ymin = RSL - RSL_er_average, ymax = RSL + RSL_er_average
     ), alpha = 0.2) +
    geom_line(data = subset(mod_output_df,data_type_id =="ProxyData"), aes(x = Age*1000, y = RSL), colour = "purple3") +
    geom_ribbon(data=subset(mod_output_df,data_type_id =="ProxyData"),
                aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000),fill = "purple3",alpha=0.2)+
    geom_ribbon(data=subset(mod_output_df,data_type_id =="ProxyData"),
                aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000),fill = "purple3",alpha=0.3)+
    
    facet_wrap(~SiteName)+
       xlab("Year (CE)") +
    ylab("RSL (m)") +
    theme_bw()+
    theme(plot.title = element_text(size=15),
          axis.title=element_text(size=12,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=10))+
    theme(strip.text.x = element_text(size = 10),
          strip.background =element_rect(fill=c("white")))+
    theme(legend.position=c(0.95,0.01),
          legend.justification = c(1,0),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title=element_blank(),
          legend.margin=margin(c(1,1,1,1)))
  plot_result
  ggsave(paste0("fig/",save_option,"/",
                unique(mod_output_df$ID), "proxy.pdf", sep = ""),
         plot_result, width = 10, height = 6)

  # plot_result_2 <- ggplot() +
  #   geom_point(data = subset(SL_df,data_type_id =="TideGaugeData"), aes(y = RSL, x = Age*1000), size = 0.5) +
  #   geom_rect(data = subset(SL_df,data_type_id =="TideGaugeData"), aes(
  #     xmin = Age*1000 - Age_er_average*1000, xmax = Age*1000 + Age_er_average*1000,
  #     ymin = RSL - RSL_er_average, ymax = RSL + RSL_er_average
  #   ), alpha = 0.2) +
  #   geom_line(data = subset(mod_output_df,data_type_id =="TideGaugeData"), aes(x = Age*1000, y = RSL), colour = "purple3") +
  #   geom_ribbon(data=subset(mod_output_df,data_type_id =="TideGaugeData"),
  #               aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000),fill = "purple3",alpha=0.3)+
  #   xlab("Year (CE)") +
  #   ylab("RSL (m)") +
  #   theme_bw()+
  #   theme(plot.title = element_text(size=15),
  #         axis.title=element_text(size=12,face="bold"),
  #         #axis.text=element_text(size=12),
  #         axis.text.x = element_text(size = 3),
  #         axis.text.y = element_text(size = 3),
  #         legend.text=element_text(size=10))+
  #         #axis.text.y=element_text(margin=margin(r=0)),
  #         #axis.text.x=element_text(margin=margin(r=0)))+
  #   theme(strip.text.x = element_text(size = 4),
  #         strip.background =element_rect(fill=c("white")),
  #         panel.spacing = unit(0,'lines'))+
  #   theme(legend.position=c(0.95,0.01),
  #         legend.justification = c(1,0),
  #         legend.spacing.y = unit(0.1, 'cm'),
  #         legend.title=element_blank(),
  #         legend.margin=margin(c(1,1,1,1)))+
  #   facet_wrap(~SiteName,scales = "free",labeller = labeller(groupwrap = label_wrap_gen(10)))
  # plot_result_2
  # ggsave(paste0("fig/",save_option,"/",
  #               unique(mod_output_df$ID), "TG.pdf", sep = ""),
  #        plot_result_2, width = 12, height = 7)
  
  
   #------------Regional Component: Spline in Time--------------
  time_component <- model_run$BUGSoutput$sims.list$regional
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
  
  write_csv(time_component_df,"data/result_csv/regional_model_fit.csv")
  
  regional_plot<-ggplot()+
    geom_line(data=time_component_df,aes(x=Age*1000,y=RSL),colour="#3b47ad")+
    geom_ribbon(data=time_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#3b47ad",alpha=0.2)+
    geom_ribbon(data=time_component_df,
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="#3b47ad",alpha=0.3)+
    #ggtitle("Regional Component")+
    #ggtitle(paste0("Regional Component",set_up_option, sep = " "))+
    theme_bw()+
    theme(plot.title = element_text(size=15),
          axis.title=element_text(size=12,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    ylab('Sea Level (m)')+
    xlab('Age(CE)')
  regional_plot
  ggsave(paste0("fig/",save_option,'/' ,
                unique(time_component_df$ID), ".pdf", sep = ""),
         regional_plot, width = 10, height = 6)
  #ggsave(regional_plot,file = "fig/regional_component.pdf", width = 10, height = 6)
  
  cat("Regional Component plotted")
  
  #------------Regional Grid Component--------------
  # Right for full data set but not for 4 sites
  time_grid_component <- model_run$BUGSoutput$sims.list$regional_grid
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
  
  write_csv(time_grid_component_df,"data/result_csv/regional_grid_model_fit.csv")
  
  #----- Plotting all posterior for regional component-----
  time_component_grid <- model_run$BUGSoutput$sims.list$regional_grid
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

  
  regional_grid_plot<-ggplot()+
    geom_line(data=time_grid_component_df,aes(x=Age*1000,y=RSL),colour="#3b47ad")+
    geom_ribbon(data=time_grid_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#3b47ad",alpha=0.2)+
    geom_ribbon(data=time_grid_component_df,
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="#3b47ad",alpha=0.3)+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    ylab('Sea Level (m)')+
    xlab('Year (CE)')
  regional_grid_plot
  ggsave(paste0("fig/",save_option,'/' ,
                unique(time_grid_component_df$ID), ".pdf", sep = ""),
         regional_grid_plot, width = 10, height = 6)


  regional_grid_plot_post<-ggplot()+
    geom_line(data=time_grid_component_df,aes(x=Age*1000,y=RSL),colour="#3b47ad")+
    geom_ribbon(data=time_grid_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#3b47ad",alpha=0.2)+
    geom_ribbon(data=time_grid_component_df,
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="#3b47ad",alpha=0.3)+
    geom_line(data  = post_regional_long_df,aes(Age*1000,value,
                                                colour=post_run),alpha = 0.5)+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    scale_colour_manual(values=c("#999999", "#999999","#999999","#999999","#999999",
                                 "#999999","#999999","#999999","#999999","#999999"))+
    annotate(geom="text", x=1500, y=-0.19, label="10 Posterior Samples",
             color="#999999")+
    ylab('Sea Level (m)')+
    xlab('Year (CE)')
  regional_grid_plot_post
  ggsave(paste0("fig/",save_option,'/' ,
                unique(time_grid_component_df$ID), "post.pdf", sep = ""),
         regional_grid_plot_post, width = 10, height = 6)
  
  
  cat("Regional Grid Component plotted")
  
  #------Linear Local Component: Random Effect------
  GIA_component <- model_run$BUGSoutput$sims.list$GIA_correction
  RSL_mod_NI<-apply(GIA_component,2,mean)
  RSL_mod_lwr_NI<-apply(GIA_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(GIA_component,2,quantile,probs=0.975)
  GIA_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
                               Age = SL_df$Age,
                               SiteName = SL_df$SiteName,
                               data_type_id = SL_df$data_type_id,
                               ID = "Linear Local Component")
  names(GIA_component_df)<-c("RSL","upr","lwr","Age","SiteName","ID","data_type_id")
  
  write_csv(GIA_component_df,"data/result_csv/GIA_model_fit.csv")
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
  #   xlab('Age(CE)')
  # GIA_plot
  # ggsave(paste0("fig/",save_option,"/",
  #               unique(GIA_component_df$ID), ".pdf", sep = ""),
  #        GIA_plot, width = 10, height = 6)
  # #ggsave(GIA_plot,file = "fig/GIA_component.pdf", width = 10, height = 6)
  
  cat("Linear Local Component plotted")
  
  #---GIA + yoffset---
  GIA_yoffset_component <- model_run$BUGSoutput$sims.list$GIA_yoffset
  RSL_mod_NI<-apply(GIA_yoffset_component,2,mean)
  RSL_mod_lwr_NI<-apply(GIA_yoffset_component,2,quantile,probs=0.025)
  RSL_mod_upr_NI<-apply(GIA_yoffset_component,2,quantile,probs=0.975)
  GIA_yoffset_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
                                       Age = SL_df$Age,
                                       SiteName = SL_df$SiteName,
                                       ID = "Linear Local Component and site-specific vertical offset",
                                       data_type_id = SL_df$data_type_id)
  names(GIA_yoffset_component_df)<-c("RSL","upr","lwr","Age","SiteName","ID", "data_type_id")
  write_csv(GIA_yoffset_component_df,"data/result_csv/gia_yoffset_model_fit.csv")
  
  #------Non-Linear Local Component: Spline in Space Time------
  space_time_component <- model_run$BUGSoutput$sims.list$local
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
  
  write_csv(space_time_component_df,"data/result_csv/local_model_fit.csv")
  
  local_plot1<-ggplot()+
    geom_line(data=subset(space_time_component_df,data_type_id =="ProxyData"),
              aes(x=Age*1000,y=RSL),colour="#ad4c14")+
    geom_ribbon(data=subset(space_time_component_df,data_type_id =="ProxyData"),
                aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#ad4c14",alpha=0.2)+
    geom_ribbon(data=subset(space_time_component_df,data_type_id =="ProxyData"),
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),fill="#ad4c14",alpha=0.3)+
    geom_hline(yintercept = 0)+
    theme_bw()+
    ylab('Sea Level (m)')+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    facet_wrap(~SiteName) +
    theme(strip.text.x = element_text(size = 10),
          strip.background =element_rect(fill=c("white")))+
    xlab('Age(CE)')
  local_plot1
  ggsave(paste0("fig/",save_option,"/",
                unique(space_time_component_df$ID), "proxy.pdf", sep = ""),
         local_plot1, width = 10, height = 6)
  
  # local_plot<-ggplot()+
  #   geom_line(data=subset(space_time_component_df,data_type_id =="TideGaugeData"),aes(x=Age*1000,y=RSL),colour="#ad4c14")+
  #   geom_ribbon(data=subset(space_time_component_df,data_type_id =="TideGaugeData"),aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#ad4c14",alpha=0.3)+
  #   geom_hline(yintercept = 0)+
  #   theme_bw()+
  #   ylab('Sea Level (m)')+
  #   theme(plot.title = element_text(size=22),
  #         axis.title=element_text(size=14,face="bold"),
  #         axis.text=element_text(size=12),
  #         legend.text=element_text(size=12))+
  #   facet_wrap(~SiteName,scales = "free") +
  #   theme(strip.text.x = element_text(size = 14),
  #         strip.background =element_rect(fill=c("white")))+
  #   xlab('Age(CE)')
  # local_plot
  # #ggsave(local_plot,file = "fig/local_component.pdf", width = 10, height = 6)
  # ggsave(paste0("fig/",save_option,"/",
  #               unique(space_time_component_df$ID), "TG.pdf", sep = ""),
  #        local_plot, width = 10, height = 6)
  # cat("Non-Linear Local Component plotted")
  
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
  
  #----Separate Components on one plot with CI proxy----
  all_components_CI_plot <- ggplot()+
      # Local
    geom_line(data=subset(space_time_component_df,data_type_id =="ProxyData"),
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=subset(space_time_component_df,data_type_id =="ProxyData"),
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
    # Linear Local Component + y offset
    geom_line(data=subset(GIA_yoffset_component_df,data_type_id == "ProxyData"),
              aes(x=Age*1000,y=RSL,colour = ID))+
    geom_ribbon(data=subset(GIA_yoffset_component_df,data_type_id == "ProxyData"),
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
    
    # regional Component
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
    #ggtitle(paste0("All Components on 1 plot with",set_up_option, sep = " "))+
    #ggtitle(paste0("Decomposition of the RSL field",set_up_option, sep = " "))+
    theme(strip.text.x = element_text(size = 9),
          strip.background =element_rect(fill=c("white")))+
    scale_fill_manual(name = '',values = c("#5bac06","#ad4c14","#3b47ad","purple3"),
                      guide = guide_legend(override.aes = list(alpha = 0.1)))+
    #guides(fill = guide_legend(override.aes = list(shape = c(19,19,19,19,NA,NA),
    #                                               linetype = c(0,0,0,0, 0,1))))+
    scale_colour_manual(name = '',values = c("#5bac06","#ad4c14","#3b47ad","purple3"))+
    
    ylab('RSL (m)')+
    facet_wrap(~SiteName) +
    xlab('Year (CE)')+
    theme(legend.position=c(0.95,-0.05),
          legend.justification = c(1,0),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title=element_blank(),
          legend.margin=margin(c(1,1,1,1)))
  all_components_CI_plot
  #ggsave(all_components_CI_plot,file = "fig/all_components_one_plot.pdf", width = 10, height = 6)
  ggsave(paste0("fig/",save_option,"/","all_component_plot.pdf", sep = ""),
         all_components_CI_plot, width = 10, height = 6)
  
  # #----Separate Components on one plot with CI for TG----
  # all_components_CI_plotTG <- ggplot()+
  #   # Local
  #   geom_line(data=subset(space_time_component_df,data_type_id =="TideGaugeData"),
  #             aes(x=Age*1000,y=RSL,colour = ID))+
  #   geom_ribbon(data=subset(space_time_component_df,data_type_id =="TideGaugeData"),
  #               aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
  #   # Linear Local Component + y offset
  #   geom_line(data=subset(GIA_yoffset_component_df,data_type_id =="TideGaugeData"),
  #             aes(x=Age*1000,y=RSL,colour = ID))+
  #   geom_ribbon(data=subset(GIA_yoffset_component_df,data_type_id =="TideGaugeData"),
  #               aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
  # 
  #   # Regional Component
  #   geom_line(data=subset(time_component_df,data_type_id == "TideGaugeData"),
  #              aes(x=Age*1000,y=RSL,colour = ID))+
  #   geom_ribbon(data=subset(time_component_df,data_type_id == "TideGaugeData"),
  #                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
  #   # Total Model
  #   geom_line(data=subset(mod_output_df,data_type_id =="TideGaugeData"),
  #             aes(x=Age*1000,y=RSL,colour = ID))+
  #   geom_ribbon(data=subset(mod_output_df,data_type_id =="TideGaugeData"),
  #               aes(ymin=lwr,ymax=upr,x=Age*1000,fill = ID),alpha=0.3)+
  #   theme_bw()+
  #   theme(strip.text.x = element_text(size = 4),
  #         strip.background =element_rect(fill=c("white")))+
  #   scale_fill_manual(name = '',values = c("#5bac06","#ad4c14","#3b47ad","purple3"),
  #                     guide = guide_legend(override.aes = list(alpha = 0.1)))+
  #   #guides(fill = guide_legend(override.aes = list(shape = c(19,19,19,19,NA,NA),
  #   #                                               linetype = c(0,0,0,0, 0,1))))+
  #   scale_colour_manual(name = '',values = c("#5bac06","#ad4c14","#3b47ad","purple3"))+
  #   
  #   ylab('RSL (m)')+
  #   facet_wrap(~SiteName,scales = "free") +
  #   xlab('Year (CE)')+
  #   theme(legend.position=c(0.95,0.01),
  #         legend.justification = c(1,0),
  #         legend.spacing.y = unit(0.1, 'cm'),
  #         legend.title=element_blank(),
  #         legend.margin=margin(c(1,1,1,1)))
  # all_components_CI_plotTG
  # ggsave(paste0("fig/",save_option,"/","all_component_plot_TG.pdf", sep = ""),
  #        all_components_CI_plotTG, width = 12, height = 7)
  
  cat("All Components plotted")
}


