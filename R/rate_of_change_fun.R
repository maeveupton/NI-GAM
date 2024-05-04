#-----------Calculating Rates of Change-------
rate_of_change_fun<- function(SL_df, model_run,save_loc,save_csv,save_option,basis_fun_list){
  #-----Get posterior samples for SL-----
  h_intercept_post <- model_run$BUGSoutput$sims.list$intercept
  b_l_post <- model_run$BUGSoutput$sims.list$b_l
  b_r_post <- model_run$BUGSoutput$sims.list$b_r
  b_g_post <- model_run$BUGSoutput$sims.list$b_g
  n_iterations <- nrow(b_l_post)
  
  # Derivative of Regional Component using a grid
  Age_grid <- seq(min(SL_df$Age),max(SL_df$Age), by = 0.005)
  
  h <- 0.001 # A small number for doing the derivative
  B_regional_grid <- bs_bbase(Age_grid,
                              xl=min(SL_df$Age),xr=max(SL_df$Age),nseg = 20)
  B_new_high_reg <-bs_bbase(Age_grid + h, xl = min(SL_df$Age),
                            xr = max(SL_df$Age), nseg = 20)
  B_new_low_reg <-bs_bbase(Age_grid - h, xl = min(SL_df$Age),
                           xr = max(SL_df$Age), nseg = 20)
  B_deriv_reg <- (B_new_high_reg - B_new_low_reg) / (2 * h)
  
  
  regional_pred_mat <- matrix(NA,nrow=length(Age_grid),ncol = n_iterations)
  for (i in 1:n_iterations){
    regional_pred_mat[,i] <- B_deriv_reg %*% b_r_post[i,]
  }
  time_component_mean <- apply(regional_pred_mat,1,mean,na.rm=TRUE)
  time_component_lwr <- apply(regional_pred_mat,1,quantile,probs=0.025)
  time_component_upr <- apply(regional_pred_mat,1,quantile,probs=0.975)
  time_component_lwr_50 <- apply(regional_pred_mat,1,quantile,probs=0.25)
  time_component_upr_50 <- apply(regional_pred_mat,1,quantile,probs=0.75)
  
  time_component_df<-data.frame(time_component_mean,
                                time_component_lwr,
                                time_component_upr,
                                time_component_lwr_50,
                                time_component_upr_50,
                                Age = Age_grid,
                                ID = "Rate of Change of Regional Component")
  names(time_component_df)<-c("RSL","lwr","upr","lwr_50","upr_50", "Age","ID")
  
  
  #-----Plotting Rate of Change for Regional------
  rate_of_change_total_plot<-ggplot()+
    geom_line(data=time_component_df,aes(x=Age*1000,y=RSL,colour = "mean"))+
    geom_ribbon(data=time_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = "CI95"),alpha=0.2)+
    geom_ribbon(data=time_component_df,
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000,fill = "CI50"),alpha=0.3)+
    theme_bw()+
    geom_hline(yintercept = 0)+
    ylab('Rate of Change (mm/yr)')+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=8),
          legend.text=element_text(size=12))+
    scale_x_continuous(breaks = round(seq(min(basis_fun_list$Age_grid*1000)-7,
                                          max(basis_fun_list$Age_grid*1000)+7, by = 100),1))+
    ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
    ggplot2::scale_fill_manual("",
                               values = c(
                                 "CI95" = ggplot2::alpha("dodgerblue1", 0.3),
                                 "CI50" = ggplot2::alpha("dodgerblue1", 0.4)
                               ),
                               labels = c(
                                 CI95 = paste0("95% Credible Interval"),
                                 CI50 = paste0("50% Credible Interval")
                               )
    ) +
    ggplot2::scale_colour_manual("",
                                 values = c(
                                   "mean" = "dodgerblue1"
                                 ),
                                 labels = c("Posterior Mean Fit")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(
        alpha = c(0.4, 0.2),
        size = 0.75
      )),
      colour = ggplot2::guide_legend(override.aes = list(
        linetype = c( 1),
        shape = c(NA),
        size = 1
      ))
    )+ 
    xlab('Year (CE)')
  rate_of_change_total_plot
  ggsave(rate_of_change_total_plot,file = "fig/all_sites/rate_of_change_regional.pdf", width = 10, height = 6)
  

  first_deriv_calc <- function(t_new){
    # Create the regional basis functions - 23 col
    B_r <- bs_bbase(t_new,xl=min(SL_df$Age),
                    xr=max(SL_df$Age), nseg = 20)
    # Now the local basis functions
    B_time <- bs_bbase(t_new,xl=min(SL_df$Age),
                       xr=max(SL_df$Age), nseg = 6, deg = 3)
    B_space_1 <- bs_bbase(SL_df$Latitude,xl=min(SL_df$Latitude),
                          xr=max(SL_df$Latitude), nseg = 6, deg = 3) 
    B_space_2 <- bs_bbase(SL_df$Longitude,xl=min(SL_df$Longitude),
                          xr=max(SL_df$Longitude), nseg = 6, deg = 3) 
    
    B_l_full <- matrix(NA, ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1), 
                       nrow = nrow(SL_df))
    regional_knots_loc <- rep(NA, ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1))
    count <- 1
    for (i in 1:ncol(B_time)) {
      for (j in 1:ncol(B_space_1)) {
        for (k in 1:ncol(B_space_2)) {
          regional_knots_loc[count] <- i
          B_l_full[, count] <- B_time[, i] * B_space_1[, j] * B_space_2[, k]
          count <- count + 1
        }
      }
    }
    
    # Get rid of all the columns which are just zero - 220 col
    B_l <- B_l_full[,-which(colSums(B_l_full) < 0.1)]
    
    # Dummy matrix for intercept & GIA
    B_h <- dbarts::makeModelMatrixFromDataFrame(data.frame(SL_df$SiteName))
    B_g <- B_h * t_new
    # Basis function matrix with B_local & B_regional
    output_B_tot <- cbind(B_h=B_h,B_g=B_g,
                          B_r = B_r,B_l = B_l)
    
    return(output_B_tot)
  }
  #-------Now create derivatives----
  h = 0.001
  t = SL_df$Age
  
  first_deriv_step1 = first_deriv_calc(t+h)
  first_deriv_step2 = first_deriv_calc(t-h)
  new_deriv_output = (first_deriv_step1 - first_deriv_step2)/(2*h)
  
  
  #--------------New basis function for site specific vertical offset----------
  B_h_deriv<- new_deriv_output[,1:87]
  h_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  for (i in 1:n_iterations){
    h_pred_mat[,i] <- B_h_deriv %*% h_intercept_post[i,]
  }
  
  h_component_mean<- apply(h_pred_mat,1,mean,na.rm=TRUE)
  h_component_lwr<- apply(h_pred_mat,1,quantile,probs=0.025)
  h_component_upr<-  apply(h_pred_mat,1,quantile,probs=0.975)
  
  #-----------New Basis Functions for Random linear local component-Site Specific slope-------
  B_deriv_g_re <- new_deriv_output[,88:174]
  g_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  for (i in 1:n_iterations){
    g_pred_mat[,i] <- B_deriv_g_re %*% b_g_post[i,]
  }
  
  g_component_mean<- apply(g_pred_mat,1,mean,na.rm=TRUE)
  g_component_lwr<- apply(g_pred_mat,1,quantile,probs=0.025)
  g_component_upr<-  apply(g_pred_mat,1,quantile,probs=0.975)
  
  #-- New basis functions for linear local component-Site Specific vertical offset---
  # Loop over all the iterations
  g_h_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol =n_iterations)
  for (i in 1:n_iterations){
    g_h_pred_mat[,i] <- B_h_deriv %*%h_intercept_post[i,]+B_deriv_g_re %*% b_g_post[i,]
  }
  
  g_h_component_mean<- rowMeans(g_h_pred_mat)
  g_h_component_lwr<- apply(g_h_pred_mat,1,quantile,probs=0.025)
  g_h_component_upr<-  apply(g_h_pred_mat,1,quantile,probs=0.975)
  g_h_component_df<-data.frame(g_h_component_mean,
                               g_h_component_lwr,
                               g_h_component_upr,
                               Age = SL_df$Age,
                               SiteName = SL_df$SiteName,
                               ID = "Linear Local Component and site specific vertical offset")
  names(g_h_component_df)<-c("RSL","lwr","upr","Age","SiteName","ID")
  
  #------------New Basis Functions for spline in time----------------
  B_deriv_t <- new_deriv_output[,175:197]
  #B_deriv_t <- new_deriv_output
  # Loop over all the iterations
  regional_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  for (i in 1:n_iterations){
    regional_pred_mat[,i] <- B_deriv_t %*% b_r_post[i,]
  }
  time_component_mean <- apply(regional_pred_mat,1,mean,na.rm=TRUE)
  time_component_lwr <- apply(regional_pred_mat,1,quantile,probs=0.025)
  time_component_upr <- apply(regional_pred_mat,1,quantile,probs=0.975)
  time_component_lwr_50 <- apply(regional_pred_mat,1,quantile,probs=0.25)
  time_component_upr_50 <- apply(regional_pred_mat,1,quantile,probs=0.75)
  
  time_component_df<-data.frame(time_component_mean,
                                time_component_lwr,
                                time_component_upr,
                                time_component_lwr_50,
                                time_component_upr_50,
                                Age = SL_df$Age,
                                SiteName = SL_df$SiteName,
                                ID = "Regional Component")
  names(time_component_df)<-c("RSL","lwr","upr",
                              "lwr_50","upr_50", "Age","SiteName","ID")
  

  #------------New Basis Functions in Space time----------------
  B_deriv_st <- new_deriv_output[,198:530]#459]
  # Loop over all the iterations
  local_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  for (i in 1:n_iterations){
    local_pred_mat[,i] <- B_deriv_st %*% b_l_post[i,]
  }
  space_time_component_mean<- apply(local_pred_mat,1,mean,na.rm=TRUE)
  space_time_component_lwr<- apply(local_pred_mat,1,quantile,probs=0.025)
  space_time_component_upr<- apply(local_pred_mat,1,quantile,probs=0.975)
  space_time_component_lwr_50 <- apply(local_pred_mat,1,quantile,probs=0.25)
  space_time_component_upr_50 <- apply(local_pred_mat,1,quantile,probs=0.75)
  space_time_component_df<-data.frame(space_time_component_mean,
                                      space_time_component_lwr,
                                      space_time_component_upr,
                                      space_time_component_lwr_50,
                                      space_time_component_upr_50,
                                      Age = SL_df$Age,
                                      SiteName = SL_df$SiteName,
                                      ID = "Local Component",
                                      data_type_id = SL_df$data_type_id)
  names(space_time_component_df)<-c("RSL","lwr","upr","lwr_50",
                                    "upr_50", "Age","SiteName","ID","data_type_id")
  
  #----Total model run-----
  # Loop over all the iterations
  mu_tot_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  for (i in 1:n_iterations){
    mu_tot_pred_mat[,i] <- B_h_deriv %*% h_intercept_post[i,] +
      B_deriv_g_re %*% b_g_post[i,] +
      B_deriv_t %*% b_r_post[i,] + B_deriv_st %*% b_l_post[i,]
  }
  mu_deriv_mean<- apply(mu_tot_pred_mat,1,mean,na.rm=TRUE)
  mu_deriv_lwr<- apply(mu_tot_pred_mat,1,quantile,probs=0.025)
  mu_deriv_upr<- apply(mu_tot_pred_mat,1,quantile,probs=0.975)
  mu_deriv_lwr_50<- apply(mu_tot_pred_mat,1,quantile,probs=0.25)
  mu_deriv_upr_50<- apply(mu_tot_pred_mat,1,quantile,probs=0.75)
  mu_deriv_component_df<-data.frame(mu_deriv_mean,
                                    mu_deriv_lwr,
                                    mu_deriv_upr,
                                    mu_deriv_lwr_50,
                                    mu_deriv_upr_50,
                                    Age = SL_df$Age,
                                    SiteName = SL_df$SiteName,
                                    data_type_id = SL_df$data_type_id,
                                    ID = "Total Component")
  names(mu_deriv_component_df)<-c("RSL",
                                  "lwr","upr","lwr_50","upr_50",
                                  "Age","SiteName","data_type_id", "ID")
  
  
  #-----Plotting Rate of Change for Total component proxy------
  rate_of_change_total_plot<-ggplot()+
    geom_line(data=subset(mu_deriv_component_df,data_type_id =="ProxyData"),
              aes(x=Age*1000,y=RSL,colour = "mean"))+
    geom_ribbon(data=subset(mu_deriv_component_df,data_type_id =="ProxyData"),
                aes(ymin=lwr,ymax=upr,x=Age*1000,fill = "CI95"),alpha=0.3)+
    geom_ribbon(data=subset(mu_deriv_component_df,data_type_id =="ProxyData"),
                aes(ymin=lwr_50,ymax=upr_50,x=Age*1000,fill = "CI50"),alpha=0.6)+
    ylab('Rate of Change (mm/yr)')+
    theme_bw()+
    theme(axis.title=element_text(size=10,face="bold"),
          axis.text=element_text(size=8),
          legend.text=element_text(size=9),
          strip.text.x = element_text(size = 8),
          strip.background =element_rect(fill=c("white")))+
    ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
    ggplot2::scale_fill_manual("",
                               values = c(
                                 "CI95" = ggplot2::alpha("purple3", 0.3),
                                 "CI50" = ggplot2::alpha("purple3", 0.6)
                               ),
                               labels = c(
                                 CI95 = paste0("95% Credible Interval"),
                                 CI50 = paste0("50% Credible Interval")
                               )
    ) +
    ggplot2::scale_colour_manual("",
                                 values = c(
                                   "mean" = "purple3"
                                 ),
                                 labels = c("Posterior Mean Fit")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(
        alpha = c( 0.6, 0.3),
        size = 0.75
      )),
      colour = ggplot2::guide_legend(override.aes = list(
        linetype = c( 1),
        shape = c( NA),
        size = 1
      ))
    ) +
    facet_wrap(~SiteName) +
    xlab('Year (CE)')+
    theme(legend.position=c(0.95,-0.05),
          legend.justification = c(1,0),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title=element_blank(),
          legend.margin=margin(c(1,1,1,1)))
  rate_of_change_total_plot
  ggsave(rate_of_change_total_plot,file = "fig/all_sites/rate_of_change_totalproxy.pdf", width = 10, height = 6)
  
  
  
  #---Rate Total plot for TG----
  mu_deriv_component_df_TG <- mu_deriv_component_df %>% 
    filter(data_type_id == "TideGaugeData") %>% 
    mutate(Age_range = max(Age*1000) - min(Age*1000))
  
  plot_rate_result_TG <- ggplot() +
    facet_wrap(~SiteName)+
    geom_line(data =  subset(mu_deriv_component_df_TG,Age_range >= 50), 
              aes(x = Age*1000, y = RSL,colour="mean")) +
    geom_ribbon(data= subset(mu_deriv_component_df_TG,Age_range >= 50),
                aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill="CI95"),alpha=0.3)+
    geom_ribbon(data= subset(mu_deriv_component_df_TG,Age_range >= 50),
                aes(y = RSL,ymin=lwr_50,ymax=upr_50,x=Age*1000,fill="CI50"),alpha=0.4)+
    xlab("Year (CE)") +
    ylab("Rate of Change (mm/year)") +
    theme_bw()+
    theme(axis.title=element_text(size=10,face="bold"),
          axis.text=element_text(size=8),
          legend.text=element_text(size=9),
          strip.text.x = element_text(size = 8),
          strip.background =element_rect(fill=c("white")))+
    #ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
    ggplot2::scale_fill_manual("",
                               values = c(
                                 "CI95" = ggplot2::alpha("purple3", 0.3),
                                 "CI50" = ggplot2::alpha("purple3", 0.4)
                               ),
                               labels = c(
                                 CI95 = paste0("95% Credible Interval"),
                                 CI50 = paste0("50% Credible Interval")
                               )
    ) +
    ggplot2::scale_colour_manual("",
                                 values = c(
                                   "mean" = "purple3"
                                 ),
                                 labels = c( "Posterior Mean Fit")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(
        alpha = c( 0.3, 0.4),
        size = 0.75
      )),
      colour = ggplot2::guide_legend(override.aes = list(
        linetype = c(1),
        shape = c(NA),
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
  plot_rate_result_TG
  ggsave(paste0("fig/",save_option,"/",
                unique(mu_deriv_component_df_TG$ID), "TG.pdf", sep = ""),
         plot_rate_result_TG, width = 10, height = 6)
  
  
  #-----Plotting Rate of Change for Total component 4 sites------
  mu_deriv_component_df_4 <-mu_deriv_component_df %>% 
    filter(SiteName %in% c( "Placentia,\n Newfoundland",
                            "East River Marsh,\n Connecticut",
                            "Cedar Island,\n North Carolina",
                            "Swan Key,\n Florida"))
  # --Filtering for 4 sites---
  SL_df_4 <-SL_df %>%   filter(SiteName %in% c("East River Marsh,\n Connecticut",
                                             "Swan Key,\n Florida",
                                             "Placentia,\n Newfoundland",
                                             "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  
  #------Ordering the sites original Data set-------
  all_data_sites<-factor(paste0(SL_df_4$Longitude, SL_df_4$Longitude), labels = 1:4)#For ordering sites
  SL_df_4 <- cbind(SL_df_4,all_data_sites=all_data_sites)
  order_sites <- SL_df_4 %>%  group_by(SiteName,all_data_sites) %>%
    dplyr::summarise(n = n()) %>%
    arrange(all_data_sites)
  mu_deriv_component_df_4$SiteName <- factor(mu_deriv_component_df_4$SiteName,levels = unique(order_sites$SiteName))

  rate_of_change_total_plot_4<-ggplot()+
    geom_line(data=mu_deriv_component_df_4,aes(x=Age*1000,y=RSL,colour = "mean"))+
    geom_ribbon(data=mu_deriv_component_df_4,aes(ymin=lwr,
                                                 ymax=upr,x=Age*1000,fill="CI95"),alpha=0.2)+
    geom_ribbon(data=mu_deriv_component_df_4,aes(ymin=lwr_50,
                                                 ymax=upr_50,x=Age*1000,fill = "CI50"),alpha=0.3)+
    theme_bw()+
    ylab('Rate of Change (mm/yr)')+
    #geom_hline(yintercept = 0)+
    theme(
          axis.title=element_text(size=12,face="bold"),
          axis.text=element_text(size=10),
          legend.text=element_text(size=10))+
    facet_wrap(~SiteName) +
    theme(strip.text.x = element_text(size = 12),
          strip.background =element_rect(fill=c("white")))+
    xlab('Year (CE)')+
    ggplot2::theme(legend.box = "horizontal", legend.position = "bottom") +
    ggplot2::scale_fill_manual("",
                               values = c(
                                 "CI95" = ggplot2::alpha("purple3", 0.2),
                                 "CI50" = ggplot2::alpha("purple3", 0.3)
                               ),
                               labels = c(
                                 CI95 = paste0("95% Credible Interval"),
                                 CI50 = paste0("50% Credible Interval")
                               )
    ) +
    ggplot2::scale_colour_manual("",
                                 values = c(
                                   "mean" = "purple3"
                                 ),
                                 labels = c("Posterior Mean Fit")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(
        alpha = c(0.3, 0.2),
        size = 0.75
      )),
      colour = ggplot2::guide_legend(override.aes = list(
        linetype = c(1),
        shape = c(NA),
        size = 1
      ))
    ) 
  
  rate_of_change_total_plot_4
  
  ggsave(rate_of_change_total_plot_4,
         file = "fig/4_sites/rate_of_change_total_4_sites.pdf", width = 10, height = 6)
  
  cat("Rate of change completed")
}
