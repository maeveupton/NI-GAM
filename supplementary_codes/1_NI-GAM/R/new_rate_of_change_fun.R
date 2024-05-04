#-----------Calculating Rates of Change-------
rate_of_change_fun<- function(SL_df, model_run,save_loc,save_csv){
  #-----Get posterior samples for SL-----
  y_offset_intercept_post <- model_run$BUGSoutput$sims.list$intercept
  b_local_post <- model_run$BUGSoutput$sims.list$b_local
  b_regional_post <- model_run$BUGSoutput$sims.list$b_regional
  b_GIA_post <- model_run$BUGSoutput$sims.list$b_GIA
  n_iterations <- nrow(b_local_post)
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
    regional_pred_mat[,i] <- B_deriv_reg %*% b_regional_post[i,]
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
  write_csv(time_component_df,"data/result_csv/regional_rate.csv")
  
  
  #-----Plotting Rate of Change for Regional------
  rate_of_change_total_plot<-ggplot()+
    geom_line(data=time_component_df,aes(x=Age*1000,y=RSL),colour="#3b47ad")+
    geom_ribbon(data=time_component_df,
                aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#3b47ad",alpha=0.2)+
    geom_ribbon(data=time_component_df,aes(ymin=lwr_50,ymax=upr_50,x=Age*1000),
                fill="#3b47ad",alpha=0.3)+
    theme_bw()+
    geom_hline(yintercept = 0)+
    ylab('Rate of Change (mm/yr)')+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    theme(strip.text.x = element_text(size = 14),
          strip.background =element_rect(fill=c("white")))+
    xlab('Year(CE)')
  rate_of_change_total_plot
  ggsave(rate_of_change_total_plot,file = "fig/rate_of_change_regional.pdf", width = 10, height = 6)
  
  basis_fun_list$B_local %>% dim()
  
  # Derivative of Total Component using a time grid
  #----New Prediction Dataframe-----
  sites <- SL_df %>%
    dplyr::select(Longitude, Latitude, SiteName) %>%
    unique()
  times <- SL_df %>%
    pull(Age) %>%
    unique() %>%
    rep(nrow(sites))
  sites <- sites[rep(seq_len(nrow(sites)), each = length(SL_df %>% pull(Age) %>% unique())), ]
  to_pred <- tibble(
    Age = times,
    Longitude = sites$Longitude,
    Latitude = sites$Latitude,
    SiteName = sites$SiteName
  )
  
  
  B_regional_2 <- bs_bbase(to_pred$Age + h,xl=min(SL_df$Age),
                           xr=max(SL_df$Age), nseg = 6, deg = 2)
  B_space_1 <- bs_bbase(to_pred$Latitude,
                        xl=min(SL_df$Latitude),xr=max(SL_df$Latitude),
                        nseg = 6, deg = 2)
  B_space_2 <- bs_bbase(to_pred$Longitude,
                        xl=min(SL_df$Longitude),xr=max(SL_df$Longitude), nseg = 6, deg = 2) 
  
  B_local_deriv_full <- matrix(NA, ncol = ncol(B_regional_2) * ncol(B_space_1) * ncol(B_space_2), nrow = nrow(to_pred))
  regional_knots_loc <- rep(NA, ncol = ncol(B_regional_2) * ncol(B_space_1) * ncol(B_space_2))
  count <- 1
  for (i in 1:ncol(B_regional_2)) {
    for (j in 1:ncol(B_space_1)) {
      for (k in 1:ncol(B_space_2)) {
        regional_knots_loc[count] <- i
        B_local_deriv_full[, count] <- B_regional_2[, i] * B_space_1[, j] * B_space_2[, k]
        count <- count + 1
      }
    }
  }
  
  # Get rid of all the columns which are just zero
  #B_local_deriv <- B_local_deriv[,-which(colSums(B_local_deriv) < 0.1)]
  B_local_deriv_upr <- B_local_deriv_full[,-remove_col_index]
  
  first_deriv_calc <- function(t_new){
    # #--New grid to predict onto for derivatives--
    # new_deriv_df <- tibble(
    #   Longitude= as.numeric(SL_df$Longitude),
    #   Latitude = as.numeric(SL_df$Latitude),
    #   Age = t_new,
    #   SiteName = SL_df$SiteName)
    # N_deriv <- nrow(new_deriv_df)
    # #--Write new_df--
    # write_csv(new_deriv_df,"new_1st_deriv_df.csv")
    
    # Create the regional basis functions - 23 col
    B_regional <- bs_bbase(t_new,xl=min(SL_df$Age),
                           xr=max(SL_df$Age), nseg = 20)
    # Now the local basis functions
    B_regional_2 <- bs_bbase(t_new,xl=min(SL_df$Age),
                             xr=max(SL_df$Age), nseg = 6, deg = 2)
    B_space_1 <- bs_bbase(SL_df$Latitude,xl=min(SL_df$Latitude),
                          xr=max(SL_df$Latitude), nseg = 6, deg = 2) 
    B_space_2 <- bs_bbase(SL_df$Longitude,xl=min(SL_df$Longitude),
                          xr=max(SL_df$Longitude), nseg = 6, deg = 2) 
    
    B_local <- matrix(NA, ncol = ncol(B_regional_2) * ncol(B_space_1) * ncol(B_space_1), 
                      nrow = nrow(SL_df))
    regional_knots_loc <- rep(NA, ncol = ncol(B_regional_2) * ncol(B_space_1) * ncol(B_space_1))
    count <- 1
    for (i in 1:ncol(B_regional_2)) {
      for (j in 1:ncol(B_space_1)) {
        for (k in 1:ncol(B_space_2)) {
          regional_knots_loc[count] <- i
          B_local[, count] <- B_regional_2[, i] * B_space_1[, j] * B_space_2[, k]
          count <- count + 1
        }
      }
    }
    
    # Get rid of all the columns which are just zero - 220 col
    B_local <- B_local[,-which(colSums(B_local) < 0.1)]#Check this because order matters
    
    # Dummy matrix for intercept & GIA
    B_intercept <- dbarts::makeModelMatrixFromDataFrame(data.frame(SL_df$SiteName))
    B_GIA <- B_intercept * t_new
    # Basis function matrix with B_local & B_regional
    output_B_tot <- cbind(B_intercept=B_intercept,B_GIA=B_GIA,
                          B_regional = B_regional,B_local = B_local)
    
    return(output_B_tot)
    #return(B_regional)
    # # #----Deriv----
    # # mu_deriv <- intercept_post[SL_df$SiteName] + B_local%*%colMeans(b_local_post) + B_regional %*% colMeans(b_regional_post) +
    # #   + b_GIA_post[SL_df$SiteName]*(t_new)
    # # return(mu_deriv)
    # 
    #----Deriv----
    # mu_deriv <- matrix(NA,nrow=nrow(SL_df),ncol = n_iterations)
    # for (i in 1:n_iterations){
    #     mu_deriv[,i] <- intercept_post[SL_df$SiteName[i]]+
    #       B_local %*% b_local_post[i,]+
    #       B_regional %*% b_regional_post[i,] + b_GIA_post[SL_df$SiteName[i]]*(t_new)
    #   }
    # return(mu_deriv)
  }
  #-------Now create derivatives----
  h = 0.01
  t = SL_df$Age
  
  first_deriv_step1 = first_deriv_calc(t+h)
  first_deriv_step2 = first_deriv_calc(t-h)
  new_deriv_output = (first_deriv_step1 - first_deriv_step2)/(2*h)
  
  
  
  #--------------New basis function for site specific intercept/yoffset----------
  B_yoffset_deriv<- new_deriv_output[,1:21]
  yoffset_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  for (i in 1:n_iterations){
    yoffset_pred_mat[,i] <- B_yoffset_deriv %*% y_offset_intercept_post[i,]
  }
  
  yoffset_component_mean<- apply(yoffset_pred_mat,1,mean,na.rm=TRUE)
  yoffset_component_lwr<- apply(yoffset_pred_mat,1,quantile,probs=0.025)
  yoffset_component_upr<-  apply(yoffset_pred_mat,1,quantile,probs=0.975)
  
  #-----------New Basis Functions for Random GIA Term Site Specific slope-------
  B_deriv_GIA_re <- new_deriv_output[,22:42]
  GIA_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  for (i in 1:n_iterations){
    GIA_pred_mat[,i] <- B_deriv_GIA_re %*% b_GIA_post[i,]
  }
  
  GIA_component_mean<- apply(GIA_pred_mat,1,mean,na.rm=TRUE)
  GIA_component_lwr<- apply(GIA_pred_mat,1,quantile,probs=0.025)
  GIA_component_upr<-  apply(GIA_pred_mat,1,quantile,probs=0.975)
  
  #-- New basis functions for GIA + y offset---
  # Loop over all the iterations
  GIA_yoffset_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol =n_iterations)
  for (i in 1:n_iterations){
    GIA_yoffset_pred_mat[,i] <- B_yoffset_deriv %*%y_offset_intercept_post[i,]+B_deriv_GIA_re %*% b_GIA_post[i,]
  }
  
  GIA_yoffset_component_mean<- rowMeans(GIA_yoffset_pred_mat)
  GIA_yoffset_component_lwr<- apply(GIA_yoffset_pred_mat,1,quantile,probs=0.025)
  GIA_yoffset_component_upr<-  apply(GIA_yoffset_pred_mat,1,quantile,probs=0.975)
  GIA_yoffset_component_df<-data.frame(GIA_yoffset_component_mean,
                                       GIA_yoffset_component_lwr,
                                       GIA_yoffset_component_upr,
                                       Age = SL_df$Age,
                                       SiteName = SL_df$SiteName,
                                       ID = "y offset + GIA Component")
  names(GIA_yoffset_component_df)<-c("RSL","lwr","upr","Age","SiteName","ID")
  
  #------------New Basis Functions for spline in time----------------
  B_deriv_t <- new_deriv_output[,43:65]
  #B_deriv_t <- new_deriv_output
  # Loop over all the iterations
  regional_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  for (i in 1:n_iterations){
    regional_pred_mat[,i] <- B_deriv_t %*% b_regional_post[i,]
  }
  time_component_mean <- apply(regional_pred_mat,1,mean,na.rm=TRUE)
  time_component_lwr <- apply(regional_pred_mat,1,quantile,probs=0.025)
  time_component_upr <- apply(regional_pred_mat,1,quantile,probs=0.975)
  
  time_component_df<-data.frame(time_component_mean,
                                time_component_lwr,
                                time_component_upr,
                                Age = SL_df$Age,
                                SiteName = SL_df$SiteName,
                                ID = "Regional Component")
  names(time_component_df)<-c("RSL","lwr","upr", "Age","SiteName","ID")
  
  #----- Plotting all posterior for regional component-----
  time_component <- model_run$BUGSoutput$sims.list$regional
  post_regional <- as.data.frame(time_component[1:5,])
  rownames(post_regional) <- paste0("sample_",1:nrow(post_regional))
  post_regional_long_df <- pivot_longer(data = as.data.frame(cbind(Age = SL_df$Age,t(post_regional))),
                                        cols = starts_with("sample_"))
  # Plotting corresponding derivative
  regional_deriv_mat <- as.data.frame(regional_pred_mat[,1:5])
  colnames(regional_deriv_mat) <- paste0("sample_",1:nrow(post_regional))
  regional_deriv_mat_df <- pivot_longer(data = as.data.frame(cbind(Age = SL_df$Age,regional_deriv_mat)),
                                        cols = starts_with("sample_"))
  
  regional_post_deriv_plot <- 
    ggplot(data  = post_regional_long_df,aes(Age*1000,value*10))+
    geom_line()+
    geom_line(data  = regional_deriv_mat_df,aes(Age*1000,value),colour = "red")+
    facet_wrap(~ name)
  regional_post_deriv_plot
  ggsave(paste0("fig/",save_option,'/',"regional_deriv_post_plot.pdf", sep = ""),
         regional_post_deriv_plot, width = 10, height = 6)
  
  #------------New Basis Functions in Space time----------------
  B_deriv_st <- new_deriv_output[,66:322]
  # Loop over all the iterations
  local_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  for (i in 1:n_iterations){
    local_pred_mat[,i] <- B_deriv_st %*% b_local_post[i,]
  }
  space_time_component_mean<- apply(local_pred_mat,1,mean,na.rm=TRUE)
  space_time_component_lwr<- apply(local_pred_mat,1,quantile,probs=0.025)
  space_time_component_upr<- apply(local_pred_mat,1,quantile,probs=0.975)
  space_time_component_df<-data.frame(space_time_component_mean,
                                      space_time_component_lwr,
                                      space_time_component_upr,
                                      Age = SL_df$Age,
                                      SiteName = SL_df$SiteName,
                                      ID = "Local Component")
  names(space_time_component_df)<-c("RSL","lwr","upr", "Age","SiteName","ID")
  
  #----Total model run-----
  # Loop over all the iterations
  mu_tot_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  for (i in 1:n_iterations){
    mu_tot_pred_mat[,i] <- B_yoffset_deriv %*%y_offset_intercept_post[i,] +
      B_deriv_GIA_re %*% b_GIA_post[i,] +
      B_deriv_t %*% b_regional_post[i,] + B_deriv_st %*% b_local_post[i,]
  }
  mu_deriv_mean<- apply(mu_tot_pred_mat,1,mean,na.rm=TRUE)
  mu_deriv_lwr<- apply(mu_tot_pred_mat,1,quantile,probs=0.025)
  mu_deriv_upr<- apply(mu_tot_pred_mat,1,quantile,probs=0.975)
  mu_deriv_component_df<-data.frame(mu_deriv_mean,
                                    mu_deriv_lwr,
                                    mu_deriv_upr,
                                    Age = SL_df$Age,
                                    SiteName = SL_df$SiteName,
                                    ID = "Total Component")
  names(mu_deriv_component_df)<-c("RSL",
                                  "lwr","upr",
                                  "Age","SiteName","ID")
  
  #-----Plotting Rate of Change for Total component------
  rate_of_change_total_plot<-ggplot()+
    geom_line(data=mu_deriv_component_df,aes(x=Age*1000,y=RSL),colour="#ad4c14")+
    geom_ribbon(data=mu_deriv_component_df,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#ad4c14",alpha=0.3)+
    ggtitle("Rate of Change for Total Component")+
    theme_bw()+
    ylab('Rate of Change (mm/yr)')+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    facet_wrap(~SiteName) +
    theme(strip.text.x = element_text(size = 14),
          strip.background =element_rect(fill=c("white")))+
    xlab('Age(CE)')
  rate_of_change_total_plot
  ggsave(rate_of_change_total_plot,file = "fig/rate_of_change_total.pdf", width = 10, height = 6)
  
  #-----Plotting Rate of Change for Total component 4 sites------
  mu_deriv_component_df_4 <-mu_deriv_component_df %>% 
    filter(SiteName %in% c( "Big River Marsh,\n Newfoundland",
                            "East River Marsh,\n Connecticut",
                            "Cedar Island,\n North Carolina",
                            "Little Manatee River,\n Florida"))
  #------Ordering the sites original Data set-------
  all_data_sites<-factor(paste0(SL_df$Longitude, SL_df$Longitude), labels = 1:4)#For ordering sites
  SL_df <- cbind(SL_df,all_data_sites=all_data_sites)
  order_sites <- SL_df %>%  group_by(SiteName,all_data_sites) %>%
    dplyr::summarise(n = n()) %>%
    arrange(all_data_sites)
  mu_deriv_component_df_4$SiteName <- factor(mu_deriv_component_df_4$SiteName,levels = unique(order_sites$SiteName))
  
  rate_of_change_total_plot_4<-ggplot()+
    geom_line(data=mu_deriv_component_df_4,aes(x=Age*1000,y=RSL),colour="#ad4c14")+
    geom_ribbon(data=mu_deriv_component_df_4,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#ad4c14",alpha=0.3)+
    ggtitle("Site Specific Rate of Change")+
    theme_bw()+
    ylab('Rate of Change (mm/yr)')+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    facet_wrap(~SiteName) +
    theme(strip.text.x = element_text(size = 14),
          strip.background =element_rect(fill=c("white")))+
    xlab('Age(CE)')
  rate_of_change_total_plot_4
  
  ggsave(rate_of_change_total_plot_4,
         file = "fig/rate_of_change_total_4_sites.pdf", width = 10, height = 6)
  
  
  
  # #----Total component----
  # mu_new <- GIA_yoffset_component_mean +time_component_mean + space_time_component_mean
  # mu_lwr <- GIA_yoffset_component_mean +time_component_lwr +space_time_component_lwr
  # mu_upr <- GIA_yoffset_component_mean +time_component_upr + space_time_component_upr
  # 
  # #-----Create data frame for plotting-----
  # mod_output_df<-data.frame(mu_new,
  #                           mu_lwr,
  #                           mu_upr,
  #                           SL_df$Age,
  #                           SL_df$SiteName,
  #                           ID = "Total Posterior Model")
  # names(mod_output_df)<-c("RSL","lwr","upr",
  #                         "Age",
  #                         "SiteName",
  #                         "ID")
  
  
  # #---- Reording sites for plots----
  # all_data_sites<-factor(paste0(SL_df$Longitude, SL_df$Longitude), labels = 1:16)#For ordering sites
  # SL_df <- cbind(SL_df,all_data_sites)
  # order_sites <- SL_df %>%  group_by(SiteName,all_data_sites) %>%
  #   dplyr::summarise(n = n()) %>%
  #   arrange(all_data_sites)
  # SL_df$SiteName <- factor(SL_df$SiteName,levels = unique(order_sites$SiteName))
  
  #-----Plotting Rate of Change for regional component------
  rate_of_change_regional_plot<-ggplot()+
    geom_line(data=time_component_df,aes(x=Age*1000,y=RSL),colour="#3b47ad")+
    geom_ribbon(data=time_component_df,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#3b47ad",alpha=0.3)+
    ggtitle("Rate of Change for Regional Component")+
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))+
    ylab('Rate of Change (mm/yr)')+
    #facet_wrap(~SiteName)+
    xlab('Age(CE)')
  rate_of_change_regional_plot
  ggsave(rate_of_change_regional_plot,file = "fig/rate_of_change_regional.pdf", width = 10, height = 6)
  
  # #-----Plotting Rate of Change for Local component------
  # rate_of_change_local_plot<-ggplot()+
  #   geom_line(data=space_time_component_df,aes(x=Age*1000,y=RSL),colour="#ad4c14")+
  #   geom_ribbon(data=space_time_component_df,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#ad4c14",alpha=0.3)+
  #   ggtitle("Rate of Change for Local Component")+
  #   theme_bw()+
  #   ylab('Rate of Change (mm/yr)')+
  #   theme(plot.title = element_text(size=22),
  #         axis.title=element_text(size=14,face="bold"),
  #         axis.text=element_text(size=12),
  #         legend.text=element_text(size=12))+
  #   facet_wrap(~SiteName) +
  #   theme(strip.text.x = element_text(size = 14),
  #         strip.background =element_rect(fill=c("white")))+
  #   xlab('Age(CE)')
  # rate_of_change_local_plot
  # ggsave(rate_of_change_local_plot,file = "fig/rate_of_change_local.pdf", width = 10, height = 6)
  # 
  # #-----Plotting Rate of Change for GIA+yoffset component------
  # rate_of_change_GIA_plot<-ggplot()+
  #   geom_line(data=GIA_yoffset_component_df,aes(x=Age*1000,y=RSL),colour="red")+
  #   geom_ribbon(data=GIA_yoffset_component_df,aes(ymin=lwr,
  #                                                 ymax=upr,x=Age*1000),colour="red",alpha=0.1)+
  #   ggtitle("Rate of Change for y offset + GIA Component")+
  #   theme_bw()+
  #   ylab('Rate of Change (mm/yr)')+
  #   facet_wrap(~SiteName)+
  #   xlab('Age(CE)')
  # rate_of_change_GIA_plot
  # ggsave(rate_of_change_GIA_plot,file = "fig/rate_of_change_GIA_yoffset.pdf", width = 10, height = 6)
  
  
  cat("Rate of change completed")
}
