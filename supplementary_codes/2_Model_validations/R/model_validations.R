model_validations <- function(SL_df,model_run,test_set){
  #-----Get posterior samples for SL-----
  y_offset_intercept_post <- model_run$BUGSoutput$sims.list$intercept
  b_local_post <- model_run$BUGSoutput$sims.list$b_local
  b_regional_post <- model_run$BUGSoutput$sims.list$b_regional
  b_GIA_post <- model_run$BUGSoutput$sims.list$b_GIA
  n_iterations <- nrow(b_local_post)
  
  # Prediction dataframe :   
  head(test_set)
  
  # Predicted basis functions
  B_p_regional <- bs_bbase(test_set$Age,xl=min(SL_df$Age),xr=max(SL_df$Age), nseg = 20)
  
  
  # Now the local basis functions
  B_regional_2 <- bs_bbase(test_set$Age,xl=min(SL_df$Age),xr=max(SL_df$Age), nseg = 6, deg = 2)
  B_space_1 <- bs_bbase(test_set$Latitude,xl=min(SL_df$Latitude),xr=max(SL_df$Latitude), nseg = 6, deg = 2) # Put ranges on these with xl and xr for later ease of interpolation
  B_space_2 <- bs_bbase(test_set$Longitude,xl=min(SL_df$Longitude),xr=max(SL_df$Longitude), nseg = 6, deg = 2) # Put ranges on these with xl and xr for later ease of interpolation
  
  B_p_local <- matrix(NA, ncol = ncol(B_regional_2) * ncol(B_space_1) * ncol(B_space_1), nrow = nrow(test_set))
  regional_knots_loc <- rep(NA, ncol = ncol(B_regional_2) * ncol(B_space_1) * ncol(B_space_1))
  count <- 1
  for (i in 1:ncol(B_regional_2)) {
    for (j in 1:ncol(B_space_1)) {
      for (k in 1:ncol(B_space_2)) {
        regional_knots_loc[count] <- i
        B_p_local[, count] <- B_regional_2[, i] * B_space_1[, j] * B_space_2[, k]
        count <- count + 1
      }
    }
  }
  
  # Get rid of all the columns which are just zero
  #B_p_local <- B_p_local[,-which(colSums(B_p_local) < 0.1)]
  # Dummy matrix for intercept & GIA
  B_intercept_p <- dbarts::makeModelMatrixFromDataFrame(data.frame(test_set$SiteName))
  B_GIA_p <- B_intercept_p * test_set$Age
  
  # #--------------New basis function for site specific intercept/yoffset----------
  # yoffset_pred_mat <- matrix(NA,nrow=nrow(test_set),ncol = 10)
  # for (i in 1:10){
  #   yoffset_pred_mat[,i] <- B_intercept_p %*% y_offset_intercept_post[i,]
  # }
  # 
  # yoffset_component_mean<- apply(yoffset_pred_mat,1,mean,na.rm=TRUE)
  # yoffset_component_lwr<- apply(yoffset_pred_mat,1,quantile,probs=0.025)
  # yoffset_component_upr<-  apply(yoffset_pred_mat,1,quantile,probs=0.975)
  # 
  # #-----------New Basis Functions for Random GIA Term Site Specific slope-------
  # GIA_pred_mat <- matrix(NA,nrow=nrow(test_set),ncol = 10)
  # for (i in 1:10){
  #   GIA_pred_mat[,i] <- B_GIA_p %*% b_GIA_post[i,]
  # }
  # 
  # GIA_component_mean<- apply(GIA_pred_mat,1,mean,na.rm=TRUE)
  # GIA_component_lwr<- apply(GIA_pred_mat,1,quantile,probs=0.025)
  # GIA_component_upr<-  apply(GIA_pred_mat,1,quantile,probs=0.975)
  # 
  # #-- New basis functions for GIA + y offset---
  # # Loop over all the iterations
  # GIA_yoffset_pred_mat <- matrix(NA,nrow=nrow(test_set),ncol =10)
  # for (i in 1:10){
  #   GIA_yoffset_pred_mat[,i] <- B_intercept_p %*%y_offset_intercept_post[i,]+B_GIA_p %*% b_GIA_post[i,]
  # }
  # 
  # GIA_yoffset_component_mean<- rowMeans(GIA_yoffset_pred_mat)
  # GIA_yoffset_component_lwr<- apply(GIA_yoffset_pred_mat,1,quantile,probs=0.025)
  # GIA_yoffset_component_upr<-  apply(GIA_yoffset_pred_mat,1,quantile,probs=0.975)
  # GIA_yoffset_component_df<-data.frame(GIA_yoffset_component_mean,
  #                                      GIA_yoffset_component_lwr,
  #                                      GIA_yoffset_component_upr,
  #                                      Age = test_set$Age,
  #                                      SiteName = test_set$SiteName,
  #                                      ID = "y offset + GIA Component")
  # names(GIA_yoffset_component_df)<-c("RSL","lwr","upr","Age","SiteName","ID")
  # 
  # #------------New Basis Functions for spline in time----------------
  # # Loop over all the iterations
  # regional_pred_mat <- matrix(NA,nrow=nrow(test_set),ncol = 10)
  # for (i in 1:10){
  #   regional_pred_mat[,i] <- B_p_regional %*% b_regional_post[i,]
  # }
  # time_component_mean <- apply(regional_pred_mat,1,mean,na.rm=TRUE)
  # time_component_lwr <- apply(regional_pred_mat,1,quantile,probs=0.025)
  # time_component_upr <- apply(regional_pred_mat,1,quantile,probs=0.975)
  # 
  # time_component_df<-data.frame(time_component_mean,
  #                               time_component_lwr,
  #                               time_component_upr,
  #                               Age = test_set$Age,
  #                               SiteName = test_set$SiteName,
  #                               ID = "Regional Component")
  # names(time_component_df)<-c("RSL","lwr","upr", "Age","SiteName","ID")
  # 
  # 
  # #------------New Basis Functions in Space time----------------
  # # Loop over all the iterations
  # local_pred_mat <- matrix(NA,nrow=nrow(test_set),ncol = 10)
  # for (i in 1:10){
  #   local_pred_mat[,i] <- B_p_local %*% b_local_post[i,]
  # }
  # space_time_component_mean<- apply(local_pred_mat,1,mean,na.rm=TRUE)
  # space_time_component_lwr<- apply(local_pred_mat,1,quantile,probs=0.025)
  # space_time_component_upr<- apply(local_pred_mat,1,quantile,probs=0.975)
  # space_time_component_df<-data.frame(space_time_component_mean,
  #                                     space_time_component_lwr,
  #                                     space_time_component_upr,
  #                                     Age = test_set$Age,
  #                                     SiteName = test_set$SiteName,
  #                                     ID = "Local Component")
  # names(space_time_component_df)<-c("RSL","lwr","upr", "Age","SiteName","ID")
  
  #----Total model run-----
  # Loop over all the iterations
  # mu_tot_pred_mat <- matrix(NA,nrow=nrow(new_deriv_output),ncol = n_iterations)
  # for (i in 1:n_iterations){
  #   mu_tot_pred_mat[,i] <- B_intercept_p %*%y_offset_intercept_post[i,] +
  #     B_GIA_p %*% b_GIA_post[i,] +
  #     B_p_regional %*% b_regional_post[i,] + B_p_local %*% b_local_post[i,]
  # }
  # mu_deriv_mean<- apply(mu_tot_pred_mat,1,mean,na.rm=TRUE)
  # mu_deriv_lwr<- apply(mu_tot_pred_mat,1,quantile,probs=0.025)
  # mu_deriv_upr<- apply(mu_tot_pred_mat,1,quantile,probs=0.975)
  # mu_deriv_component_df<-data.frame(mu_deriv_mean,
  #                                   mu_deriv_lwr,
  #                                   mu_deriv_upr,
  #                                   Age = test_set$Age,
  #                                   SiteName = test_set$SiteName,
  #                                   ID = "Total Component")
  # names(mu_deriv_component_df)<-c("RSL",
  #                                 "lwr","upr",
  #                                 "Age","SiteName","ID")
  
  #-----Plotting True vs predicted------ 
  mu_pred <- B_intercept_p[test_set$SiteName] + B_p_regional %*% colMeans(b_regional_post) +
    + b_GIA_post[test_set$SiteName]*(test_set$Age) + B_p_local %*% colMeans(b_local_post)
  
  true <- SL_df$RSL
  pred <- mu_deriv_component_df$RSL
  pred_true_df <- data.frame(true,pred)
  plot(pred,true)
  
  true_pred_plot<-ggplot()+
    #geom_line(data=mu_deriv_component_df,aes(x=Age*1000,y=RSL),colour="#ad4c14")+
    geom_point (data = pred_true_df, aes(x = true, y = pred))+
    #geom_ribbon(data=mu_deriv_component_df,aes(ymin=lwr,ymax=upr,x=Age*1000),fill="#ad4c14",alpha=0.3)+
    ggtitle("True vs Predicted")+
    theme_bw()+
    theme(plot.title = element_text(size=22),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12),
          legend.text=element_text(size=12))
    #facet_wrap(~SiteName) +
  true_pred_plot
  ggsave(true_pred_plot,file = "fig/true_pred.pdf", width = 10, height = 6)
  
  

  
  
  # residual = observed - estimated
  #my_loo_normal_95_PI$residual <- my_loo_normal_95_PI$observed_proportion - my_loo_normal_95_PI$median_p
  model_residual <- model_run$BUGSoutput$sims.list$residuals
  model_residual
  # RMSE, MSE --------------------------------------
  model_residual_MSE <- mean(model_residual$residual^2) 
  model_residual_MSE
  #model_residual_MSE_sites <- model_residual %>% group_by(SiteName) %>% summarise(MSE = mean(residual^2)) 
  model_residual_RMSE <- sqrt(mean(model_residual^2))
  model_residual_RMSE
  #model_residual_RMSE_sites <- my_loo_normal_95_PI %>% group_by(sector_category) %>% summarise(RMSE = sqrt(mean(residual^2))) 

  # ## Prediction intervals ---------------------------
  # se_prop <- as.matrix(loo_test_data[,c("Public.SE", "Commercial_medical.SE")])
  # match_method_test <- loo_test_data$index_method
  # match_country_test <- loo_test_data$index_country
  # match_years_test <- loo_test_data$index_year
  # n_sims <- dim(post_samps$P)[1] # 6000
  # y.rep <- array (NA, dim = c(n_sims, 3, max(match_method_test),max(match_country_test),max(match_years_test)))
  # for(i in 1:n_sims) {
  #   for(j in 1:length(match_country_test)) {
  #     print(i)
  #     print(j)
  #     y.rep[i,1,match_method_test[j],match_country_test[j],match_years_test[j]] <- rnorm(1, 
  #                                                                                        mean = post_samps$P[i,1,match_method_test[j],match_country_test[j],match_years_test[j]],
  #                                                                                        sd = as.vector(unlist(se_prop[j,1])))
  #     y.rep[i,2,match_method_test[j],match_country_test[j],match_years_test[j]] <- rnorm(1, 
  #                                                                                        mean = post_samps$P[i,2,match_method_test[j],match_country_test[j],match_years_test[j]],
  #                                                                                        sd = as.vector(unlist(se_prop[j,2])))
  #     y.rep[i,3,match_method_test[j],match_country_test[j],match_years_test[j]] <- 1 - (y.rep[i,1,match_method_test[j],match_country_test[j],match_years_test[j]] + y.rep[i,2,match_method_test[j],match_country_test[j],match_years_test[j]])
  #   }
  # }
  # 
  # # Median PI width  
  # my_loo_normal_95_PI <- my_loo_normal_95_PI %>% mutate(ci_width = upper_q_P - lower_q_P)
  # my_loo_normal_95_PI %>% group_by(sector_category) %>% summarise(Median=median(ci_width, na.rm=TRUE)) 
  # # Mean error
  # mean_error <- my_loo_normal_95_PI %>% group_by(sector_category) %>% summarise(mean_error = mean(residual)) 
  # # Median error
  # med_error <- my_loo_normal_95_PI %>% group_by(sector_category) %>% summarise(med_error = median(residual)) 
  # # Median absolute error
  # medabs_error <- my_loo_normal_95_PI %>% group_by(sector_category) %>% summarise(medabs_error = median(abs(residual))) 
}