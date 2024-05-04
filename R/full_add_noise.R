add_noise<- function(jags_data,SL_df,save_csv,model_run){
  #-----Get posterior samples for SL-----
  intercept_post <- model_run$BUGSoutput$sims.list$intercept
  n_iterations <- nrow(intercept_post)
  b_local_post <- matrix(0,n_iterations,jags_data$n_knots_local)#model_run$BUGSoutput$sims.list$b_local
  b_regional_post <- model_run$BUGSoutput$sims.list$b_regional
  b_GIA_post <- model_run$BUGSoutput$sims.list$b_GIA
  
  pred_mean_calc <- function(t_new){
    # Create the regional basis functions
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
    
    # Get rid of all the columns which are just zero
    B_local <- B_local[,-which(colSums(B_local) < 0.1)]
    
    #----Deriv----
    return(intercept_post[SL_df$SiteName] + B_local%*%colMeans(b_local_post) + B_regional %*% colMeans(b_regional_post) +
             + b_GIA_post[SL_df$SiteName]*(t_new))
    
  }
  #-------Now create derivatives----
  h <- 0.01
  t <- SL_df$Age
  deriv <- (pred_mean_calc(t+h) - pred_mean_calc(t-h))/(2*h)
  plot(t,deriv)
  SL_df_test <- cbind(SL_df, deriv=deriv)
  ggplot(SL_df_test)+
    geom_line(aes(x = Age*1000,y = deriv))+
    facet_wrap(~SiteName)+
    ylab("Rate of change (mm/yr)")+
    xlab("Age (CE)")
  #-----Add this new term in - this is the extra standard deviation on each term----
  SL_df$extra <- sqrt(deriv^2 %*% SL_df$Age_er_average^2)[,1] # 0 years for tide gauge
  SL_df$extra[which(SL_df$data_type_id == "TideGaugeData")] <- 0 ## Added by NC
  #----------Writing new dataframe with noisy extra column------
  write_csv(SL_df,save_csv)
  cat("Noise Added to data frame")
  return(SL_df)
  
}