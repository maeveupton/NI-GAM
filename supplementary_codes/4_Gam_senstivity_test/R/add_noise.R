add_noise<- function(jags.data,jags.file,SL_df,save_csv,model_run,nseg_t){
  #-----Get posterior samples for SL-----
  intercept_post <- model_run$BUGSoutput$sims.list$intercept
  b_r_post <- model_run$BUGSoutput$sims.list$b_r
  b_g_post <- model_run$BUGSoutput$sims.list$b_g
  
  pred_mean_calc <- function(t_new){
    # Create the regional basis functions
    B_r <- bs_bbase(t_new,xl=min(SL_df$Age),
                           xr=max(SL_df$Age), nseg = nseg_t)#20)
    #----Deriv----
    return(intercept_post[SL_df$SiteName] + B_r %*% colMeans(b_r_post) +
             + b_g_post[SL_df$SiteName]*(t_new))
    
  }
  #-------Now create derivatives----
  h <- 0.01
  t <- SL_df$Age
  deriv <- (pred_mean_calc(t+h) - pred_mean_calc(t-h))/(2*h)
  
 #-----Add this new term in - this is the extra standard deviation on each term----
  SL_df$NI_var_term <- sqrt(deriv^2 %*% SL_df$Age_er_average^2)[,1] 

  #----------Writing new dataframe with noisy extra column------
  write_csv(SL_df,save_csv)
  cat("Noise Added to data frame")
  return(SL_df)
  
}