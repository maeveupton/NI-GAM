linear_reg_rates<- function(SL_df){#, save_loc){
  SL_df_filter <- SL_df %>%  
    filter(!Age > 1.800) # Ignoring recent human influences to SL rise
  
  # Doing linear regression on rest of data
  SL_df_lm <- SL_df_filter %>% 
    group_by(SiteName) %>% 
    mutate(
      data_lm_slope = round(lm(RSL ~ Age)$coefficients[[2]], 2),
      data_lm_slope_err = summary(lm(RSL ~ Age))$coefficients[2, 4]
    )
  
 
  # Table of GIA rate vs lm rate from proxy data
 lm_slopes <- SL_df_lm %>%
   dplyr::select(SiteName, data_lm_slope,data_lm_slope_err) %>%
    unique()
 
   #write_csv(GIA_vs_lm_slopes,sav_loc)
  cat("Linear Regression on the data excluding tide gauges \n")
  return(lm_slopes)
  
}