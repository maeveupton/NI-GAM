run_noise_spline <- function(SL_df,
                       basis_fun_list,
                       save_option,model_run_1
                       ){
  jags_file <- "model/JAGS_mod_2_noise.jags"
  
  jags_data <- list(
    extra = SL_df$extra,
    y = SL_df$RSL,
    b_r_value = model_run_1$BUGSoutput$mean$b_regional,
    b_r_sd_value = model_run_1$BUGSoutput$sd$b_regional,
    y0_value = model_run_1$BUGSoutput$mean$intercept,
    y0_sd_value = model_run_1$BUGSoutput$sd$intercept,
    sigma_known = SL_df$RSL_er_average,
    n = nrow(SL_df),
    site = SL_df$SiteName,
    age = SL_df$Age,
    n_sites = length(unique(SL_df$SiteName)),
    B_regional = basis_fun_list$B_regional,
    B_regional_grid = basis_fun_list$B_regional_grid,
    Age_grid = basis_fun_list$Age_grid,
    n_knots_regional = ncol(basis_fun_list$B_regional),
    B_local = basis_fun_list$B_local,
    n_knots_local = ncol(basis_fun_list$B_local),
    GIA_slope = SL_df %>%
      group_by(SiteName) %>%
      slice(1) %>%
      dplyr::select(data_lm_slope) %>%
      pull(),
    GIA_slope_sd = SL_df %>%
      group_by(SiteName) %>%
      slice(1) %>%
      dplyr::select(data_lm_slope_err) %>%
      pull()
  )
  
  #----Parameters to look at----
  jags_par <- c("mu", "sigma_res", "b_regional","GIA_correction", "local",
                "regional","y_offset","regional_grid",
                "GIA_yoffset","intercept",
                "b_local","sigma_local","b_GIA", "sigma_site")
  
  
  #----Run JAGS----
  model_run_2 <- jags.parallel(
    data = jags_data,
    parameters.to.save = jags_par,
    model.file = jags_file)

  saveRDS(model_run_2,file = paste0("output/model_run_regional_noise", save_option, ".rds"))
  return("Stage 3: Model Completed with Noise")
  
}