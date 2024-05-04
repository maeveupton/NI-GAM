run_noise_spline <- function(SL_df,
                       basis_fun_list,
                       save_option,model_run_1
                       ){
  jags_file <- "model/JAGS_mod_2_noise.jags"
  
  jags_data <- list(
    NI_var_term = SL_df$NI_var_term,
    y = SL_df$RSL,
    age = SL_df$Age,
    b_r_value = model_run_1$BUGSoutput$mean$b_r,
    b_r_sd_value = model_run_1$BUGSoutput$sd$b_r,
    h_value = model_run_1$BUGSoutput$mean$intercept,
    h_sd_value = model_run_1$BUGSoutput$sd$intercept,
    sigma_known = SL_df$RSL_er_average,
    n = nrow(SL_df),
    site = SL_df$SiteName,
    n_sites = length(unique(SL_df$SiteName)),
    B_r = basis_fun_list$B_r,
    B_r_grid = basis_fun_list$B_r_grid,
    n_knots_r = ncol(basis_fun_list$B_r),
    B_l = basis_fun_list$B_l,
    n_knots_l = ncol(basis_fun_list$B_l),
    known_rate = SL_df %>%
      group_by(SiteName) %>%
      slice(1) %>%
      dplyr::select(data_lm_slope) %>%
      pull(),
    known_rate_err = SL_df %>%
      group_by(SiteName) %>%
      slice(1) %>%
      dplyr::select(data_lm_slope_err) %>%
      pull()
  )
  
  #----Parameters to look at----
  jags_par <- c("mu", "sigma_res", "b_r","g_z_x", "l",
                "r","h_z_x","r_grid",
                "g_h_z_x","intercept",
                "b_l","sigma_l","b_g")
  
  
  #----Run JAGS----
  model_run_2 <- jags.parallel(
    data = jags_data,
    n.iter=5000,
    parameters.to.save = jags_par,
    model.file = jags_file)

  saveRDS(model_run_2,file = paste0("output/model_run_regional_noise", save_option, ".rds"))
  return("Stage 3: Model Completed with Noise")
  
}