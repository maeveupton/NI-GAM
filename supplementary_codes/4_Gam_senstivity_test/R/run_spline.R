run_spline <- function(SL_df,
                       basis_fun_list,
                       save_option
                       ){
  jags_file <- "model/JAGS_mod_1_no_noise.jags"
  
  jags_data <- list(
    y = SL_df$RSL,
    sigma_known = SL_df$RSL_er_average,
    n = nrow(SL_df),
    site = SL_df$SiteName,
    age = SL_df$Age,
    n_sites = length(unique(SL_df$SiteName)),
    B_r = basis_fun_list$B_r,
    B_r_grid = basis_fun_list$B_r_grid,
    n_knots_r = ncol(basis_fun_list$B_r),
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
  jags_par <- c("mu", "sigma_res", "b_r","h_z_x","g_h_z_x",
                "g_z_x","r","intercept",
                "sigma_r","sigmasq_all","r_grid",
                "b_g", "sigma_h","residuals_mod1")
  
  
  #----Run JAGS----
  model_run_1 <- jags.parallel(
    data = jags_data,
    parameters.to.save = jags_par,
    model.file = jags_file,
    n.chains=3, n.iter=5000, n.burnin=1000,n.thin = 5)

  saveRDS(model_run_1,file = paste0("output/model_run_no_local_no_noise", save_option, ".rds"))
  return("Stage 1: Model Completed without Noise")
  
}