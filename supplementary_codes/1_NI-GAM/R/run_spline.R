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
    B_regional = basis_fun_list$B_regional,
    B_regional_grid = basis_fun_list$B_regional_grid,
    Age_grid = basis_fun_list$Age_grid,
    n_knots_regional = ncol(basis_fun_list$B_regional),
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
  jags_par <- c("mu", "sigma_res", "b_regional","y_offset","GIA_yoffset",
                "GIA_correction","regional","intercept",
                "sigma_regional","sigmasq_all","regional_grid",
                "b_GIA", "sigma_site","residuals")
  
  
  #----Run JAGS----
  model_run_1 <- jags.parallel(
    data = jags_data,
    parameters.to.save = jags_par,
    model.file = jags_file)

  saveRDS(model_run_1,file = paste0("output/model_run_no_local_no_noise", save_option, ".rds"))
  return("Stage 1: Model Completed without Noise")
  
}