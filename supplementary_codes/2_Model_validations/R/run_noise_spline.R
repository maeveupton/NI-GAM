run_noise_spline <- function(training_set,
                             test_set,
                             training_basis_fun_list,
                             test_basis_fun_list,
                             save_option,
                             model_run_1
                       ){
  jags_file <- "model/JAGS_mod_2_noise.jags"
  
  
  jags_data <- list(
    extra = training_set$extra,
    extra_pred = test_set$extra,
    b_r_value = model_run_1$BUGSoutput$mean$b_regional,
    b_r_sd_value = model_run_1$BUGSoutput$sd$b_regional,
    y0_value = model_run_1$BUGSoutput$mean$intercept,
    y0_sd_value = model_run_1$BUGSoutput$sd$intercept,
    y = training_set$RSL,
    sigma_known = training_set$RSL_er_average,
    sigma_known_pred = test_set$RSL_er_average,
    n = nrow(training_set),
    n_pred = nrow(test_set),
    site_pred = test_set$SiteName,
    age_pred = test_set$Age,
    site = training_set$SiteName,
    age = training_set$Age,
    n_sites = length(unique(training_set$SiteName)),
    n_sites_pred = length(unique(test_set$SiteName)),
    B_regional = training_basis_fun_list$B_regional,
    B_regional_pred = test_basis_fun_list$B_regional,
    n_knots_regional = ncol(training_basis_fun_list$B_regional),
    B_local = training_basis_fun_list$B_local,
    B_local_pred = test_basis_fun_list$B_local,
    n_knots_local = ncol(training_basis_fun_list$B_local),
    GIA_slope = training_set %>%
      group_by(SiteName) %>%
      slice(1) %>%
      dplyr::select(data_lm_slope) %>%
      pull(),
    GIA_slope_sd = training_set %>%
      group_by(SiteName) %>%
      slice(1) %>%
      dplyr::select(data_lm_slope_err) %>%
      pull()
  )
  
  
  #----Parameters to look at----
  jags_par <- c("mu", "sigma_res", "b_regional",
                "GIA_correction", "local", "regional","y_offset",
                "GIA_yoffset","intercept","residuals",#"regional_grid",
                "b_local",
                "sigmasq_all_pred","y_pred","mu_pred",
                "y_offset_pred", "GIA_correction_pred","regional_pred","local_pred",
                "sigma_local","residuals_pred",
                "b_GIA", "sigma_site")
  
  
  #----Run JAGS----
  model_run_2 <- jags.parallel(
    data = jags_data,
    parameters.to.save = jags_par,
    model.file = jags_file)

  saveRDS(model_run_2,file = paste0("output/noise_full_model/model_run_regional_noise_", save_option, ".rds"))
  return("Stage 3: Model Completed with Noise")
  
}