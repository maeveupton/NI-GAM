run_spline <- function(training_set,
                       test_set,
                       training_basis_fun_list,
                       test_basis_fun_list,
                       save_option
                       ){
  jags_file <- "model/JAGS_mod_1_no_noise.jags"
  
  jags_data <- list(
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
  jags_par <- c("mu", "sigma_res", "b_regional","y_offset","GIA_yoffset",
                "GIA_correction", "regional","intercept",
                "sigma_regional", "sigmasq_all","sigmasq_all_pred","y_pred",
                "mu_pred","y_offset_pred","GIA_correction_pred","regional_pred",
                "b_GIA", "sigma_site","residuals","residuals_pred")
  
  
  #----Run JAGS----
  model_run_1 <- jags.parallel(
    data = jags_data,
    parameters.to.save = jags_par,
    model.file = jags_file)

  saveRDS(model_run_1,file = paste0("output/no_noise_no_local/model_run_no_noise_", save_option, ".rds"))
  return("Stage 1: Model Completed without Noise")
  
}