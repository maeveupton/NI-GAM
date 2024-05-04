spline_basis_function<- function(SL_df,interval){
  #-----------------------Basis functions-----------------------
  # Create the regional basis functions
  B_regional <- bs_bbase(SL_df$Age,xl=min(interval$Age),xr=max(interval$Age), nseg = 20)
  
  # Now the local basis functions
  B_regional_2 <- bs_bbase(SL_df$Age,xl=min(interval$Age),xr=max(interval$Age), nseg = 6, deg = 2)
  B_space_1 <- bs_bbase(SL_df$Latitude,xl=min(interval$Latitude),xr=max(interval$Latitude), nseg = 6, deg = 2) # Put ranges on these with xl and xr for later ease of interpolation
  B_space_2 <- bs_bbase(SL_df$Longitude,xl=min(interval$Longitude),xr=max(interval$Longitude), nseg = 6, deg = 2) # Put ranges on these with xl and xr for later ease of interpolation
  
  B_local <- matrix(NA, ncol = ncol(B_regional_2) * ncol(B_space_1) * ncol(B_space_1), nrow = nrow(SL_df))
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
  #B_local <- B_local[,-which(colSums(B_local) < 0.1)]
  
  # # ---- BASIS Function for derivative of regional----
  # # Check to sort the time points
  # Age_grid <- seq(min(SL_df$Age),max(SL_df$Age), by = 0.01)# every 100 years
  # B_regional_grid <- bs_bbase(Age_grid,
  #                             xl=min(SL_df$Age),xr=max(SL_df$Age),nseg = 20)
  
  # All Basis Functions
  basis_list <- list(
    #Age_grid = Age_grid,
    B_regional=B_regional, 
    B_local=B_local
    #B_regional_grid=B_regional_grid
    )
  cat("Basis functions created for the regional and local splines")
  return(basis_list)
  
  
}