spline_basis_function<- function(SL_df){
  #-----------------------Basis functions-----------------------
  # Create the regional basis functions
  B_regional <- bs_bbase(SL_df$Age,xl=min(SL_df$Age),xr=max(SL_df$Age), nseg = 20)
  
  # Now the local basis functions
  B_regional_2 <- bs_bbase(SL_df$Age,xl=min(SL_df$Age),xr=max(SL_df$Age), nseg = 6, deg = 2)
  B_space_1 <- bs_bbase(SL_df$Latitude,xl=min(SL_df$Latitude),xr=max(SL_df$Latitude), nseg = 6, deg = 2) # Put ranges on these with xl and xr for later ease of interpolation
  B_space_2 <- bs_bbase(SL_df$Longitude,xl=min(SL_df$Longitude),xr=max(SL_df$Longitude), nseg = 6, deg = 2) # Put ranges on these with xl and xr for later ease of interpolation
  
  B_local_full <- matrix(NA, ncol = ncol(B_regional_2) * ncol(B_space_1) * ncol(B_space_1), nrow = nrow(SL_df))
  regional_knots_loc <- rep(NA, ncol = ncol(B_regional_2) * ncol(B_space_1) * ncol(B_space_1))
  count <- 1
  for (i in 1:ncol(B_regional_2)) {
    for (j in 1:ncol(B_space_1)) {
      for (k in 1:ncol(B_space_2)) {
        regional_knots_loc[count] <- i
        B_local_full[, count] <- B_regional_2[, i] * B_space_1[, j] * B_space_2[, k]
        count <- count + 1
      }
    }
  }

  # Get rid of all the columns which are just zero
  B_local <- B_local_full[,-which(colSums(B_local_full) < 0.1)]
  
  # Find the index here that you remove then use this in the derivative
  remove_col_index <- which(colSums(B_local_full) < 0.1)
  
  # ---- BASIS Function for derivative of regional----
  # NOT Right
  # Check to sort the time points
  Age_grid <- seq(min(SL_df$Age),max(SL_df$Age), by = 0.01)# every 100 years
  B_regional_grid <- bs_bbase(Age_grid,
                              xl=min(SL_df$Age),xr=max(SL_df$Age),nseg = 20)
  
  # All Basis Functions
  basis_list <- list(
    Age_grid = Age_grid,
    B_regional=B_regional, 
    B_local=B_local,
    B_regional_grid=B_regional_grid)
  cat("Basis functions created for the regional and local splines")
  return(basis_list)
  
  # # Now the local basis functions
  # B_regional_grid_2 <- bs_bbase(Age_grid,xl=min(SL_df$Age),xr=max(SL_df$Age), deg = 20)
  # B_space_grid_1 <- bs_bbase(SL_df$Latitude,xl=min(SL_df$Latitude),xr=max(SL_df$Latitude), nseg = 6, deg = 2) # Put ranges on these with xl and xr for later ease of interpolation
  # B_space_grid_2 <- bs_bbase(SL_df$Longitude,xl=min(SL_df$Longitude),xr=max(SL_df$Longitude), nseg = 6, deg = 2) # Put ranges on these with xl and xr for later ease of interpolation
  # 
  # B_local_grid <- matrix(NA, ncol = ncol(B_regional_grid_2) * ncol(B_space_grid_1) * ncol(B_space_grid_1), nrow = nrow(SL_df))
  # regional_knots_loc <- rep(NA, ncol = ncol(B_regional_grid_2) * ncol(B_space_grid_1) * ncol(B_space_grid_1))
  # count <- 1
  # for (i in 1:ncol(B_regional_grid_2)) {
  #   for (j in 1:ncol(B_space_grid_1)) {
  #     for (k in 1:ncol(B_space_grid_2)) {
  #       regional_knots_loc[count] <- i
  #       B_local[, count] <- B_regional_grid_2[, i] * B_space_grid_1[, j] * B_space_grid_2[, k]
  #       count <- count + 1
  #     }
  #   }
  # }
  # 
  # # Get rid of all the columns which are just zero
  # B_local_grid <- B_local_grid[,-which(colSums(B_local_grid) < 0.1)]
  
  
}