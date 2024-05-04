spline_basis_function<- function(SL_df,nseg_t,nseg_st){
  #-----------------------Basis functions-----------------------
  # Create the regional basis functions
  B_r <- bs_bbase(SL_df$Age,xl=min(SL_df$Age),xr=max(SL_df$Age), nseg = nseg_t)
  # Create the regional grid basis functions
  Age_grid <- seq(min(SL_df$Age),max(SL_df$Age), by = 0.01)# every 100 years
  B_r_grid <- bs_bbase(Age_grid,
                              xl=min(SL_df$Age),xr=max(SL_df$Age),nseg = nseg_t)
  # Plot the B-splines
  plot(Age_grid, B_r_grid[,1], ylim = range(B_r_grid), type = 'n',
       ylab = "B_r")
  for(j in 1:ncol(B_r_grid)) {
    lines(Age_grid, B_r_grid[,j], col = j)
  }

  # Now the local basis functions
  B_time <- bs_bbase(SL_df$Age,xl=min(SL_df$Age),xr=max(SL_df$Age), nseg = nseg_st, deg = 2)
  B_space_1 <- bs_bbase(SL_df$Latitude,xl=min(SL_df$Latitude),xr=max(SL_df$Latitude), nseg = nseg_st, deg = 2) # Put ranges on these with xl and xr for later ease of interpolation
  B_space_2 <- bs_bbase(SL_df$Longitude,xl=min(SL_df$Longitude),xr=max(SL_df$Longitude), nseg = nseg_st, deg = 2) # Put ranges on these with xl and xr for later ease of interpolation
  
  B_l_full <- matrix(NA, ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1), nrow = nrow(SL_df))
  regional_knots_loc <- rep(NA, ncol = ncol(B_time) * ncol(B_space_1) * ncol(B_space_1))
  count <- 1
  for (i in 1:ncol(B_time)) {
    for (j in 1:ncol(B_space_1)) {
      for (k in 1:ncol(B_space_2)) {
        regional_knots_loc[count] <- i
        B_l_full[, count] <- B_time[, i] * B_space_1[, j] * B_space_2[, k]
        count <- count + 1
      }
    }
  }

  # Get rid of all the columns which are just zero
  B_l <- B_l_full[,-which(colSums(B_l_full) < 0.1)]
  
  # Find the index here that you remove then use this in the derivative
  remove_col_index <- which(colSums(B_l_full) < 0.1)
  
  
  # All Basis Functions
  basis_list <- list(
    Age_grid = Age_grid,
    B_r=B_r, 
    B_l=B_l,
    B_r_grid=B_r_grid)
  cat("Basis functions created for the regional and local splines")
  return(basis_list)
}