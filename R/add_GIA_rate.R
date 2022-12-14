add_GIA_rate<- function(SL_df,save_loc){
  #--------------------GIA DATA from Peltier Website ICE5G--------------
  download.file('https://www.dropbox.com/s/a8m4kzmcxh54mmt/dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc?dl=1',
                destfile = "dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc",
                method = "libcurl", mode = "wb")
  nc_data <- nc_open("dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc")#ICE5G: better fit for data

  #---Roudning to 1 decimal point to reduce number of spatial options--
  dat_lon <- round(SL_df$Longitude,1)
  dat_lat <- round(SL_df$Latitude,1)
  
  ### get lon and lat from GIA model output
  lon <- round(ncvar_get(nc_data, "Lon"),1)
  lat <- round(ncvar_get(nc_data, "Lat"),1)
  
  ### needs to be sorted for the match.closest function() below
  ### note need the index for the unsorted coordinate for pulling the correct SL rate later
  gia_lat <- tibble(index = 1:length(lat),lat) %>% arrange(lat)
  gia_lon <- tibble(index = 1:length(lon),lon) %>% arrange(lon)
  
  #--- Matching closest long & lat values---
  lat_index <- gia_lat$index[match.closest(dat_lat, gia_lat$lat)]
  lon_index <- gia_lon$index[match.closest(360+dat_lon, gia_lon$lon)] # change data lon to degrees east
  
  #---GIA rates---
  SL <- ncvar_get(nc_data, "Dsea_250")
  #---repicating to match dim of SL_df---
  GIA_slope <- rep(NA, nrow(SL_df))
  for(i in 1:nrow(SL_df)){
    GIA_slope[i] <- (SL[lon_index[i],lat_index[i]])
  }
  
  #---Combining GIA with other dataset--
  SL_df <- cbind(SL_df,ICE5_GIA_slope = GIA_slope)#mm/yr
  
  write_csv(SL_df,file = save_loc)
  return(SL_df)
  cat("Add GIA rate to dataset")
}