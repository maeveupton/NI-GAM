add_GIA_rate<- function(SL_df,save_loc){
  #--------------------GIA DATA from Peltier Website ICE6G--------------
  download.file('https://www.dropbox.com/s/wiy8hs74as6fv9l/dsea.1grid_O512.nc?dl=1',
                destfile = "dsea.1grid_O512.nc",
                method = "libcurl", mode = "wb")
  nc_data <- nc_open('dsea.1grid_O512.nc')#ICE6G -- gives weird result for sites but Tanghua uses it for GIA uncertainty
  
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
  SL_df <- cbind(SL_df,ICE6_GIA_slope = GIA_slope)#mm/yr
  
  #SL_df <- cbind(SL_df,GIA_slope_mmyr = GIA_slope)
  #SL_df <- SL_df %>% mutate(GIA_slope = GIA_slope_mmyr/1000)  # GIA slope in m/yr
  
  # #---Putting RSL into mm-----
  # SL_df <- SL_df %>% mutate(RSL = RSL*1000) %>% 
  #   mutate(RSL_er_max = RSL_er_max*1000) %>% 
  #   mutate(RSL_er_min = RSL_er_min*1000) %>% 
  #   mutate(RSL_er_average = RSL_er_average*1000)
  
  
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
  
  
  #---Including GIA rate & uncertainties from Tanghua paper--
  GIA_rates_Tanghua <- read.table(file = "https://www.dropbox.com/s/izw9xaetynmpfj8/rslt_Mean_Global.txt?dl=1", header = TRUE)
  colnames(GIA_rates_Tanghua) = c("Longitude","Latitude","GIA_rate")
  GIA_rate_uncertainty_Tanghua<- read.table(file = "https://www.dropbox.com/s/m4gmlxp7p86y31a/rslt_2Sigma_Global.txt?dl=1", header = TRUE)
  colnames(GIA_rate_uncertainty_Tanghua) = c("Longitude","Latitude","GIA_rate_2_uncertainty")
  
  GIA_rate_Tanghua_df <- cbind(GIA_rates_Tanghua,
                                    GIA_rate_2_uncertainty=GIA_rate_uncertainty_Tanghua$GIA_rate_2_uncertainty)
  
  GIA_rate_Tanghua_df <- GIA_rate_Tanghua_df %>% 
    mutate(GIA_rate_uncertainty = GIA_rate_2_uncertainty/2) %>% 
    mutate(Longitude = round(Longitude,1)) %>% 
    mutate(Latitude = round(Latitude,1)) %>% 
    unite("LongLat",Latitude:Longitude,remove=FALSE) %>%  #Uniting 2 columns
    mutate(site = sprintf("%02d", as.integer(as.factor(LongLat))))
  
  #---Matching the long lat from Tanghua to the SL data---
  # SL_df <- left_join(SL_df,GIA_rate_Tanghua_df,by = c("LongLat" = "LongLat",
  #                                                     "Longitude"="Longitude",
  #                                                     "Latitude"="Latitude"))
  
  
  gia_lat_new <- tibble(index = 1:length(GIA_rate_Tanghua_df$Latitude),site = GIA_rate_Tanghua_df$site,
                        Latitude = GIA_rate_Tanghua_df$Latitude) %>% arrange(Latitude)
  gia_lon_new <- tibble(index = 1:length(GIA_rate_Tanghua_df$Longitude),site = GIA_rate_Tanghua_df$site,
                        Longitude=GIA_rate_Tanghua_df$Longitude) %>% arrange(Longitude)
  
  #--- Matching closest long & lat values---
  lat_index_new <- gia_lat_new$index[match.closest(dat_lat, gia_lat_new$Latitude)]
  lon_index_new <- gia_lon_new$index[match.closest(dat_lon, gia_lon_new$Longitude)] 
  long_lat_tang<-data.frame(Latitude = GIA_rate_Tanghua_df$Latitude[lat_index_new], 
                          Longitude = GIA_rate_Tanghua_df$Longitude[lon_index_new])
  new_tang_df <- long_lat_tang %>%  
    unite("LongLat",Latitude:Longitude,remove=FALSE)   #Uniting 2 columns

  new_tang_df <- left_join(new_tang_df,GIA_rate_Tanghua_df,by = c("LongLat","Longitude","Latitude"))
  
  #---Combining GIA with other dataset--
  SL_df <- cbind(SL_df,
                 GIA_slope = new_tang_df$GIA_rate,
                 sd_GIA_slope = new_tang_df$GIA_rate_uncertainty)#mm/yr
  
  # #---------------Detrending the data -----------
  # SL_df <- SL_df %>% group_by(SiteName) %>%
  #   #--Using Kopps reference point 2000---
  #   mutate(SL = (2000-Age)*GIA_slope + RSL) # m/yr
  
  #---------------Detrending the data -----------
  SL_df <- SL_df %>% group_by(SiteName) %>%
    #--Using Kopps reference point 2000---
    mutate(SL = (2 - Age)*GIA_slope + RSL) # mm/yr
  
  #---Plot detrended data---
  detrend_plot <- ggplot()+
    geom_point(data = SL_df, aes(x = Age, y = RSL),colour = "red")+
    geom_point(data = SL_df, aes(x = Age, y = SL ),colour= "blue")+
    facet_wrap(~SiteName)+
    geom_text(data = SL_df,mapping = aes(x = 1, y = -3, label = GIA_slope))+
    theme_bw()+
    ggtitle("Detrended using Peltier's Dataset")
  detrend_plot
  ggsave(detrend_plot, filename = "fig/detrended_data.pdf",width = 10, height = 6)
  # #-----t - t0--
  # t_ref <- SL_df %>%
  #   dplyr::select(Age,SiteName) %>% group_by(SiteName) %>%
  #   #summarise(Age = Age,Age_ref_point=Age - 2000)
  #   summarise(Age = Age,Age_ref_point=Age - 2)
  # 
  # SL_df$Age_ref_point <- t_ref$Age_ref_point
  
  write_csv(SL_df,file = save_loc)
  return(SL_df)
  cat("Add GIA rate to dataset")
}