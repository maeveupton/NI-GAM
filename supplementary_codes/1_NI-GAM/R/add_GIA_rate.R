add_GIA_rate<- function(SL_df,save_loc){
  #-----GIA from Caron 2018------
  caron_GIA <- read_csv("https://www.dropbox.com/s/zpfqa73jthht3qr/GIA_Caron2018.csv?dl=1")
  
  caron_GIA_df <- caron_GIA %>% 
    mutate(Latitude = as.numeric(colatitude) - 90) %>% 
    mutate(Longitude = ifelse(longitude >=180,longitude -360,longitude)) %>% 
    filter(between(Longitude,min(SL_df$Longitude)-5,max(SL_df$Longitude)+5) & between(Latitude,min(SL_df$Latitude)-5,max(SL_df$Latitude)+5)) %>% 
    unite("LongLat",Latitude:Longitude,remove=FALSE) %>%  #Uniting 2 columns to make column ID
    mutate(grid_site = sprintf("%02d", as.integer(as.factor(LongLat)))) 
    
  #plot(caron_GIA_df$Longitude,caron_GIA_df$Latitude)
  #points(SL_df$Longitude,SL_df$Latitude,col="red")
  
  # 87 sites- need 87 corresponding GIA rates
  SL_df_lon_lat_unique <- SL_df %>% dplyr::select(Longitude,Latitude,SiteName)%>%
    unique() #%>% dplyr::select(!SiteName)
  caron_GIA_df_unique <- caron_GIA_df %>% dplyr::select(grid_site,Longitude,Latitude) %>% unique
  
  #---Distance Matrix for GIA grid to dataset---
  mat.distance<- distm(SL_df_lon_lat_unique[,1:2],caron_GIA_df_unique[,2:3],
                       fun = distGeo)
  mat.distance <-  as.data.frame(mat.distance)
  #--finding row mins & corresponding tidal gauge--
  rownames(mat.distance) = SL_df_lon_lat_unique$SiteName#as.character(SL_df_lon_lat$SiteName %>% unique)
  colnames(mat.distance) = caron_GIA_df_unique$grid_site
  #--finding row mins --
  dist_data_grid <- t(sapply(seq(nrow(mat.distance)), function(z) {
    js <- order(mat.distance[z,])[1:2]
    c(rownames(mat.distance)[z], colnames(mat.distance)[js[1]], mat.distance[z,js[1]])
  }))
  
  dist_data_grid_df <- as.data.frame(dist_data_grid)
  colnames(dist_data_grid_df) <- c("nearest_site", "grid_site", "min_dist1")
  caron_GIA_df_filter <- caron_GIA_df %>%  filter(grid_site %in%dist_data_grid_df$grid_site)
  #----GIA rate for each site----
  dist_data_grid_df_GIA <- dist_data_grid_df %>% 
    dplyr::select(!min_dist1) %>% 
    left_join(.,caron_GIA_df_filter,by = "grid_site")# VLM is the GIA rate
  SL_data_GIA_grid <- dist_data_grid_df_GIA %>% 
    rename(SiteName = nearest_site,GIA_caron = VLM,GIA_caron_sd = VLM_sd) %>% 
    dplyr::select(SiteName,GIA_caron,GIA_caron_sd)
  
  #----Joining Caron GIA rates back with big SL_df
  SL_df <- SL_df %>% left_join(.,SL_data_GIA_grid,by = "SiteName")

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
  
  write_csv(SL_df,file = save_loc)
  return(SL_df)
  cat("Add GIA rate to dataset")
}
  #---Including GIA rate & uncertainties from Tanghua paper--
  # GIA_rates_Tanghua <- read.table(file = 
  #                                   "https://www.dropbox.com/s/izw9xaetynmpfj8/rslt_Mean_Global.txt?dl=1", header = TRUE)
  # colnames(GIA_rates_Tanghua) = c("Longitude","Latitude","GIA_rate")
  # GIA_rate_uncertainty_Tanghua<- read.table(file = 
  #                                             "https://www.dropbox.com/s/m4gmlxp7p86y31a/rslt_2Sigma_Global.txt?dl=1", header = TRUE)
  # colnames(GIA_rate_uncertainty_Tanghua) = c("Longitude","Latitude","GIA_rate_2_uncertainty")
  # 
  # GIA_rate_Tanghua_df <- cbind(GIA_rates_Tanghua,
  #                                   GIA_rate_2_uncertainty=GIA_rate_uncertainty_Tanghua$GIA_rate_2_uncertainty)
  # 
  # GIA_rate_Tanghua_df <- GIA_rate_Tanghua_df %>% 
  #   mutate(GIA_rate_uncertainty = GIA_rate_2_uncertainty/2) %>% 
  #   mutate(Longitude = round(Longitude,1)) %>% 
  #   mutate(Latitude = round(Latitude,1)) %>% 
  #   unite("LongLat",Latitude:Longitude,remove=FALSE) %>%  #Uniting 2 columns
  #   mutate(site = sprintf("%02d", as.integer(as.factor(LongLat))))
  # 
  # #---Matching the long lat from Tanghua to the SL data---
  # # SL_df <- left_join(SL_df,GIA_rate_Tanghua_df,by = c("LongLat" = "LongLat",
  # #                                                     "Longitude"="Longitude",
  # #                                                     "Latitude"="Latitude"))
  # 
  # 
  # gia_lat_new <- tibble(index = 1:length(GIA_rate_Tanghua_df$Latitude),site = GIA_rate_Tanghua_df$site,
  #                       Latitude = GIA_rate_Tanghua_df$Latitude) %>% arrange(Latitude)
  # gia_lon_new <- tibble(index = 1:length(GIA_rate_Tanghua_df$Longitude),site = GIA_rate_Tanghua_df$site,
  #                       Longitude=GIA_rate_Tanghua_df$Longitude) %>% arrange(Longitude)
  # 
  # #--- Matching closest long & lat values---
  # lat_index_new <- gia_lat_new$index[match.closest(dat_lat, gia_lat_new$Latitude)]
  # lon_index_new <- gia_lon_new$index[match.closest(dat_lon, gia_lon_new$Longitude)] 
  # long_lat_tang<-data.frame(Latitude = GIA_rate_Tanghua_df$Latitude[lat_index_new], 
  #                         Longitude = GIA_rate_Tanghua_df$Longitude[lon_index_new])
  # new_tang_df <- long_lat_tang %>%  
  #   unite("LongLat",Latitude:Longitude,remove=FALSE)   #Uniting 2 columns
  # 
  # new_tang_df <- left_join(new_tang_df,GIA_rate_Tanghua_df,by = c("LongLat","Longitude","Latitude"))
  # 
  # #---Combining GIA with other dataset--
  # SL_df <- cbind(SL_df,
  #                GIA_slope = new_tang_df$GIA_rate,
  #                sd_GIA_slope = new_tang_df$GIA_rate_uncertainty)#mm/yr
  
  # #---------------Detrending the data -----------
  # SL_df <- SL_df %>% group_by(SiteName) %>%
  #   #--Using Kopps reference point 2000---
  #   mutate(SL = (2000-Age)*GIA_slope + RSL) # m/yr
  
  # #---------------Detrending the data -----------
  # SL_df <- SL_df %>% group_by(SiteName) %>%
  #   #--Using Kopps reference point 2000---
  #   mutate(SL = (2 - Age)*GIA_slope + RSL) # mm/yr
  # 
  # #---Plot detrended data---
  # detrend_plot <- ggplot()+
  #   geom_point(data = SL_df, aes(x = Age, y = RSL),colour = "red")+
  #   geom_point(data = SL_df, aes(x = Age, y = SL ),colour= "blue")+
  #   facet_wrap(~SiteName)+
  #   geom_text(data = SL_df,mapping = aes(x = 1, y = -3, label = GIA_slope))+
  #   theme_bw()+
  #   ggtitle("Detrended using Peltier's Dataset")
  # detrend_plot
  #ggsave(detrend_plot, filename = "fig/detrended_data.pdf",width = 10, height = 6)
