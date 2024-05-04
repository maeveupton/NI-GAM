clean_tidal_gauge_data <- function(SL_df,save_loc,path_to_data){
  #---Annual Tidal Gauge data----
  annual_tidal_gauge_data <- read_csv(path_to_data)
  annual_tidal_gauge_data_df <- annual_tidal_gauge_data %>% 
    dplyr::select(Age,RSL, Latitude,Longitude,name,RSL_offset,Age_epoch_id) %>%
    rename(SiteName = name) %>% 
    mutate(RSL = RSL/1000)  # from mm --> m
  
  #---Reordering by group---
  annual_tidal_gauge_tidy_df <- annual_tidal_gauge_data_df %>% 
    group_by(SiteName) %>% 
    arrange(SiteName, .by_group = TRUE)
  #----Decadal Averages------
  decadal_averages_TG <-
    annual_tidal_gauge_tidy_df %>%
    mutate(decade = (Age - 1) %/% 10) %>%
    group_by(decade,SiteName) %>%
    summarise(decade_meanRSL = mean(RSL),
              #sd_TG = sd(RSL),
              Age=max(Age),
              rows_site = n())#Age=min(Age)
  #---Using standard deviation of RSL as uncertainty----
  decadal_averages_TG <- decadal_averages_TG %>%
    group_by(SiteName) %>% 
    mutate(sd_TG = sd(decade_meanRSL))
  
  #----- New df with decadal averages for tide gauges-----
  tidal_gauge_average_10_df <- merge(decadal_averages_TG,annual_tidal_gauge_tidy_df)
  
  #---Rsl & Age error for tidal gauge data----
  tidal_gauge_full_df <-tidal_gauge_average_10_df %>% 
    mutate(Age_er_min = 5,#years --> half a year/half a decade
           Age_er_max = 5) %>%
    mutate(sd_TG = ifelse(is.na(sd_TG),0.001,sd_TG)) %>% 
    group_by(SiteName) %>%
    mutate(RSL_er_min =   sd_TG,
           RSL_er_max =  sd_TG)
  
  tidal_gauge_full_df <- tidal_gauge_full_df %>%
    mutate(Age = Age/1000) %>% 
    mutate(Age_er_max = Age_er_max/1000) %>% 
    mutate(Age_er_min = Age_er_min/1000) %>% 
    mutate(Age_er_average = (Age_er_max+Age_er_min)/2) %>% 
    mutate(RSL_er_average = (RSL_er_min + RSL_er_max)/2) %>% 
    mutate(RSL_annual = RSL) %>% 
    mutate(RSL = decade_meanRSL)
  
  #----Range of Longitude and Lat from Data for East Coast of North America----
  Long_range <- SL_df$Longitude %>% range()
  Lat_range <- SL_df$Latitude %>% range()
  
  #---Removing the global and focusing on east coast of US
  decadal_NA_TG<-tidal_gauge_full_df %>%
    filter(between(Longitude,min(Long_range)-5,max(Long_range)+5) & between(Latitude,min(Lat_range)-5,max(Lat_range)+5))
  
  #--Removing sites with only 2 points (20 years of data)---
  decadal_NA_TG_df <- decadal_NA_TG %>%
    group_by(SiteName) %>%
    #filter(n() > 3) %>% 
    filter(n() >= 2) %>% 
    dplyr::select(!decade,decade_meanRSL,RSL_annual)
  
  #-----Uniting original dataset and model run to give a site index to model_result data set-----
  SL_site_df<- SL_df %>% 
    mutate(Longitude = round(Longitude,1)) %>% 
    mutate(Latitude = round(Latitude,1)) %>% 
    unite("LongLat",Latitude:Longitude,remove = FALSE) %>%  #Uniting 2 columns
    mutate(site = sprintf("%02d", as.integer(as.factor(LongLat)))) %>% 
    mutate(data_type_id = "ProxyData") %>% 
    group_by(SiteName) %>% 
    mutate(Longitude = first(Longitude),
           Latitude = first(Latitude))
  
  SL_tide_site_df<- decadal_NA_TG_df %>% 
    #dplyr::select(!all_tidal_data_sites) %>% 
    mutate(Longitude = round(Longitude,1)) %>% 
    mutate(Latitude = round(Latitude,1)) %>% 
    unite("LongLat",Latitude:Longitude,remove=FALSE) %>%  #Uniting 2 columns
    mutate(site = sprintf("%02d", as.integer(as.factor(LongLat)))) %>% 
    mutate(data_type_id = "TideGaugeData")
  
  #------Joining proxy dataframe to Tide gauges data----
  SL_tide_proxy <- bind_rows(SL_site_df,SL_tide_site_df)
  
  #------Joining proxy sites to gauges based on shortest distance----
  SL_proxy_unique <- SL_site_df %>%
    dplyr::select(SiteName, Longitude, Latitude,data_type_id) %>%
    unique()
  SL_tide_unique <-  SL_tide_site_df %>% 
    dplyr::select(SiteName, Longitude, Latitude,data_type_id) %>%
    unique()
  
  #---Distance Matrix for each site to each other---
  mat.distance<- distm(SL_proxy_unique[,2:3],SL_tide_unique[,2:3],
                       fun = distGeo)
  mat.distance <-  as.data.frame(mat.distance)
  #--finding row mins & corresponding tidal gauge--
  rownames(mat.distance) = SL_proxy_unique$SiteName
  colnames(mat.distance) = SL_tide_unique$SiteName
  #--finding row mins & corresponding tidal gauge--
  dist_TG_proxy <- t(sapply(seq(nrow(mat.distance)), function(z) {
    js <- order(mat.distance[z,])[1:4]
    c(rownames(mat.distance)[z], colnames(mat.distance)[js[1]], mat.distance[z,js[1]],
      colnames(mat.distance)[js[2]], mat.distance[z,js[2]],
      colnames(mat.distance)[js[3]], mat.distance[z,js[3]],
      colnames(mat.distance)[js[4]], mat.distance[z,js[4]])
  }))
  
  dist_TG_proxy <- as.data.frame(dist_TG_proxy)
  colnames(dist_TG_proxy) <- c("nearest_proxy_site", "SiteName1", "min_dist1",
                               "SiteName_duplicate2", "min_dist2","SiteName_duplicate3",
                               "min_dist3","SiteName_duplicate4","min_dist4")
  # Sorting the minimum distances from lowest to highest
  dist_TG_proxy <- dist_TG_proxy %>% arrange(desc(min_dist1))
  # TG near the proxy sites & TG longer than 150 years (New York(The Battery))
  all_nearest_TG <- dist_TG_proxy %>% 
    dplyr::select(!c(nearest_proxy_site)) %>% 
    pivot_longer(cols = starts_with("SiteName"),
                 values_to = "SiteName") %>% 
    dplyr::select(!name) %>% 
    pivot_longer(cols = starts_with("min_dist"),
                 values_to = "MinimumDistance")
  
  # 1 degree away from proxy is 111.1km
  all_nearest_TG_closest <- all_nearest_TG %>% filter(MinimumDistance>111100)
  
  # Joining the selected TG sites back with the original data
  join_new_index_tide_df <- SL_tide_site_df %>%
    filter(SiteName %in% all_nearest_TG_closest$SiteName)
  
  
  #--There will be NAs were the proxy data doesn't have a corresponding index--
  SL_df_tide_proxy <- plyr::rbind.fill(SL_site_df,
                                       join_new_index_tide_df)#stacking rows
  
  
  # Ensuring the SiteName is a factor
  SL_df <- SL_df_tide_proxy %>% 
    dplyr::select(!c(RSL_annual,Age_epoch_id,
                     RSL_offset,sd_TG,rows_site,decade_meanRSL,
                     Indicator,Basin,Age_2_er_max,Age_2_er_min))%>% 
    mutate(SiteName = as.factor(SiteName))
  write_csv(SL_df,save_loc)
  cat("Tidal Gauge Data added to overall dataset")
  return(SL_df)
}
