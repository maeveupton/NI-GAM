clean_tidal_gauge_data <- function(SL_df,save_loc,path_to_data,andy_TG_site){
  
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
  
  SL_tide_site_df<- decadal_NA_TG %>% 
    #dplyr::select(!all_tidal_data_sites) %>% 
    mutate(Longitude = round(Longitude,1)) %>% 
    mutate(Latitude = round(Latitude,1)) %>% 
    unite("LongLat",Latitude:Longitude,remove=FALSE) %>%  #Uniting 2 columns
    mutate(site = sprintf("%02d", as.integer(as.factor(LongLat)))) %>% 
    mutate(data_type_id = "TideGaugeData")
  
  #------Joining proxy dataframe to Tide gauges data----
  SL_tide_proxy_df <- bind_rows(SL_site_df,SL_tide_site_df)
  
  #--Removing sites with only 3 points (30 years of data)---
  SL_tide_proxy_df <- SL_tide_proxy_df %>%
    group_by(SiteName) %>%
    #filter(n() > 3) %>% 
    filter(n() >= 2) %>% 
    dplyr::select(!decade,decade_meanRSL,RSL_annual)
  
  
  #------Joining proxy sites to gauges based on shortest distance----
  SL_proxy_unique <- SL_site_df %>%
    dplyr::select(SiteName, Longitude, Latitude,data_type_id) %>%
    unique()
  SL_tide_unique <-  SL_tide_proxy_df %>% #SL_tide_site_df %>%
    dplyr::select(SiteName, Longitude, Latitude,data_type_id) %>%
    unique()
  
  
  # If you want more tide gauges
  all_nearest_TG <- dist_TG_proxy %>% 
    dplyr::select(!c(nearest_proxy_site,dup_check,dup_check2)) %>% 
    pivot_longer(cols = starts_with("SiteName"),
                  values_to = "SiteName") %>% 
    dplyr::select(!name) %>% 
    pivot_longer(cols = starts_with("min_dist"),
                 values_to = "MinimumDistance")
  # Include Andy's recommended TGs
  all_nearest_TG_andy <- union(unique(all_nearest_TG$SiteName),andy_TG_site_list$AndyTG_sites)
  
  # Joining the selected TG sites back with the original data
  join_new_index_tide_df <- SL_tide_site_df %>% filter(SiteName %in% all_nearest_TG_andy)

  
  #--There will be NAs were the proxy data doesn't have a corresponding index--
  SL_df_tide_proxy <- plyr::rbind.fill(SL_site_df,
                                       join_new_index_tide_df)#stacking rows
  
  
  # #--Replace NA with site names--
  # SL_df_tide_proxy_tidy<- SL_df_tide_proxy %>%
  #   mutate(nearest_proxy_site = coalesce(nearest_proxy_site,SiteName)) %>%
  #   rename(TideGaugeSite = SiteName) %>%
  #   rename(SiteName=nearest_proxy_site)#duplicate columns with same name
  # 
  # #--Changing long & lat values to be the same for sites --
  # overall_SL_tide_df <- SL_df_tide_proxy_tidy %>%
  #   group_by(SiteName) %>%
  #   mutate(Longitude = first(Longitude),
  #          Latitude = first(Latitude)) #assign the first code value for each group.
  
  #--Writing new dataframe----
  SL_df <- SL_df_tide_proxy %>%  
    dplyr::select(-c(decade,decade_meanRSL))
  
  # #---Removing Tide gauge data from certain sites-----
  # #---Find the min age in the tide gauge record to see when to join tide with proxy---
  # min_Age_tide <- SL_df %>% 
  #   filter(data_type_id == "TideGaugeData") %>% 
  #   group_by(SiteName) %>% 
  #   summarise(min_TG_Age = min(Age))
  # SL_df <- left_join(SL_df,min_Age_tide, by = "SiteName")
  # 
  # # Find max proxy age to see how far TG is away from Proxy in time
  # max_Age_proxy <- SL_df %>% 
  #   filter(data_type_id == "ProxyData") %>% 
  #   group_by(SiteName) %>% 
  #   summarise(max_Age_proxy = max(Age))
  # 
  # SL_df <- left_join(SL_df,max_Age_proxy, by = "SiteName")
  # 
  # # How far is the TG from the proxy in terms of Age?
  # SL_df <- SL_df %>% 
  #   group_by(SiteName) %>% 
  #   mutate(Age_diff = abs(min_TG_Age - max_Age_proxy))
  # 
  # # Removing sites with large temporal distances between TG and Proxy
  # # 150 years separation between proxy & TG
  # SL_df$No_TG <- ifelse(SL_df$Age_diff >= 0.15,TRUE,FALSE)
  # #---Remove TG data for proxy with a gap between TG & proxy greater than 100 years---
  # remove_TG <- SL_df %>%  filter(No_TG=="TRUE" & data_type_id == "TideGaugeData")
  # SL_df <- anti_join(SL_df,remove_TG)
  # Ensuring the SiteName is a factor
  SL_df <- SL_df %>% mutate(SiteName = as.factor(SiteName))
  write_csv(SL_df,save_loc)
  cat("Tidal Gauge Data added to overall dataset")
  return(SL_df)
}
