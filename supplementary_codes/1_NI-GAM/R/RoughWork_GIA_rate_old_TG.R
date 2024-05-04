# #----GIA boxplot------
# data_GIA_slope = SL_df %>% 
#   group_by(SiteName) %>% 
#   dplyr::select(data_lm_slope,data_lm_slope_err) %>%
#   unique() %>% arrange(SiteName)  
# 
# post_GIA<- model_run$BUGSoutput$sims.list$b_GIA
# post_GIA_df <- data.frame(post_GIA)
# colnames(post_GIA_df) <- data_GIA_slope$SiteName
# post_GIA_df_long <- post_GIA_df %>% 
#   pivot_longer(cols = everything())
# colnames(post_GIA_df_long) <- c("SiteName","GIA_rate_post")

# ICE6_G_rates <- SL_df %>% 
#   group_by(SiteName) %>% 
#   dplyr::select(ICE6_GIA_slope) %>%
#   unique() %>% arrange(SiteName) 
# ICE5_G_rates <- SL_df %>% 
#   group_by(SiteName) %>% 
#   dplyr::select(ICE5_GIA_slope) %>%
#   unique() %>% arrange(SiteName)
# 
# Output results in a table
#GIA_table <- data.frame(data_GIA_slope,)

# GIA_boxplot <- ggplot(data = post_GIA_df_long, aes(x = SiteName, y = GIA_rate_post)) +
#   geom_boxplot() +
#   theme_bw()+
#   geom_point(data = data_GIA_slope, 
#               aes(x = SiteName,y = data_lm_slope,colour = "Data GIA Slopes"),shape=18,size=3)+
#   
#   geom_point(data = ICE5_G_rates, 
#               aes(x = SiteName,y = ICE5_GIA_slope, colour="ICE5G GIA Rate"),shape=18,size=3)+
#   # # geom_point(data = ICE6_G_rates, 
#   #            aes(x = SiteName,y = ICE6_GIA_slope,colour=" ICE6G GIA Rate"),shape=18,size=3)+
#   xlab("") +
#   ylab("GIA Rate (mm/yr)")+
#   #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7))+
#   ggtitle("GIA rate estimates")+
#   guides(colour=guide_legend(title=""))+
#   coord_flip()
# GIA_boxplot
# ggsave(paste0("fig/",save_option,"/", "GIA_boxplot.pdf", sep = ""),
#        GIA_boxplot, width = 10, height = 6)

#ggsave(GIA_boxplot,file ="fig/gia_boxplot_update.pdf", width = 10, height = 6)




#----- Plotting all posterior for regional component-----
time_component <- model_run$BUGSoutput$sims.list$regional
post_regional <- as.data.frame(time_component[1:10,])
rownames(post_regional) <- paste0("sample_",1:nrow(post_regional))
post_regional_long_df <- pivot_longer(data = as.data.frame(cbind(Age = SL_df$Age,t(post_regional))),
                                      cols = starts_with("sample_"))

# regional_post_plot <- ggplot(data  = post_regional_long_df,aes(Age,value))+
#   geom_point()
# ggsave(paste0("fig/",save_option,'/',"regional_post_plot.pdf", sep = ""),
#        regional_post_plot, width = 10, height = 6)


# #---TG don't have GIA values from LR use closest proxy site----
# slopes_proxy <- SL_df %>% group_by(SiteName) %>%
#   filter(data_type_id == "ProxyData") %>% 
#   summarise(slopes = unique(data_lm_slope), error = unique(data_lm_slope_err))
# 
# # Joining proxy sites to gauges based on shortest distance
# SL_proxy_unique <- SL_df %>%
#   filter(data_type_id == "ProxyData") %>% 
#   dplyr::select(SiteName, Longitude, Latitude) %>%
#   unique()
# SL_tide_unique <-  SL_df %>% 
#   filter(data_type_id == "TideGaugeData") %>% 
#   dplyr::select(SiteName, Longitude, Latitude) %>%
#   unique()
# 
# #---Distance Matrix for each site to each other---
# mat.distance<- distm(SL_tide_unique[,2:3],SL_proxy_unique[,2:3],
#                      fun = distGeo)
# mat.distance <-  as.data.frame(mat.distance)
# #--finding row mins & corresponding tidal gauge--
# rownames(mat.distance) = SL_tide_unique$SiteName
# colnames(mat.distance) = SL_proxy_unique$SiteName
# #--finding row mins & corresponding tidal gauge--
# dist_TG_proxy <- t(sapply(seq(nrow(mat.distance)), function(z) {
#   js <- order(mat.distance[z,])[1:1]
#   c(rownames(mat.distance)[z], colnames(mat.distance)[js[1]], mat.distance[z,js[1]])
# }))
# dist_TG_proxy <- as.data.frame(dist_TG_proxy)
# colnames(dist_TG_proxy) <- c( "tide_gauge_site","nearest_proxy_site", "min_dist")
# # Sorting the minimum distances from lowest to highest
# dist_TG_proxy_df <- dist_TG_proxy %>% 
#   arrange(desc(min_dist)) %>% 
#   dplyr::rename(SiteName = tide_gauge_site) %>% 
#   dplyr::select(-min_dist)
# 
# dist_TG_proxy_GIA_df <-dist_TG_proxy_df %>% 
#   mutate(GIA_TG_new = ifelse(match(nearest_proxy_site,slopes_proxy$SiteName),slopes_proxy$slopes,0)) %>% 
#   mutate(GIA_err_TG_new = ifelse(match(nearest_proxy_site,slopes_proxy$SiteName),slopes_proxy$error,0)) %>% 
#   dplyr::select(-nearest_proxy_site)
# # Joining the GIA rates from proxy to TG sites
# SL_df<-left_join(SL_df,dist_TG_proxy_GIA_df, by="SiteName")
# SL_df <- SL_df %>% 
#   mutate(data_lm_slope = ifelse(data_type_id=="TideGaugeData",GIA_TG_new,data_lm_slope),
#          data_lm_slope_err = ifelse(data_type_id=="TideGaugeData",GIA_err_TG_new,data_lm_slope))
# SL_df <- SL_df %>% mutate(SiteName = as.factor(SiteName))



clean_tidal_gauge_data <- function(SL_df,save_loc,path_to_data,andy_TG_site){
  path_to_data = "https://raw.githubusercontent.com/maeveupton/tide_gauge_data_PSMSL/main/annual_SL_tide_df.csv"
  andy_TG_site = "https://www.dropbox.com/s/jjcimcc6q5kw8sj/TG_proxy_sites_new.csv?dl=1"
  
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
  #---Standard deviation of the data----
  sd_TG_data <- annual_tidal_gauge_tidy_df %>% 
    group_by(SiteName) %>% 
    summarise(sd_TG = sd(RSL)) 
  #----Decadal Averages------
  decadal_averages_TG <-
    annual_tidal_gauge_tidy_df %>%
    mutate(decade = (Age - 1) %/% 10) %>%
    group_by(decade,SiteName) %>%
    summarise(decade_meanRSL = mean(RSL),Age=min(Age))
  
  #----- New df with decadal averages for tide gauges-----
  tidal_gauge_full_df <- merge(decadal_averages_TG,annual_tidal_gauge_tidy_df)
  
  #---Rsl & Age error for tidal gauge data----
  tidal_gauge_full_df$Age_er_min <- 5#0.5#0.0001 #years --> half a year/half a decade
  tidal_gauge_full_df$Age_er_max <- 5#0.50.0001 # years
  tidal_gauge_full_df$RSL_er_min <- 0.02# 2 cm#0.005 #mm- m
  tidal_gauge_full_df$RSL_er_max <- 0.02#2 cm#0.005 #mm - m
  
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
    filter(n() > 3) %>% 
    dplyr::select(!decade,decade_meanRSL,RSL_annual)
  
  
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
  # Finding the duplicates
  dist_TG_proxy$dup_check <- as.factor(duplicated(dist_TG_proxy$SiteName1))
  # Replacing Duplicates
  dist_TG_proxy <- dist_TG_proxy %>% 
    mutate(SiteName = ifelse(dup_check==TRUE,
                             as.character(SiteName_duplicate2),
                             as.character(SiteName1)))
  # More duplicates created
  dist_TG_proxy$dup_check2 <- as.factor(duplicated(dist_TG_proxy$SiteName))
  # Replacing new duplicates
  dist_TG_proxy <- dist_TG_proxy %>% 
    mutate(SiteName = ifelse(dup_check2==TRUE,
                             as.character(SiteName_duplicate2),
                             as.character(SiteName)))
  
  # Using only the nearest_proxy_site & SiteName columns to join with rest of data
  dist_TG_proxy_unique <- dist_TG_proxy %>% 
    dplyr::select(nearest_proxy_site,SiteName)
  # Check if the nearest sites is best TG was selected based on Andy's recommendation list
  andy_TG_site_list <- read_csv(andy_TG_site) %>% 
    dplyr::select(SiteName, AndyTG_sites) %>% 
    rename(nearest_proxy_site = SiteName) %>% 
    mutate(nearest_proxy_site = as.factor(nearest_proxy_site))
  dist_TG_proxy_unique_andy <- left_join(dist_TG_proxy_unique,andy_TG_site_list, by ="nearest_proxy_site")
  # If you want sites near tide gauges
  # # Make the column of Andy's TG sites to be the ones we want --> Will need to change this in the future
  # dist_TG_proxy_df <- dist_TG_proxy_unique_andy %>% 
  #   mutate(SiteName = coalesce(SiteName,AndyTG_sites))
  #   #dplyr::select(nearest_proxy_site,AndyTG_sites) %>% 
  #   #rename(SiteName = AndyTG_sites)
  
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
  #join_new_index_tide_df <- left_join(dist_TG_proxy_df,SL_tide_site_df,by = "SiteName")
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
