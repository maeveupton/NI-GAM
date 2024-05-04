clean_tidal_gauge_data <- function(SL_df,save_loc,sediment_average_TG = 10){
  # Using data from PSMSL website for annual tide gauge data----------------------------------
  # Set up the URL for downloading the data
  url <- "https://psmsl.org/data/obtaining/rlr.annual.data/rlr_annual.zip"
  
  # Create a temporary file
  temp_file <- tempfile()
  
  # Download the file and save it to the temporary file
  utils::download.file(url,
                       destfile = temp_file,
                       quiet = TRUE)
  
  # Back up if an error
  # # Create a temporary file
  # temp_file <- paste0(tempdir(),'/tg',format(Sys.Date(), "%Y%m%d"))
  # # Download the file and save it to the temporary file
  # if(file.exists(temp_file) | override_download) {
  #   utils::download.file(url,
  #                        destfile = temp_file,
  #                        quiet = TRUE)
  # }
  
  # Unzip the data file to a temporary directory
  temp_dir <- tempfile()
  utils::unzip(temp_file, exdir = temp_dir)
  
  
  ### ------------Loop to open all RSL & Age data files------------
  read_plus <- function(flnm) {
    # fread quicker way to read in & allows for ; to be used
    data.table::fread(flnm, sep = ";") %>%
      # allows you to include the file name as id
      dplyr::mutate(filename = flnm)
  }
  # Warnings: there are some files without data
  suppressWarnings(
    temp_SL <-
      list.files(
        path = file.path(temp_dir, "rlr_annual", "data"),
        pattern = "*.rlrdata",
        full.names = T
      ) %>%
      purrr::map_df(~ read_plus(.)) %>%
      dplyr::tibble()
  )
  
  colnames(temp_SL) <- c("Age", "RSL", "flag_attention_1", "flag_attention_2", "id")
  temp_SL$id <- stringr::str_extract(basename(temp_SL$id), "[0-9]+")
  
  # Access the individual data files within the 'rlr_annual' folder
  file_path <- file.path(temp_dir, "rlr_annual", "filelist.txt")
  file_list <- utils::read.csv(file_path, stringsAsFactors = FALSE, header = F, sep = ";")
  colnames(file_list) <- c(
    "id", "Latitude", "Longitude", "name",
    "coastline", "stationcode", "stationflag"
  )
  
  # Removing white space in the name of each site
  file_list$name <- stringr::str_trim(file_list$name, side = "both")
  file_list$stationflag <- stringr::str_trim(file_list$stationflag, side = "both")
  
  # Data from the PSMSL website
  data_TG <- temp_SL %>%
    # Pulling out the file number from string so that it matches the name from other files
    dplyr::mutate(id = stringr::str_extract(basename(temp_SL$id), "[0-9]+")) %>%
    # Cases where bad data was collected
    dplyr::filter(!RSL == -99999) # %>%
  # dplyr::group_by(id) %>%
  # 2000-2018 used as the tidal epoch
  # dplyr::mutate(Age_epoch_id = ifelse(dplyr::between(Age, 2000, 2018), TRUE, FALSE))
  
  # Removing offset based on the location---
  # Offset value is the mean of RSL over the tidal epoch
  # Setting 2000-2018 as the tidal epoch
  # Age_epoch_ref <- data_TG %>%
  #  dplyr::select(RSL, Age_epoch_id) %>%
  #  dplyr::filter(Age_epoch_id == TRUE) %>%
  #  dplyr::summarise(RSL_offset = unique(mean(RSL)))
  # data_TG <- merge(data_TG, Age_epoch_ref, by = "id", all = TRUE)
  # Cases where no data between 2000-2018 set the offset to 7000
  # data_TG$RSL_offset[is.na(data_TG$RSL_offset)] <- 7000
  
  # Updating the RSL to the shifted RSL value using PSMSL instructions
  data_TG$RSL <- data_TG$RSL - 7000 # data_TG$RSL_offset #
  
  #--Joining SL data with location names--
  annual_SL_tide_df <- base::merge(data_TG, file_list, by = "id", all = TRUE)
  #-- Removing sites which have a station flag raised as they are poor sites---
  annual_SL_tide_df <- annual_SL_tide_df %>%
    dplyr::filter(!stationflag == "Y") %>%
    tidyr::drop_na()
  
  # Remove the temporary file and directory
  unlink(temp_file)
  unlink(temp_dir, recursive = TRUE)
  
  # Annual Tidal Gauge data----
  annual_tidal_gauge_data_df <- annual_SL_tide_df %>%
    dplyr::select(
      Age, RSL, Latitude, name,
      # RSL_offset, Age_epoch_id,
      Longitude
    ) %>%
    dplyr::rename(SiteName = name) %>%
    # from mm --> m
    dplyr::mutate(RSL = RSL / 1000) %>%
    # Reordering by group
    dplyr::group_by(SiteName) %>%
    dplyr::arrange(SiteName, .by_group = TRUE) %>%
    dplyr::arrange(Age) %>%
    dplyr::mutate(data_type_id = "TideGaugeData")
  
  
  # Set the window size for the moving average (in this case, 10 years)
  # Rate of sedimentation for proxies when using continuous cores
  window_size <- sediment_average_TG
  
  # # Version A: Create a new column for the decade based on the midpoint of the rolling window
  # annual_tidal_gauge_data_df$rolling_avg <- zoo::rollapply(annual_tidal_gauge_data_df$RSL,
  #   width = window_size,
  #   FUN = mean,
  #   align = "right",
  #   fill = NA
  # )
  # # Calculate the decadal averages based on the rolling average
  # decadal_averages_TG <- annual_tidal_gauge_data_df %>% tidyr::drop_na()
  
  # Version B: Decadal Averages using simple method------
  decadal_averages_TG <-
    annual_tidal_gauge_data_df %>%
    dplyr::mutate(decade = (Age - 1) %/% window_size) %>%
    dplyr::group_by(decade, SiteName) %>%
    dplyr::summarise(
      # decade_meanRSL = mean(RSL)#,
      rolling_avg = mean(RSL),
      Age = max(Age) # ,
      # rows_site = dplyr::n()
    ) # Age=min(Age)
  
  # Using standard deviation of RSL over the decade as uncertainty----
  decadal_averages_TG <- decadal_averages_TG %>%
    dplyr::group_by(SiteName) %>%
    # dplyr::mutate(sd_TG = sd(decade_meanRSL))
    dplyr::mutate(sd_TG = sd(rolling_avg))
  
  #----- New df with decadal averages for tide gauges-----
  tidal_gauge_average_10_df <- base::merge(
    decadal_averages_TG,
    annual_tidal_gauge_data_df
  )
  
  #---Rsl & Age error for tidal gauge data----
  tidal_gauge_full_df <- tidal_gauge_average_10_df %>%
    dplyr::mutate(
      Age_er_min = 2.5, # years --> half a year/half a decade, 2.5 years either side.
      Age_er_max = 2.5
    ) %>%
    # dplyr::mutate(sd_TG = ifelse(is.na(sd_TG), 0.001, sd_TG)) %>%
    dplyr::group_by(SiteName) %>%
    dplyr::mutate(
      RSL_er_min =  0.001,
      RSL_er_max =  0.001 #sd_TG
    )
  
  tidal_gauge_full_df <- tidal_gauge_full_df %>%
    dplyr::mutate(Age = Age / 1000) %>%
    mutate(Age_er_max = Age_er_max/1000) %>% 
    mutate(Age_er_min = Age_er_min/1000) %>% 
    mutate(RSL_er_average = (RSL_er_min + RSL_er_max)/2) %>% 
    dplyr::mutate(RSL_annual = RSL) %>%
    # Change for Version A
    dplyr::mutate(RSL = rolling_avg) %>%
    # Change for Version B
    # dplyr::mutate(RSL = decade_meanRSL) %>%
    dplyr::select(!c(
      decade, # Age_epoch_id,RSL_offset,
      rolling_avg,
      RSL_annual
    ))
  
  # No user option here -> this is a must: Removing sites with only 2 points (20 years of data)-----
  decadal_TG_df <-
    tidal_gauge_full_df %>%
    dplyr::group_by(SiteName) %>%
    dplyr::filter(dplyr::n() > 2)
  
  
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
  
  SL_tide_site_df<- decadal_TG_df %>% 
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
  mat.distance <-  as.matrix(mat.distance)
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
    dplyr::select(!c(#RSL_annual,Age_epoch_id,RSL_offset,rows_site,decade_meanRSL,
                     sd_TG,
                     Indicator,Basin,Age_2_er_max,Age_2_er_min))%>% 
    mutate(SiteName = as.factor(SiteName))
  write_csv(SL_df,save_loc)
  cat("Tidal Gauge Data added to overall dataset")
  return(SL_df)
}