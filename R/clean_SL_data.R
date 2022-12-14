clean_SL_data <- function(SL_df,SL_data,save_loc){
  
  # Correcting some spelling mistakes in the dataset
  SL_df<-SL_data %>% 
    mutate(Region=factor(x=SL_data$Region,levels = unique(SL_data$Region))) %>% 
    mutate(Site = replace(Site, Site == "Sand Point Russian", "Sand Point")) %>% 
    mutate(Site = replace(Site, Site == "Sand Point VC", "Sand Point")) %>%    
  #--- Creating Cedar Island site as they are close by--
  # For small scale local factors could separate them again.
    mutate(Site = replace(Site, Site == "Tump Point", "Cedar Island")) %>% 
    mutate(Site = replace(Site, Site == "Sand Hill Point", "Cedar Island")) %>% 
  #--- Combining Les Sillons & Bassin: Andy recommendation--
    mutate(Site = replace(Site, Site == "Bassin", "Les Sillons"))
  
  #----Remove Nydick data as its a different kind of data-----
  nydick <- SL_df %>% filter(Reference == "Nydick et al. (1995)")
  SL_df <- anti_join(SL_df,nydick)
  #---Joining the Region & site----
  SL_df <- SL_df %>% 
    mutate(SiteName = as.factor(paste0(Site,",","\n"," ", Region)))
  
  # Filtering out Brazilian data as they are similar to index points
  Brazilian_sites <- SL_df %>% filter(Region %in% 
                                        c("Santa Catarina Brazil","Rio de Janeiro Brazil"))
  
  # Updated SL dataframe
  SL_df <- anti_join(SL_df, Brazilian_sites)
  
  # Filtering out Hammock as they are similar to index points
  hammock_site <- SL_df %>% filter(Site %in% 
                                        c("Hammock River"))
  
  # Updated SL dataframe
  SL_df <- anti_join(SL_df, hammock_site)
  
  # Changing long & lat values at each "long" record site to be the same
  SL_df <- SL_df %>%
    group_by(SiteName) %>% 
    mutate(Longitude = first(Longitude),
           Latitude = first(Latitude))
  
  # Some sites missing uncertainties
  SL_df <- SL_df %>%
    group_by(SiteName) %>%
    #--Correcting when Age uncertainties = 0 with 1 years----
    mutate(Age_2_er_max = ifelse(Age_2_er_max == 0,1,Age_2_er_max)) %>% 
    mutate(Age_2_er_min = ifelse(Age_2_er_min == 0,1,Age_2_er_min)) %>% 
    #Age uncertainty is 2 standard deviations so divide by 2
    mutate(Age_er_max = Age_2_er_max/2) %>% 
    mutate(Age_er_min = Age_2_er_min/2) %>% 
    mutate(Age_er_average = (Age_er_max+Age_er_min)/2) %>% # Age error not same on both sides
    mutate(RSL_er_average = (RSL_er_min + RSL_er_max)/2) #RSL error not same on both sides
  
  #--Removing sites with only index points---
  #SL_df <- SL_df %>%  group_by(SiteName) %>% filter(n() > 5) 
  
  #----Setting up the uncertainities & divide Age by 1000----
  SL_df <- SL_df %>%
    mutate(Age = Age/1000) %>%
    mutate(Age_er_max = Age_er_max/1000) %>% 
    mutate(Age_er_min = Age_er_min/1000) %>% 
    mutate(Age_er_average = Age_er_average/1000)
  
  # Just looking at North Atlantic, specifically North America (22 sites)
  #SL_df <- SL_df %>% filter(Basin == "North Atlantic")#Not all sites have Basin filled in
  SL_df <- SL_df %>% filter(Region %in% c("Connecticut","Magdelen Islands",
                                          "Newfoundland","New Jersey","New York",
                                          "North Carolina","Nova Scotia",
                                          "Rhode Island","Massachusetts","Quebec",
                                          "Maine","Florida"))
  
  #--Truncating the time---- (22 Sites)
  SL_df <- SL_df %>%  filter(Age > -1)
  
  write_csv(SL_df,file = save_loc)
  cat("SL dataset for North America is cleaned")
  return(SL_df)
  
}

