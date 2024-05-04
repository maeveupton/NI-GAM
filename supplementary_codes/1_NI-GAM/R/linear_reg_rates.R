linear_reg_rates<- function(SL_df){#, save_loc){
  SL_df_filter <- SL_df %>%  
    #filter(!Age > 1.850) # Ignoring recent human influences to SL rise
    filter(!Age > 1.800) # Ignoring recent human influences to SL rise
  
  # Doing linear regression on rest of data
  SL_df_lm <- SL_df_filter %>% 
    group_by(SiteName) %>% 
    mutate(
      data_lm_slope = round(lm(RSL ~ Age)$coefficients[[2]], 2),
      data_lm_slope_err = summary(lm(RSL ~ Age))$coefficients[2, 4]
    )
  
  plot_lm <- ggplot() +
    facet_wrap(~SiteName)+
    geom_point(data = SL_df_filter, aes(y = RSL, x = Age*1000,colour = data_type_id), size = 0.5) +
    geom_rect(data = SL_df_filter, aes(
      xmin = Age*1000 - Age_er_min*1000, xmax = Age*1000 + Age_er_max*1000,
      #xmin = Age - Age_2_er_min, xmax = Age + Age_2_er_max,
      ymin = RSL - RSL_er_min, ymax = RSL + RSL_er_max
    ), alpha = 0.2) +
    xlab("Age (CE)") +
    geom_smooth(data = SL_df_filter, aes(y = RSL, x = Age*1000),colour = 'black',method = "lm",se = TRUE)+
    ylab("RSL (m MTL)") +
    theme_bw()+
    labs(colour = "")+
    theme(strip.background =element_rect(fill="aquamarine3"))+
    ggtitle("Linear Regression for Proxy Data for east coast of North America")
  plot_lm
  ggsave(plot_lm, file = "fig/lm_SL.pdf", width = 10, height = 6)

  
  # Table of GIA rate vs lm rate from proxy data
 lm_slopes <- SL_df_lm %>%
   dplyr::select(SiteName, data_lm_slope,data_lm_slope_err) %>%
    #dplyr::select(SiteName, Longitude,Latitude, data_lm_slope,data_lm_slope_err) %>%
    unique()
 
   #write_csv(GIA_vs_lm_slopes,sav_loc)
  return(lm_slopes)
  cat("Linear Regression on the data excluding tide gauges \n")
  
}
# # Matching GIA rate from proxy to TG data
# SL_df_TG_unique <- SL_df %>%  filter(data_type_id == "TideGaugeData") %>% 
#   dplyr::select(SiteName,Longitude,Latitude) %>% 
#   unique()
# mat.distance<- distm(lm_slopes[,2:3],SL_df_TG_unique[,2:3],
#                      fun = distGeo)
# mat.distance <-  as.data.frame(mat.distance)
# #--finding row mins & corresponding tidal gauge--
# rownames(mat.distance) = lm_slopes$SiteName#as.character(SL_df_lon_lat$SiteName %>% unique)
# colnames(mat.distance) = SL_df_TG_unique$SiteName
# #--finding row mins --
# dist_data_grid <- t(sapply(seq(nrow(mat.distance)), function(z) {
#   js <- order(mat.distance[z,])[1:2]
#   c(rownames(mat.distance)[z], colnames(mat.distance)[js[1]], mat.distance[z,js[1]])
# }))
# 
# dist_data_grid_df <- as.data.frame(dist_data_grid)
# colnames(dist_data_grid_df) <- c("proxy_site", "TG_site", "min_dist1")


