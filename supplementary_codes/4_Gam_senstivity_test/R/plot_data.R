plot_data<-function(SL_df,save_name,n_sites)
{
  SL_df_proxy <- SL_df %>% filter(data_type_id == "ProxyData")
  #------Ordering the sites-------
  all_data_sites<-factor(paste0(SL_df_proxy$Longitude, SL_df_proxy$Longitude), labels = 1:20)#n_sites)#For ordering sites
  SL_df_proxy <- cbind(SL_df_proxy,all_data_sites=all_data_sites)
  order_sites <- SL_df_proxy %>%  group_by(SiteName,all_data_sites) %>%
    dplyr::summarise(n = n()) %>%
    arrange(all_data_sites)
  SL_df_proxy$SiteName <- factor(SL_df_proxy$SiteName,levels = unique(order_sites$SiteName))
  
  plot_data <- ggplot() +
    geom_rect(data = SL_df_proxy, aes(
      xmin = Age*1000 - Age_er_min*1000, xmax = Age*1000 + Age_er_max*1000,
      ymin = RSL - RSL_er_min, ymax = RSL + RSL_er_max,
      fill = "grey",
    ), alpha = 2) +
    
    geom_point(data = SL_df_proxy, 
               aes(y = RSL, x = Age*1000,colour="black"), 
               size = 0.4, stroke = 0)+
               #size = 0.001) +
    xlab("Year (CE)") +
    ylab("Relative Sea Level (m)") +
    theme_bw()+
    labs(colour = "")+
    theme(strip.text.x = element_text(size = 5),
          strip.background =element_rect(fill=c("white")))+
    scale_fill_manual('',
                      values = 'grey', 
                      labels = expression(paste("1",sigma," " ,"RSL and age uncertainty")),
                      guide = guide_legend(override.aes = list(alpha = 0.7))) +
    scale_colour_manual(values = c("black"),
                        labels=c("Proxy Records"))+
    facet_wrap(~SiteName)+
    #ggtitle("Tide Gauge and Proxy Data for East Coast of North America")+
    theme(legend.position=c(1,0),legend.justification = c(1,0), legend.box = "horizontal")+
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(plot.title = element_text(size=18,face="bold"),
          axis.title=element_text(size=12,face="bold"),
          legend.text=element_text(size=10))
  plot_data
  ggsave(plot_data, file = save_name, width = 10, height = 6)
  
  #----Plotting the TG------
  SL_df_TG <- SL_df %>% filter(data_type_id == "TideGaugeData")
  #------Ordering the TG sites-------
  all_data_sites<-factor(paste0(SL_df_TG$Longitude, SL_df_TG$Longitude), labels = 1:52)
  SL_df_TG <- cbind(SL_df_TG,all_data_sites=all_data_sites)
  order_sites <- SL_df_TG %>%  group_by(SiteName,all_data_sites) %>%
    dplyr::summarise(n = n()) %>%
    arrange(all_data_sites)
  SL_df_TG$SiteName <- factor(SL_df_TG$SiteName,levels = unique(order_sites$SiteName))

  
  # Filtering for TG sites with most amount of data----
  SL_df_TG <- SL_df_TG %>% 
    group_by(SiteName) %>% 
    mutate(Age_range = max(Age*1000) - min(Age*1000))
  # Only plotting sites with more than 100 years data----
  plot_TG_data_100 <- ggplot(subset(SL_df_TG,Age_range >= 100)) +
    geom_rect(aes(
      xmin = Age*1000 - Age_er_min*1000, xmax = Age*1000 + Age_er_max*1000,
      ymin = RSL - RSL_er_min, ymax = RSL + RSL_er_max,
      fill = "grey",
    ), alpha = 3) +
    geom_point(aes(y = RSL, x = Age*1000,
                   colour="black"), size = 0.5) +
    xlab("Year (CE)") +
    ylab("Relative Sea Level (m)") +
    theme_bw()+
    labs(colour = "")+
    theme(strip.text.x = element_text(size = 7),
          strip.background =element_rect(fill=c("white")))+
    scale_fill_manual('',
                      values = 'grey', 
                      labels = expression(paste("1",sigma," " ,"RSL and age uncertainty")),
                      guide = guide_legend(override.aes = list(alpha = 0.7))) +
    scale_colour_manual(values = c("black"),
                        labels=c("Tide Gauge Data"))+
    theme(legend.position=c(1,0),legend.justification = c(1,0), legend.box = "horizontal")+
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(plot.title = element_text(size=18,face="bold"),
          axis.title=element_text(size=12,face="bold"),
          legend.text=element_text(size=8),
          legend.title=element_blank(),
          legend.margin=margin(c(1,1,1,1)))+
    facet_wrap(~SiteName)
  plot_TG_data_100
  ggsave(plot_TG_data_100, file = "fig/TG_plot_1_over_100yrs.pdf", width = 10, height = 6)
  
  
  cat("Plot of raw data saved to figure folder \n")
}

plot_map<-function(SL_df,save_name)
{
  #---Map sites----
  cand_usa_states <- rnaturalearth::ne_states(country = c('United States of America',"Canada"), returnclass = 'sf')
  sites <- SL_df %>%  mutate(SiteName = as.character(SiteName))%>%
                               dplyr::select(Longitude,Latitude,SiteName) %>% unique()
  sites_4 <- SL_df %>% filter(SiteName %in% c("East River Marsh,\n Connecticut",
                                              "Swan Key,\n Florida",
                                              "Placentia,\n Newfoundland",
                                              "Cedar Island,\n North Carolina")) %>%
    dplyr::select(Longitude,Latitude,SiteName) %>% unique()

  
  NAmap_new <- 
    ggplot(data=cand_usa_states,size = 0.5, alpha = 0) +
    geom_sf(color="darkgrey",fill="darkseagreen2") +
    geom_label_repel(data=sites,aes(x=Longitude,y=Latitude,
                                    label = ifelse(SiteName %in% sites_4$SiteName, SiteName, '')),
                     min.segment.length = 0,
                     max.overlaps = 50,fontface = "bold", size = 3)+
    theme_minimal()+
    coord_sf(xlim = c(-85, -50), ylim = c(24, 52), expand = FALSE) +
        geom_point(data=sites_4,aes(x=Longitude,y=Latitude),
               size=5,colour='red',stroke = 1,shape = 1)+
    geom_point(data=SL_df,aes(x=Longitude,y=Latitude,
                              shape = data_type_id,colour = data_type_id),size = 2)+#colour="black",

    labs(x="Longitude",y="Latitude")+
    annotation_scale(location = "br", width_hint = 0.5) +
    geom_text(aes(x = -60,y=35),label = "Atlantic Ocean",
              fontface="italic",size = 6)+
    annotation_north_arrow(location = "br", which_north = "true",
                           pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                           style = north_arrow_fancy_orienteering) +
    theme(plot.title = element_text(size=18,face="bold"),
          panel.grid.major = element_line(color = "lightblue1", linetype = "dashed", size = 0.05),
          panel.background = element_rect(fill = "aliceblue"),
          axis.title=element_text(size=10,face="bold"),
          legend.position = "none",
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    scale_colour_manual(values = c("black","saddlebrown"),
                        labels=c("Proxy Record Sites","Tide-Gauge Sites"),
                        guide = TRUE)+
    scale_shape_manual(values = c(16,17),
                       labels=c("Proxy Record Sites","Tide-Gauge Sites"),
                       guide = TRUE)+
    theme(legend.position=c(0,1),legend.justification = c(0,1), 
          legend.box.background = element_rect(color="black",fill="white", size=0.5),
          legend.box = "horizontal",legend.margin = margin(6, 6, 6, 6),
          legend.title=element_blank(),
          legend.text = element_text(size = 10,face="bold"))+
    guides(color = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3)))
  NAmap_new
  ggsave(NAmap_new, file = save_name, width = 10, height = 6)
  
  cat("Plot of Map of sites saved to figure folder \n")
}

plot_4_data<-function(SL_df,save_name,n_sites)
{
  SL_df_proxy <- SL_df %>% filter(data_type_id == "ProxyData")
  #--- Filter for plotting 4 sites----
  SL_df_proxy <- SL_df_proxy %>% filter(SiteName %in% c("East River Marsh,\n Connecticut",
                                                        "Swan Key,\n Florida",
                                                        "Placentia,\n Newfoundland",
                                                        "Cedar Island,\n North Carolina")) %>%
    mutate(SiteName = as.factor(SiteName))
  
  #------Ordering the sites-------
  all_data_sites<-factor(paste0(SL_df_proxy$Longitude, SL_df_proxy$Longitude), labels = 1:4)#20)#n_sites)#For ordering sites
  SL_df_proxy <- cbind(SL_df_proxy,all_data_sites=all_data_sites)
  order_sites <- SL_df_proxy %>%  group_by(SiteName,all_data_sites) %>%
    dplyr::summarise(n = n()) %>%
    arrange(all_data_sites)
  SL_df_proxy$SiteName <- factor(SL_df_proxy$SiteName,levels = unique(order_sites$SiteName))
  
  plot_data <- ggplot() +
    geom_rect(data = SL_df_proxy, aes(
      xmin = Age*1000 - Age_er_min*1000, xmax = Age*1000 + Age_er_max*1000,
      #xmin = Age*1000 - Age_er_average*1000, xmax = Age*1000 + Age_er_average*1000,
      ymin = RSL - RSL_er_min, ymax = RSL + RSL_er_max,
      fill = "grey",
    ), alpha = 0.7) +
    
    geom_point(data = SL_df_proxy, aes(y = RSL, x = Age*1000,colour="black"), size = 0.3) +
    xlab("Year (CE)") +
    ylab("RSL (m)") +
    theme_bw()+
    labs(colour = "")+
    theme(strip.text.x = element_text(size = 5),
          strip.background =element_rect(fill=c("white")))+
    scale_fill_manual('',
                      values = 'grey', 
                      labels = expression(paste("1",sigma," " ,"RSL and age uncertainty")),
                      guide = guide_legend(override.aes = list(alpha = 0.7))) +
    scale_colour_manual(values = c("black"),
                        labels=c("Proxy Records"))+
    facet_wrap(~SiteName)+
    #ggtitle("Tide Gauge and Proxy Data for East Coast of North America")+
    theme(legend.position=c(0.95,0.01),
          legend.justification = c(1,0),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title=element_blank(),
          legend.margin=margin(c(1,1,1,1)))+
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(plot.title = element_text(size=18,face="bold"),
          axis.title=element_text(size=12,face="bold"),
          legend.text=element_text(size=8))
  plot_data
  ggsave(plot_data, file = save_name, width = 10, height = 6)
  
  cat("Plot of raw data saved to figure folder \n")
}


