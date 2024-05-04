# Comparing different number of knots for regional component:
dataframe_mod_outputs <- function(model_run,name_df){
#-----------Plotting Overall Output----------
mu_post <- model_run$BUGSoutput$sims.list$mu
#---Get estimates and uncertainty bounds--
RSL_mod<-apply(mu_post,2,mean)
RSL_mod_upr<-apply(mu_post,2,quantile,probs=0.025)
RSL_mod_lwr<-apply(mu_post,2,quantile,probs=0.975)
lwr_50<- apply(mu_post,2,quantile,probs=0.25)
upr_50<- apply(mu_post,2,quantile,probs=0.75)

#-----Create data frame for plotting for total-----
mod_output_df<-data.frame(RSL_mod,RSL_mod_upr,
                          RSL_mod_lwr,
                          lwr_50,upr_50,
                          SL_df$Age,
                          SL_df$SiteName,SL_df$Longitude,SL_df$Latitude,
                          ID = "Total Posterior Model",
                          knot_setting = name_df)
names(mod_output_df)<-c("RSL","upr","lwr","upr_50","lwr_50", "Age",
                        "SiteName","Longitude","Latitude",
                        "ID","knot_setting")
# Filter for 4 sites 
mod_output_df <- mod_output_df %>% 
  filter(SiteName %in% c(
    "East River Marsh,\n Connecticut",
    "Swan Key,\n Florida",
    "Placentia,\n Newfoundland",
    "Cedar Island,\n North Carolina")) %>%
  mutate(SiteName = as.factor(SiteName))

#------Ordering the sites total component-------
all_data_sites<-factor(paste0(mod_output_df$Longitude,
                              mod_output_df$Longitude), labels = 1:4)
mod_output_df <- cbind(mod_output_df,all_data_sites=all_data_sites)
order_sites <- mod_output_df %>%  group_by(SiteName,all_data_sites) %>%
  dplyr::summarise(n = n()) %>%
  arrange(all_data_sites)
mod_output_df$SiteName <- factor(mod_output_df$SiteName,
                                 levels = unique(order_sites$SiteName))


#------------Regional Component: Spline in Time--------------
time_component <- model_run$BUGSoutput$sims.list$r
RSL_mod_NI<-apply(time_component,2,mean)
RSL_mod_lwr_NI<-apply(time_component,2,quantile,probs=0.025)
RSL_mod_upr_NI<-apply(time_component,2,quantile,probs=0.975)
lwr_50<-apply(time_component,2,quantile,probs=0.25)
upr_50<-apply(time_component,2,quantile,probs=0.75)

time_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,RSL_mod_lwr_NI,
                              lwr_50,upr_50,
                              Age = SL_df$Age,
                              SiteName = SL_df$SiteName,
                              ID = "Regional Component",
                              knot_setting = name_df)
names(time_component_df)<-c("RSL","upr","lwr","lwr_50",
                            "upr_50", "Age","SiteName","ID","knot_setting")

time_component_df <- time_component_df %>%
  filter(SiteName %in% c(
    "East River Marsh,\n Connecticut",
    "Swan Key,\n Florida",
    "Placentia,\n Newfoundland",
    "Cedar Island,\n North Carolina")) %>%
  mutate(SiteName = as.factor(SiteName))
#------Non-Linear Local Component: Spline in Space Time------
space_time_component <- model_run$BUGSoutput$sims.list$l
RSL_mod_NI<-apply(space_time_component,2,mean)
RSL_mod_lwr_NI<-apply(space_time_component,2,quantile,probs=0.025)
RSL_mod_upr_NI<-apply(space_time_component,2,quantile,probs=0.975)
lwr_50<-apply(space_time_component,2,quantile,probs=0.25)
upr_50<-apply(space_time_component,2,quantile,probs=0.75)
space_time_component_df<-data.frame(RSL_mod_NI,RSL_mod_upr_NI,
                                    RSL_mod_lwr_NI,
                                    lwr_50,upr_50,
                                    Age = SL_df$Age,
                                    Longitude = SL_df$Longitude,
                                    Latitude = SL_df$Latitude,
                                    SiteName = SL_df$SiteName,
                                    ID = "Non-Linear Local Component",
                                    knot_setting = name_df)
names(space_time_component_df)<-c("RSL","upr", "lwr","lwr_50","upr_50", "Age",
                                  "Longitude","Latitude","SiteName","ID","knot_setting")

space_time_component_df <- space_time_component_df %>%
  filter(SiteName %in% c(
    "East River Marsh,\n Connecticut",
    "Swan Key,\n Florida",
    "Placentia,\n Newfoundland",
    "Cedar Island,\n North Carolina")) %>%
  mutate(SiteName = as.factor(SiteName))
# --Filtering for 4 sites---
SL_df <-SL_df %>%   filter(SiteName %in% c("East River Marsh,\n Connecticut",
                                           "Swan Key,\n Florida",
                                           "Placentia,\n Newfoundland",
                                           "Cedar Island,\n North Carolina")) %>%
  mutate(SiteName = as.factor(SiteName))

#------Ordering the sites original Data set-------
all_data_sites<-factor(paste0(SL_df$Longitude, SL_df$Longitude), labels = 1:4)
SL_df <- cbind(SL_df,all_data_sites=all_data_sites)
order_sites <- SL_df %>%  group_by(SiteName,all_data_sites) %>%
  dplyr::summarise(n = n()) %>%
  arrange(all_data_sites)
SL_df$SiteName <- factor(SL_df$SiteName,levels = unique(order_sites$SiteName))
space_time_component_df$SiteName <- factor(SL_df$SiteName,levels = unique(order_sites$SiteName))


  df_list <- list(time_component_df,mod_output_df,space_time_component_df,SL_df)
return(df_list)
}

