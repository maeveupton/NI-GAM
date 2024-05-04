# # Comparing knot settings------
# Non-linear local knots
model_run6 <- readRDS("output/model_run_regional_noise.rds")
n_seg6 <- dataframe_mod_outputs(model_run = model_run6,name_df = "262 knots")#42mins
model_run3 <- readRDS("output/model_run_regional_noisenseg_st = 3.rds")
n_seg3 <- dataframe_mod_outputs(model_run = model_run3,name_df = "95 knots")# 40 mins
model_run9 <- readRDS("output/model_run_regional_noisenseg_st = 9.rds")
n_seg9 <- dataframe_mod_outputs(model_run = model_run9,name_df = "445 knots")#44
#
# Plotting tot mod fit for 4 sites with different settings:
# --Filtering for 4 sites---
# SL_df_4 <-SL_df %>%   filter(SiteName %in% c("East River Marsh,\n Connecticut",
#                                            "Swan Key,\n Florida",
#                                            "Placentia,\n Newfoundland",
#                                            "Cedar Island,\n North Carolina")) %>%
#   mutate(SiteName = as.factor(SiteName))
SL_df_4<-  n_seg6[[4]]
 
knot_plot_tot_nonlin <- ggplot()+
  geom_point(data = SL_df_4, aes(y = RSL, x =  Age*1000), size = 0.5) +
  geom_rect(data = SL_df_4, aes(
    xmin = Age*1000 - Age_er_average*1000, xmax = Age*1000 + Age_er_average*1000,
    ymin = RSL - RSL_er_average, ymax = RSL + RSL_er_average
  ), alpha = 0.2) +
  geom_line(data = n_seg6[[2]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg6[[2]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  geom_line(data = n_seg3[[2]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg3[[2]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  geom_line(data = n_seg9[[2]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg9[[2]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  facet_wrap(~SiteName)+
  labs(color='Non-Linear Local Knot Setting',fill = 'Non-Linear Local Knot Setting') +
  xlab("Year (CE)") +
  ylab("RSL (m)") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 14),
        strip.background =element_rect(fill=c("white")))

knot_plot_tot_nonlin
ggsave("fig/knot_plot_tot_nonlinloc_nseg.pdf",
       knot_plot_tot_nonlin, width = 10, height = 6)

knot_plot_nonlin <- ggplot()+
  geom_line(data = n_seg6[[3]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg6[[3]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  geom_line(data = n_seg3[[3]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg3[[3]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  geom_line(data = n_seg9[[3]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg9[[3]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  xlab("Year (CE)") +
  ylab("Sea Level (m)") +
  labs(color='Non-Linear Local Knot Setting',fill = 'Non-Linear Local Knot Setting') +
  theme_bw()+
  theme(strip.text.x = element_text(size = 14),
        strip.background =element_rect(fill=c("white")))+
  facet_wrap(~SiteName)
knot_plot_nonlin
ggsave("fig/knot_plot_nonlin_nseg.pdf",
       knot_plot_nonlin, width = 10, height = 6)



# Regional knots
model_run20 <- readRDS("output/model_run_regional_noise.rds")
n_seg20 <- dataframe_mod_outputs(model_run = model_run20,name_df = "23 knots")#42mins
model_run18 <- readRDS("output/model_run_regional_noisenseg_t = 18.rds")
n_seg18 <- dataframe_mod_outputs(model_run = model_run18,name_df = "21 knots")# 40 mins
#model_run10 <- readRDS("output/model_run_regional_noisenseg_t = 10.rds")
#n_seg10 <- dataframe_mod_outputs(model_run = model_run10,name_df = "13 knots")# 40 mins
model_run50 <- readRDS("output/model_run_regional_noisenseg_t = 50.rds")
n_seg50 <- dataframe_mod_outputs(model_run = model_run50,name_df = "53 knots")#44
#model_run30 <- readRDS("output/model_run_regional_noisenseg_t = 30.rds")
#n_seg30 <- dataframe_mod_outputs(model_run = model_run30,name_df = "33 knots")#44

# Plotting tot mod fit for 4 sites with different settings:
knot_plot_tot <- ggplot()+
  geom_point(data = SL_df_4, aes(y = RSL, x =  Age*1000), size = 0.5) +
  geom_rect(data = SL_df_4, aes(
    xmin = Age*1000 - Age_er_average*1000, xmax = Age*1000 + Age_er_average*1000,
    ymin = RSL - RSL_er_average, ymax = RSL + RSL_er_average
  ), alpha = 0.2) +
  geom_line(data = n_seg20[[2]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg20[[2]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  geom_line(data = n_seg18[[2]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg18[[2]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  
  #geom_line(data = n_seg10[[2]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  #geom_ribbon(data=n_seg10[[2]],
  #            aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+

  #geom_line(data = n_seg15[[2]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  #geom_ribbon(data=n_seg15[[2]],
   #           aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  geom_line(data = n_seg50[[2]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg50[[2]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  facet_wrap(~SiteName)+
  labs(color='Regional Knot Setting',fill = 'Regional Knot Setting') +
  xlab("Year (CE)") +
  ylab("RSL (m)") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 14),
        strip.background =element_rect(fill=c("white")))

knot_plot_tot
ggsave("fig/knot_plot_tot_reg_nseg.pdf",
       knot_plot_tot, width = 10, height = 6)

knot_plot_reg <- ggplot()+
  geom_line(data = n_seg20[[1]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg20[[1]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  #geom_line(data = n_seg10[[1]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  #geom_ribbon(data=n_seg10[[1]],
  #            aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+

  geom_line(data = n_seg18[[1]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg18[[1]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  
  #geom_line(data = n_seg15[[1]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  #geom_ribbon(data=n_seg15[[1]],
  #            aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  geom_line(data = n_seg50[[1]],aes(x = Age*1000, y = RSL,colour = knot_setting))+
  geom_ribbon(data=n_seg50[[1]],
              aes(y = RSL,ymin=lwr,ymax=upr,x=Age*1000,fill = knot_setting),alpha=0.2)+
  xlab("Year (CE)") +
  ylab("Sea Level (m)") +
  labs(color='Regional Knot Setting',fill = 'Regional Knot Setting') +
  theme(strip.text.x = element_text(size = 14),
        strip.background =element_rect(fill=c("white")))+
  theme_bw()
knot_plot_reg
ggsave("fig/knot_plot_reg_nseg.pdf",
       knot_plot_reg, width = 10, height = 6)
