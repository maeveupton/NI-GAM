####----Noisy Input GAM for RSL along East Coast of North America---------
# 14th November 2022
# Data: Proxies & separate tide gauges all tg along NA coast which are 1deg from proxy & have greater than 20 years of observations
# JAGS Set Up:SP Time(Regional) + SP SpaceTime(Local) + [RandomEffect Site & Age(GIA) + RandomEffect Site(y offset)]
# Priors: y offset ,GIA centred on rates from finding the linear regression through data pre industrial revolution
# GIA rates for TG from ICE5G
# Model Fitting: Running model first without local and then obtaining value for b_r & b_h and rerun
#b_r is the coefficients of the regional basis functions & b_h is the site specific offset(varying intercept)

# Clear workspace
rm(list = ls())

#---------Set working directory--------------
setwd('/users/research/mupton/1. NIGAM_paper_2022/1. NIGAM_NA_TG_separate_more_sites_updateGIA')
#----------Load packages---------------------
library(tidyverse)
library(geosphere) #distm
library(R2jags)
library(splines)
library(cowplot) # for griding plots at end
library(dbarts)# dummy matrix creation
library(ncdf4) # package for netcdf manipulation
library("rnaturalearth")#n_states for map
library("ggrepel") #label_repel on map
library("ggspatial")#annotation_scale in map
library("xtable")#Latex tables


#-----Functions-----
source("R/clean_SL_data.R")
source("R/clean_tidal_gauge_data.R")
source("R/linear_reg_rates.R")
source("R/match.closest.R")
source("R/add_GIA_rate.R")
source("R/plot_data.R")
source("R/bs_bbase.R")
source("R/spline_basis_function.R")
source("R/run_spline.R")
source("R/run_noise_spline.R")
source("R/plot_results.R")
source("R/add_noise.R")
source("R/plot_4_sites_results.R")
source("R/acceleration_fun.R")
source("R/rate_of_change_fun.R")


#---Read in SL data from Andy's data----
SL_data <- read_csv("https://www.dropbox.com/s/43kbaoxqmnrprmv/Common%20Era%20Database%202022_a.csv?dl=1")
new_col_names <- c('Basin', 'Region','Site','Reference','Indicator','Latitude','Longitude',
                   'RSL','RSL_er_max','RSL_er_min','Age','Age_2_er_max','Age_2_er_min')
names(SL_data) <- new_col_names

#---Cleaning SL data----
SL_df <- clean_SL_data(SL_df=SL_data,SL_data=SL_data,
              save_loc="data/SL_df_clean.csv")

#--Clean tidal gauge data---
SL_df <- clean_tidal_gauge_data(SL_df=SL_df,
                       #--Data from my github account--
                       path_to_data = "https://raw.githubusercontent.com/maeveupton/tide_gauge_data_PSMSL/main/annual_SL_tide_df.csv",
                       save_loc="data/SL_tidal_gauge_df_new.csv")

#---Linear regression to find rates of data pre 1800-----
lm_data_rates <- linear_reg_rates(SL_df=SL_df)
SL_df <- left_join(SL_df,lm_data_rates, by = "SiteName")

#---Adding GIA rates-----
SL_df <- add_GIA_rate(SL_df=SL_df, save_loc = "data/SL_df_GIA.csv")

SL_df <- SL_df %>%
  mutate(
    data_lm_slope = ifelse(data_type_id=="TideGaugeData",ICE5_GIA_slope,data_lm_slope),
    data_lm_slope_err = ifelse(data_type_id=="TideGaugeData",0.3,data_lm_slope_err))# Simon
    # Not right
    #data_lm_slope = ifelse(data_type_id=="TideGaugeData",GIA_caron,data_lm_slope),
    #data_lm_slope_err = ifelse(data_type_id=="TideGaugeData",GIA_caron_sd,data_lm_slope_err))
         
SL_df <- SL_df %>% mutate(SiteName = as.factor(SiteName))
# #--Plot Data---
# plot_4_data(SL_df = SL_df,
#             save_name = "fig/SL_data_plot_4_sites_1000.pdf",
#             n_sites= length(unique(SL_df$SiteName)))
# plot_data(SL_df = SL_df,
#           save_name = "fig/SL_data_plot_all_sites_1000.pdf",
#           n_sites= length(unique(SL_df$SiteName)))
# plot_map(SL_df = SL_df,
#          save_name = "fig/map_data_all_sites.pdf")

# #---Latex tables for TG data----
# SL_df_TG <- SL_df %>% filter(data_type_id == "TideGaugeData")
# ref_TG_table <- SL_df_TG %>%
#   dplyr::select(Longitude,Latitude,SiteName,data_lm_slope) %>% unique() %>% as_tibble()
# print(xtable(ref_TG_table,type = "latex"),file = "SL_TG_data_references.tex",include.rownames=FALSE)
# 
# #---Latex tables for proxy data----
# SL_df_proxy <- SL_df %>% filter(data_type_id == "ProxyData")
# ref_proxy_table <- SL_df_proxy %>%
#   dplyr::select(SiteName,data_lm_slope,ICE5_GIA_slope,ICE6_GIA_slope,GIA_caron) %>%
#   unique() %>% as_tibble()
# print(xtable(ref_proxy_table,type = "latex"),file = "SL_proxy_data_references.tex",include.rownames=FALSE)

#---Spline Basis Functions----
basis_fun_list <- spline_basis_function(SL_df = SL_df)

set_up_option <- " "
#save_option <- " "
#save_option <- "_test_caron_GIA"
save_option <- "_ICE5G_0.3"

#--- Stage 1: Run JAGS for model without noise & without Local(6min)-----
# run_spline(SL_df = SL_df,
#            save_option= save_option,
#            basis_fun_list)

#---Check convergence----
model_run_1 <- readRDS(paste0("output/model_run_no_local_no_noise", save_option, ".rds"))
#plot(model_run_1)

# Plot the regional term
# tibble(
#   site = SL_df$SiteName,
#   age = SL_df$Age,
#   RSL = SL_df$RSL,
#   regional_mean = model_run_1$BUGSoutput$mean$regional,
#   regional_low = apply(model_run_1$BUGSoutput$sims.list$regional,
#                        2, 'quantile', 0.05),
#   regional_high = apply(model_run_1$BUGSoutput$sims.list$regional,
#                         2, 'quantile', 0.95)
# ) %>% ggplot(aes(x = age, y = regional_mean)) +
#   geom_ribbon(aes(ymin=regional_low, ymax=regional_high),
#               fill = 'red', alpha = 0.5) +
#   geom_line() 

#--- Stage 2: Adding Age uncertainty-----
SL_df <- add_noise(jags.data = jags.data,
                   jags.file = jags.file,
                   model_run = model_run_1,
                   SL_df=SL_df,save_csv = "data/SL_df_noise.csv")

#--- Stage 3: Run JAGS for model with noise (20mins)-----
# run_noise_spline(SL_df = SL_df,
#            save_option= save_option,
#            basis_fun_list,
#            model_run_1 = model_run_1)

model_run_2 <- readRDS(paste0("output/model_run_regional_noise", save_option, ".rds"))

#-----------Plotting results----------
plot_results(SL_df = SL_df,model_run = model_run_2,
             save_option = save_option,basis_fun_list = basis_fun_list,
             set_up_option=set_up_option,
             n_sites = length(unique(SL_df$SiteName)))

#------Rate of Change-------
rate_of_change_fun(SL_df = SL_df,model_run= model_run_2,
                   save_option=save_option)

#-----Plotting just 4 sites------
set_up_option <- " "
#save_option <- "caron_GIA_4_sites"
save_option <- "ICE5G_0.3_GIA_4_sites"

plot_4_sites_results(SL_df = SL_df,
                     model_run = model_run_2,
                     save_option = save_option,
                     basis_fun_list=basis_fun_list,
                     set_up_option=set_up_option)
 