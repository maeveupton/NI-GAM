####----Noisy Input GAM for RSL along Atlantic coast of North America---------

# Data: NA Atlantic coast = Proxies & tide gauges(1 deg from proxy & greater than 20 observations)

# Process Level: f(x,t) = r(t) + g(z_x) + h(z_x) + l(x,t) + epilson
#             r(t) = r = Regional component(B-spline in time)
#             g(z_x) = g_z_x = Linear Local component(Random Effect)
#             h(z_x) = h_z_x = Site Specific vertical offset(Random Effect) 
#             l(x,t) = l = Non-Linear Local component(B-spline in space time)

# Priors: Linear local component informed by linear regression through data prior to 1800 & for TG rate from ICE5G with uncertainty from Engelhart 2011

# Model Fitting:1. Model run first time without Non-linear local component. Estimate for b_r & sigma_h is obtained. 
#               2. Noisy Input corrective variance term estimated.
#               3. Model re-run with non-linear local component, corrective variance term included and priors for Regional component & site specific vertical offset are informed

# Clear workspace
rm(list = ls())

#----------Load packages--------------------
library(tidyverse)
library(geosphere) #distm
library(R2jags)
library(splines)
library(dbarts)# dummy matrix creation
library(ncdf4) # package for netcdf manipulation
library("rnaturalearth")#n_states for map
library("ggrepel") #label_repel on map
library("ggspatial")#annotation_scale in map
library(ggplot2)
library(ggtext)

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
source("R/add_noise.R")
source("R/run_noise_spline.R")
source("R/plot_results.R")
source("R/rate_of_change_fun.R")
source("R/plot_4_sites_results.R")



#---Read in SL data from Andy's data----
SL_data <- readr::read_csv("https://www.dropbox.com/s/43kbaoxqmnrprmv/Common%20Era%20Database%202022_a.csv?dl=1")
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

#---Adding GIA rates from ICE5G for TG-----
SL_df <- add_GIA_rate(SL_df=SL_df, save_loc = "data/SL_df_GIA.csv")
SL_df <- SL_df %>%
  mutate(
    data_lm_slope = ifelse(data_type_id=="TideGaugeData",ICE5_GIA_slope,data_lm_slope),
    data_lm_slope_err = ifelse(data_type_id=="TideGaugeData",0.3,data_lm_slope_err))
  
SL_df <- SL_df %>% mutate(SiteName = as.factor(SiteName))

#--Plot Data---
plot_4_data(SL_df = SL_df,
            save_name = "fig/SL_data_plot_4_sites.pdf",
            n_sites= length(unique(SL_df$SiteName)))
plot_data(SL_df = SL_df,
          save_name = "fig/SL_data_plot_all_sites.pdf",
          n_sites= length(unique(SL_df$SiteName)))
# plot_map(SL_df = SL_df,
#           save_name = "fig/map_data_all_sites.pdf")

#---Spline Basis Functions----
basis_fun_list <- spline_basis_function(SL_df = SL_df)

set_up_option <- " "
save_option <- " "

#--- Stage 1: Run JAGS for model without noise & without Local-----
run_spline(SL_df = SL_df,
           save_option= save_option,
           basis_fun_list)

#---Check convergence----
model_run_1 <- readRDS(paste0("output/model_run_no_local_no_noise", save_option, ".rds"))
plot(model_run_1)

#--- Stage 2: Adding Age uncertainty-----
SL_df <- add_noise(jags.data = jags.data,
                   jags.file = jags.file,
                   model_run = model_run_1,
                   SL_df=SL_df,save_csv = "data/SL_df_noise.csv")

#--- Stage 3: Run JAGS for model with noise-----
run_noise_spline(SL_df = SL_df,
           save_option= save_option,
           basis_fun_list,
           model_run_1 = model_run_1)

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
save_option <- "4_sites"

plot_4_sites_results(SL_df = SL_df,
                     model_run = model_run_2,
                     save_option = save_option,
                     basis_fun_list=basis_fun_list,
                     set_up_option=set_up_option)
 