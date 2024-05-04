####----Noisy Input GAM for RSL along East Coast of North America---------
# 23rd November 2022 with model validations, 10% test set and 90% training set
# Data: Proxies only
# JAGS Set Up:SP Time(Regional) + SP SpaceTime(Local) + [RandomEffect Site & Age(GIA) + RandomEffect Site(y offset)]
# Priors: y offset ,GIA centred on rates from finding the linear regression through data pre industrial revolution

# Running model first without local and then obtaining value for b_r and rerun
#b_r is the coefficients of the regional basis functions 

# Clear workspace
rm(list = ls())

#---------Set working directory--------------
setwd('/users/research/mupton/1. NIGAM_paper_2022/5c.Model_validation_separate_more_sites_proxy_only')
#----------Load packages---------------------
library(tidyverse)
library(geosphere) #distm
library(R2jags)
library(splines)
library(cowplot) # for griding plots at end
library(dbarts)# dummy matrix creation
library(ncdf4) # package for netcdf manipulation
library(dismo) # k fold set up
library(gtools) # ordering the list of files properly

#-----Functions-----
source("R/bs_bbase.R")
source("R/clean_SL_data.R")
source("R/linear_reg_rates.R")
source("R/add_noise.R")
source("R/match.closest.R")
source("R/spline_basis_function.R")
source("R/add_GIA_rate.R")
source("R/run_spline.R")
source("R/run_noise_spline.R")
#source("10_fold_CV.R")

#---Read in SL data from Andy's data----
SL_data <- read_csv("https://www.dropbox.com/s/43kbaoxqmnrprmv/Common%20Era%20Database%202022_a.csv?dl=1")
new_col_names <- c('Basin', 'Region','Site','Reference','Indicator','Latitude','Longitude',
                   'RSL','RSL_er_max','RSL_er_min','Age','Age_2_er_max','Age_2_er_min')
names(SL_data) <- new_col_names

#---Cleaning SL data----
SL_df <- clean_SL_data(SL_df=SL_data,SL_data=SL_data,
                       save_loc="data/SL_df_clean.csv")

#---Linear regression to find rates of data pre 1800-----
lm_data_rates <- linear_reg_rates(SL_df=SL_df)
SL_df <- left_join(SL_df,lm_data_rates, by = "SiteName")

#---Leave out 10% of data at random for test set-----
run_num <- 10
cat(run_num)
set.seed(14561) # full set
N <- nrow(SL_df)

SL_df <- SL_df %>%
  mutate(SiteName = factor(SiteName))# as.factor doing something weird
df_split_index <- kfold(SL_df,k=10,by=SL_df$SiteName)

# Test Set
test_set <- SL_df[df_split_index==run_num,]
write_csv(test_set,paste0("data/test_sets/test_set_",run_num,".csv"))

#Training Set
training_set <- anti_join(SL_df, test_set)
write_csv(training_set,paste0("data/training_sets/training_set_",run_num,".csv"))


# Create the basis functions for training set
training_basis_fun_list <- spline_basis_function(SL_df = training_set,interval = training_set)

# Create the regional basis functions for test set
test_basis_fun_list <- spline_basis_function(SL_df = test_set,interval = training_set)


set_up_option <- paste0("Training set & test set no local run ",run_num)
save_option <- paste0("no_local_training_test_run",run_num)

#--- Stage 1: Run JAGS for model without noise & without Local-----
run_spline(training_set,
           test_set,
           training_basis_fun_list,
           test_basis_fun_list,
           save_option)

model_run_1 <- readRDS(paste0("output/no_noise_no_local/model_run_no_noise_", save_option, ".rds"))
cat("Model 1 finished")
# Check convergence -------------------------------------------------------
plot(model_run_1)

# Plot the regional term
tibble(
  site = training_set$SiteName,
  age = training_set$Age,
  RSL = training_set$RSL,
  regional_mean = model_run_1$BUGSoutput$mean$regional,
  regional_low = apply(model_run_1$BUGSoutput$sims.list$regional,
                       2, 'quantile', 0.05),
  regional_high = apply(model_run_1$BUGSoutput$sims.list$regional,
                        2, 'quantile', 0.95)
) %>% ggplot(aes(x = age, y = regional_mean)) +
  geom_ribbon(aes(ymin=regional_low, ymax=regional_high),
              fill = 'red', alpha = 0.5) +
  geom_line() 


# Stage 2: Adding Age uncertainty---------------------------
training_set <- add_noise(jags.data = jags.data,
                   jags.file = jags.file,
                   model_run = model_run_1,
                   SL_df=training_set,save_csv = paste0("data/training_sets_noise/training_set_noise",run_num,".csv"))

test_set <- add_noise(jags.data = jags.data,
                          jags.file = jags.file,
                          model_run = model_run_1,
                          SL_df=test_set,save_csv = paste0("data/test_sets_noise/test_set_noise",run_num,".csv"))

set_up_option <- paste0("Training set & test set run ",run_num)
save_option <- paste0("training_test_run",run_num)

#--- Stage 3: Run JAGS for model with noise-----
run_noise_spline(training_set = training_set,
                 test_set = test_set,
                 model_run_1 = model_run_1,
                 training_basis_fun_list = training_basis_fun_list,
                 test_basis_fun_list= test_basis_fun_list,
                 save_option = save_option)


cat("Model 2 finished")

# if(run_num==10){
#   
# }



