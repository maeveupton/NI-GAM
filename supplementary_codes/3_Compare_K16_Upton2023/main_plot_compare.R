#---- Plotting NIGAM Upton vs K16-----
# Clear workspace
rm(list = ls())

#---------Set working directory--------------
setwd('/users/research/mupton/1. NIGAM_paper_2022/3. Compare_K16_Upton23')
#----------Load packages---------------------
library(tidyverse)
#------------SL data-----------------
SL_df <- read_csv("data/SL_df_GIA.csv")
# Match the site names
SL_df$SiteName <- gsub("\\,.*","",SL_df$SiteName)
# #------Ordering the sites-------
# all_data_sites<-factor(paste0(SL_df$Longitude, SL_df$Longitude), labels = 1:n_sites)#For ordering sites
# SL_df <- cbind(SL_df,all_data_sites=all_data_sites)
# order_sites <- SL_df %>%  group_by(SiteName,all_data_sites) %>%
#   dplyr::summarise(n = n()) %>%
#   arrange(all_data_sites)
# SL_df$SiteName <- factor(SL_df$SiteName,levels = unique(order_sites$SiteName))

#erica_file <- "ice5g_1850"
erica_file <- "linear_zero"

#---Read in K16 results for full model fit----
#K16_full_model_fit <- read.table(file = "data/RSLsFull.tsv",sep = '\t',check.names = FALSE)#10 year pred
#K16_full_model_fit <- read.table(file = "data/50_year_pred/RSLsFull.tsv",sep = '\t',check.names = FALSE)#50 year pred
K16_full_model_fit <- read.table(file = paste0("data/",erica_file,"/RSLsFull.tsv"),sep = '\t',check.names = FALSE)#GIA update
names(K16_full_model_fit) <- K16_full_model_fit[1,]
K16_full_model_fit <- K16_full_model_fit[-1,]
#--- Rotate the data----
K16_full_model_fit_long_df <- K16_full_model_fit %>% 
  pivot_longer(!c(Site,ID,lat,long,"RSL_class"),
               names_to = c("Year"),
               values_to = c("RSL"))

K16_full_model_fit_wide_df<- K16_full_model_fit_long_df %>% 
  pivot_wider(names_from = RSL_class, values_from = RSL)
#----Clean data----
K16_full_model_fit_df <- K16_full_model_fit_wide_df %>%
  rename(RSL=RSLs, sd_RSL = SDs, 
         Longitude=long,Latitude= lat,SiteName = Site) %>% 
  mutate(RSL = RSL/1000,#m
         sd_RSL = sd_RSL/1000,
         SiteName = as.factor(SiteName),
         Year = as.numeric(gsub("X.",replacement = '',Year))/1000,
         ID = "K16")
K16_full_model_fit_df <- K16_full_model_fit_df %>% 
  group_by(SiteName) %>% 
  mutate(RSL_er_up = RSL + qnorm(0.975)*sd_RSL,
         RSL_er_lwr = RSL - qnorm(0.975)*sd_RSL)

#--- Upton Full Model Fit data----
NIGAM_full_model <- read.csv("data/Upton_result/total_model_fit.csv")
NIGAM_full_model$ID <- "This study"
# Match the site names
NIGAM_full_model$SiteName <- gsub("\\,.*","",NIGAM_full_model$SiteName)

#---Full Model Fit Plot----
full_model_fit_plot <- ggplot()+
  geom_line(data =K16_full_model_fit_df, aes(x = Year*1000, 
                                             y = RSL,colour = ID))+
  geom_ribbon(data = K16_full_model_fit_df, aes(x = Year*1000,
                                        ymin = RSL_er_lwr,ymax=RSL_er_up,colour=ID),alpha = 0.2)+
  geom_line(data =NIGAM_full_model, aes(x = Age*1000, 
                                             y = RSL,colour = ID))+
  geom_ribbon(data = NIGAM_full_model, aes(x = Age*1000,
                                                ymin = lwr,ymax=upr,colour=ID),alpha = 0.2)+
  facet_wrap(~SiteName)+
  ggtitle("Total")+
  theme_bw()
full_model_fit_plot
#ggsave(full_model_fit_plot, file = "fig/total_model_compare.pdf", width = 10, height = 6)
#ggsave(full_model_fit_plot, file = "fig//50_year_pred/total_model_compare.pdf", width = 10, height = 6)
ggsave(full_model_fit_plot, file = paste0("fig/",erica_file,"/total_model_compare.pdf"), width = 10, height = 6)

#---- Read in K16 results for Common Regional-----
#K16_common_reg <- read.table(file = "data/RSLsCommon Regional.tsv",sep = '\t',check.names = FALSE)# 10 year pred
#K16_common_reg <- read.table(file = "data/50_year_pred/RSLsCommon Regional.tsv",sep = '\t',check.names = FALSE)# 50 year Pred
K16_common_reg <- read.table(file = paste0("data/",erica_file,"/RSLsCommon Regional.tsv"),sep = '\t',check.names = FALSE)# 50 year Pred
names(K16_common_reg) <- K16_common_reg[1,]
K16_common_reg <- K16_common_reg[-1,]
#--- Rotate the data----
K16_common_reg_long_df <- K16_common_reg %>% 
  pivot_longer(!c(Site,ID,lat,long,"RSL_class"),
               names_to = c("Year"),
               values_to = c("RSL"))
                                                     
K16_common_reg_wide_df<- K16_common_reg_long_df %>% 
  pivot_wider(names_from = RSL_class, values_from = RSL)
#----Clean data----
K16_common_df <- K16_common_reg_wide_df %>%
  rename(RSL=RSLs, sd_RSL = SDs, 
         Longitude=long,Latitude= lat,SiteName = Site) %>% 
  mutate(RSL = RSL/1000,#m
         sd_RSL = sd_RSL/1000,
         SiteName = as.factor(SiteName),
         Year = as.numeric(gsub("X.",replacement = '',Year))/1000,
         ID = "K16")
K16_common_df <- K16_common_df %>% 
  group_by(SiteName) %>% 
  mutate(RSL_er_up = RSL + qnorm(0.975)*sd_RSL,
         RSL_er_lwr = RSL - qnorm(0.975)*sd_RSL)
#--- Upton Regional Fit data----
NIGAM_regional <- read.csv("data/Upton_result/regional_grid_model_fit.csv")
NIGAM_regional$ID <- "This study"
# Match the site names
#NIGAM_regional$SiteName <- gsub("\\,.*","",NIGAM_regional$SiteName)

#---Common Regional Plot----
com_reg_plot <- ggplot()+
  geom_line(data =K16_common_df, aes(x = Year*1000, y = RSL*100,colour = ID))+
  geom_ribbon(data = K16_common_df, aes(x = Year*1000,
                                        ymin = RSL_er_lwr*100,ymax=RSL_er_up*100),alpha = 0.2)+
  geom_line(data =NIGAM_regional, aes(x = Age*1000, y = RSL*100,colour = ID))+
  geom_ribbon(data = NIGAM_regional, aes(x = Age*1000,
                                        ymin = lwr*100,ymax=upr*100,colour = ID),alpha = 0.2)+
  ggtitle("Regional(?)")+
  #facet_wrap(~SiteName)+
  theme_bw()
com_reg_plot
ggsave(com_reg_plot, file = paste0("fig/",erica_file,"/regional_compare.pdf"), width = 10, height = 6)
#ggsave(com_reg_plot, file = "fig/50_year_pred/regional_compare.pdf", width = 10, height = 6)
#ggsave(com_reg_plot, file = "fig/regional_compare.pdf", width = 10, height = 6)


#---- Read in K16 results for linear - GIA-----
#K16_lin_GIA <- read.table(file = "data/RSLsLinear-GIA.tsv",sep = '\t',check.names = FALSE)
#K16_lin_GIA <- read.table(file = "data/50_year_pred/RSLsLinear-GIA.tsv",sep = '\t',check.names = FALSE)
#K16_lin_GIA <- read.table(file = "data/ice5_refyr0/RSLsLinear-GIA.tsv",sep = '\t',check.names = FALSE)
K16_lin_GIA <- read.table(file = paste0("data/",erica_file,"/RSLsLinear-GIA.tsv"),sep = '\t',check.names = FALSE)
names(K16_lin_GIA) <- K16_lin_GIA[1,]
K16_lin_GIA <- K16_lin_GIA[-1,]
#--- Rotate the data----
K16_lin_GIA_long_df <- K16_lin_GIA %>% 
  pivot_longer(!c(Site,ID,lat,long,"RSL_class"),
               names_to = c("Year"),
               values_to = c("RSL"))

K16_lin_GIA_wide_df<- K16_lin_GIA_long_df %>% 
  pivot_wider(names_from = RSL_class, values_from = RSL)
#----Clean data----
K16_lin_GIA_df <- K16_lin_GIA_wide_df %>%
  rename(RSL=RSLs, sd_RSL = SDs, 
         Longitude=long,Latitude= lat,SiteName = Site) %>% 
  mutate(RSL = RSL/1000,#m
         sd_RSL = sd_RSL/1000,
         SiteName = as.factor(SiteName),
         Year = as.numeric(gsub("X.",replacement = '',Year))/1000,
         ID = "K16")
K16_lin_GIA_df <- K16_lin_GIA_df %>% 
  group_by(SiteName) %>% 
  mutate(RSL_er_up = RSL + qnorm(0.975)*sd_RSL,
         RSL_er_lwr = RSL - qnorm(0.975)*sd_RSL)

#--- Upton Linear GIA Fit data----
# CHECK This one
NIGAM_GIA <- read.csv("data/Upton_result/GIA_model_fit.csv")
NIGAM_GIA$ID <- "This study"
# Match the site names
NIGAM_GIA$SiteName <- gsub("\\,.*","",NIGAM_GIA$SiteName)

#---Linear GIA Plot----
lin_GIA_plot <- ggplot()+
  geom_line(data =K16_lin_GIA_df, aes(x = Year*1000, y = RSL,colour = ID))+
  geom_ribbon(data = K16_lin_GIA_df, aes(x = Year*1000,
                                        ymin = RSL_er_lwr,ymax=RSL_er_up),alpha = 0.2)+
  geom_line(data =NIGAM_GIA, aes(x = Age*1000, y = RSL,colour = ID))+
  geom_ribbon(data = NIGAM_GIA, aes(x = Age*1000,
                                         ymin = lwr,ymax=upr),alpha = 0.2)+
  
  facet_wrap(~SiteName)+
  ggtitle("GIA")+
  theme_bw()
lin_GIA_plot
#ggsave(lin_GIA_plot, file = "fig/GIA_compare.pdf", width = 10, height = 6)
#ggsave(lin_GIA_plot, file = "fig/50_year_pred/GIA_compare.pdf", width = 10, height = 6)
ggsave(lin_GIA_plot, file = paste0("fig/",erica_file,"/GIA_compare.pdf"), width = 10, height = 6)

#---- Read in K16 results for Non-linear Local-----
#K16_local <- read.table(file = "data/RSLsLocal Non-linear.tsv",sep = '\t',check.names = FALSE)
#K16_local <- read.table(file = "data/50_year_pred/RSLsLocal Non-linear.tsv",sep = '\t',check.names = FALSE)
K16_local <- read.table(file = paste0("data/",erica_file,"/RSLsLocal Non-linear.tsv"),sep = '\t',check.names = FALSE)
names(K16_local) <- K16_local[1,]
K16_local <- K16_local[-1,]
#--- Rotate the data----
K16_local_long_df <- K16_local %>% 
  pivot_longer(!c(Site,ID,lat,long,"RSL_class"),
               names_to = c("Year"),
               values_to = c("RSL"))

K16_local_wide_df<- K16_local_long_df %>% 
  pivot_wider(names_from = RSL_class, values_from = RSL)
#----Clean data----
K16_local_df <- K16_local_wide_df %>%
  rename(RSL=RSLs, sd_RSL = SDs, 
         Longitude=long,Latitude= lat,SiteName = Site) %>% 
  mutate(RSL = RSL/1000,#m
         sd_RSL = sd_RSL/1000,
         SiteName = as.factor(SiteName),
         Year = as.numeric(gsub("X.",replacement = '',Year))/1000,
         ID = "K16")
K16_local_df <- K16_local_df %>% 
  group_by(SiteName) %>% 
  mutate(RSL_er_up = RSL + qnorm(0.975)*sd_RSL,
         RSL_er_lwr = RSL - qnorm(0.975)*sd_RSL)

#--- Upton Local Fit data----
NIGAM_local <- read.csv("data/Upton_result/local_model_fit.csv")
NIGAM_local$ID <- "This study"
# Match the site names
NIGAM_local$SiteName <- gsub("\\,.*","",NIGAM_local$SiteName)

#---Local Plot----
local_plot <- ggplot()+
  geom_line(data =K16_local_df, aes(x = Year*1000, y = RSL*100,colour = ID))+
  geom_ribbon(data = K16_local_df, aes(x = Year*1000,
                                         ymin = RSL_er_lwr*100,ymax=RSL_er_up*100),alpha = 0.2)+
  geom_line(data =NIGAM_local, aes(x = Age*1000, y = RSL*100,colour = ID))+
  geom_ribbon(data = NIGAM_local, aes(x = Age*1000,
                                       ymin = lwr*100,ymax=upr*100),alpha = 0.2)+
  facet_wrap(~SiteName)+
  ggtitle("Local")+
  theme_bw()
local_plot
#ggsave(local_plot, file = "fig/local_compare.pdf", width = 10, height = 6)
#ggsave(local_plot, file = "fig/50_year_pred/local_compare.pdf", width = 10, height = 6)
ggsave(local_plot, file = paste("fig/",erica_file,"/local_compare.pdf"), width = 10, height = 6)
