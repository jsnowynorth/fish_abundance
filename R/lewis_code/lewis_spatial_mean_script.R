##----------------------------------------------------
## Name: Joshua North
##
## Date: 03/01/2020
##
## Project: Fish Abundance
##
## Objective: Multivariate fish abundance model. Land use covariates, GDD, secchi.
##            Separate function for each component.
## Notes:
##
##----------------------------------------------------



# load libraries ----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(readxl)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(fields)
library(viridis)
library(mvtnorm)
library(progress)
library(sparklyr)
library(stringr)
library(GGally)
library(Matrix)
library(spam)


# source functions --------------------------------------------------------
# source('/home/jntmf/data/fish/new_mod/lewis_model_mean.R')
source('R/lewis_code/lewis_model_mean.R')

# load data ---------------------------------------------------------------
# fish_dat = read_csv('/home/jntmf/data/fish/new_mod/fish_dat.csv')
fish_dat = read_csv('data/fish_dat.csv')

fish_dat = fish_dat %>% 
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME)) %>% 
  arrange(DOW, year, COMMON_NAME)

mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:15,17, 25)]
temporal_covs = colnames(fish_dat)[c(23, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]

# run model ---------------------------------------------------------------
set.seed(1)

# randomly sample
fish_tmp = fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  distinct() %>% 
  slice_sample(n = 500, replace = F) %>% 
  left_join(fish_dat) %>% 
  mutate_at(vars(DOW), ~droplevels(.))

# length(unique(fish_tmp$DOW))

# 308 lakes
# only lakes observed 4 or more times
fish_tmp = fish_dat %>% 
  group_by(DOW) %>% 
  filter(n() >= 4*12) %>% 
  ungroup() %>% 
  mutate_at(vars(DOW), ~droplevels(.))

# length(unique(fish_tmp$DOW))



save_pars_spatial_mean = create_pars(fish_tmp, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs)
# save_pars_spatial_mean = read_rds('/home/jntmf/data/fish/new_mod/curr_pars_spatial_mean.rds')
run = sampler(nits = 100000, burnin = 200, thin = 10, check_num = 50, pars = save_pars_spatial_mean)

saveRDS(run$pars, file = '/home/jntmf/data/fish/new_mod/curr_pars_spatial_mean.rds')
saveRDS(run, file = '/home/jntmf/data/fish/new_mod/full_model_spatial_mean_1.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/new_mod/full_model_spatial_mean_2.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/new_mod/full_model_spatial_mean_3.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/new_mod/full_model_spatial_mean_4.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/new_mod/full_model_spatial_mean_5.rds')


par(mfrow = c(3,3))
for(i in 1:9){
  plot(run$beta[6,i,], type = 'l')
  abline(h = 0)
}


par(mfrow = c(3,4))
for(i in 1:11){
  plot(run$phi[4,i,], type = 'l')
}