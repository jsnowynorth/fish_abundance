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
source('/home/jntmf/data/fish/code/lewis_model.R')

# load data ---------------------------------------------------------------
fish_dat = read_csv('/home/jntmf/data/fish/data/fish_dat.csv')

fish_dat = fish_dat %>% 
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME))

mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:15,17, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]

# run model ---------------------------------------------------------------

save_pars_spatial = create_pars(fish_dat, mean_covs, mean_covs_log, mean_covs_logit, catch_covs)
# save_pars_spatial = read_rds('/home/jntmf/data/fish/results/spatial/curr_pars_spatial.rds')
run = sampler(nits = 50000, burnin = 1000, thin = 10, check_num = 100, pars = save_pars_spatial)

saveRDS(run$pars, file = '/home/jntmf/data/fish/results/spatial/curr_pars_spatial.rds')
saveRDS(run, file = '/home/jntmf/data/fish/results/spatial/full_model_spatial_1.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/results/spatial/full_model_spatial_2.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/results/spatial/full_model_spatial_3.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/results/spatial/full_model_spatial_4.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/results/spatial/full_model_spatial_5.rds')

