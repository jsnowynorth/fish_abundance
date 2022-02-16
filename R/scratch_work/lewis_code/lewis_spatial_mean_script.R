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

# X = fish_tmp %>% 
#   filter(COMMON_NAME == 'bluegill') %>% 
#   select(all_of(mean_covs)) %>% 
#   mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
#   mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
#   mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
#   mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
#   mutate(secchi = secchi - mean(secchi)) %>% 
#   as.matrix()
# 
# Y = fish_tmp %>%
#   select(CPUE, COMMON_NAME, DOW, SURVEYDATE, TN) %>%
#   pivot_wider(names_from = 'COMMON_NAME', values_from = 'CPUE') %>% 
#   select(-c(DOW, SURVEYDATE, TN)) %>% 
#   as.matrix()
# 
# t(solve(t(X) %*% X) %*% t(X) %*% log(Y + 0.00001))
# 
# log(apply(Y, 2, mean))


save_pars_spatial_mean = create_pars(fish_tmp, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs)
# save_pars_spatial_mean = read_rds('/home/jntmf/data/fish/new_mod/curr_pars_spatial_mean.rds')
run = sampler(nits = 50000, burnin = 100, thin = 10, check_num = 50, pars = save_pars_spatial_mean)

# saveRDS(run$pars, file = '/Users/joshuanorth/Desktop/constrained_run.rds')
saveRDS(run$pars, file = '/home/jntmf/data/fish/new_mod/curr_pars_spatial_mean.rds')
saveRDS(run, file = '/home/jntmf/data/fish/new_mod/full_model_spatial_mean_1.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/new_mod/full_model_spatial_mean_2.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/new_mod/full_model_spatial_mean_3.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/new_mod/full_model_spatial_mean_4.rds')
# saveRDS(run, file = '/home/jntmf/data/fish/new_mod/full_model_spatial_mean_5.rds')


# nits = 75000
# burnin = 100
# thin = 10
# check_num = 50
# pars = save_pars_spatial_mean

apply(run$beta_accept, c(1,2), mean)

par(mfrow = c(2,3))
for(i in 1:6){
  plot(run$beta_0[i,], type = 'l')
}


png(paste0('/Users/joshuanorth/Desktop/phi_chains/beta_0.png'), width = 800, height = 600)
par(mfrow = c(2,3))
for(i in 1:6){
  plot(run$beta_0[i,], type = 'l')
}
dev.off()


apply(run$beta, c(1,2), mean)

par(mfrow = c(3,3))
for(i in 1:8){
  plot(run$beta[6,i,], type = 'l')
  abline(h = 0)
}

fnames = run$pars$fish_names %>% str_replace(., ' ', '_')
for(j in 1:6){
  # png(paste0('/Users/joshuanorth/Desktop/mean_chains/', fnames[j], '.png'), width = 800, height = 600)
  # png(paste0('/Users/joshuanorth/Desktop/full_chains/', fnames[j], '.png'), width = 800, height = 600)
  png(paste0('/Users/joshuanorth/Desktop/phi_chains/', fnames[j], '_beta.png'), width = 800, height = 600)
  par(mfrow = c(3,3))
  for(i in 1:8){
    plot(run$beta[j,i,], type = 'l')
    abline(h = 0)
  }
  dev.off()
}


par(mfrow = c(3,4))
for(i in 1:11){
  plot(run$phi[6,i,], type = 'l')
}

fnames = run$pars$fish_names %>% str_replace(., ' ', '_')
for(j in 1:6){
  # png(paste0('/Users/joshuanorth/Desktop/mean_chains/', fnames[j], '.png'), width = 800, height = 600)
  # png(paste0('/Users/joshuanorth/Desktop/full_chains/', fnames[j], '.png'), width = 800, height = 600)
  png(paste0('/Users/joshuanorth/Desktop/phi_chains/', fnames[j], '_phi.png'), width = 800, height = 600)
  par(mfrow = c(3,4))
  for(i in 1:11){
    plot(run$phi[j,i,], type = 'l')
  }
  dev.off()
}



par(mfrow = c(3,6))
for(i in 201:203){
  plot(run$omega[i,1,], type = 'l', main = run$pars$fish_names[1])
  plot(run$omega[i,2,], type = 'l', main = run$pars$fish_names[2])
  plot(run$omega[i,3,], type = 'l', main = run$pars$fish_names[3])
  plot(run$omega[i,4,], type = 'l', main = run$pars$fish_names[4])
  plot(run$omega[i,5,], type = 'l', main = run$pars$fish_names[5])
  plot(run$omega[i,6,], type = 'l', main = run$pars$fish_names[6])
}


par(mfrow = c(2,3))
for(i in 1:6){
  plot(run$sigma_species[i,6,], type = 'l')
}
