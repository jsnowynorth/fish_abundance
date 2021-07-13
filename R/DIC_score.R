##----------------------------------------------------
## Name: Joshua North
##
## Date: 10/07/2020
##
## Project: Fish Abundance
##
## Objective: Calculate the DIC score to compare models
##
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
# library(MCMCpack)
library(progress)
library(sparklyr)
library(stringr)
library(readxl)
library(GGally)
library(xtable)
library(Vizumap)
library(latex2exp)
library(ggcorrplot)
library(corrplot)
library(ggpubr)
library(Matrix)
library(spam)
library(abind)
library(zoo)


# load in data ------------------------------------------------------------

fish_dat <- read_csv('data/fish_dat.csv')

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

# model -------------------------------------------------------------------

source('R/lewis_code/lewis_model.R')

# spatial covs
mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:15,17, 25)]
temporal_covs = colnames(fish_dat)[c(23, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]

# load in sampler data ----------------------------------------------------

run_1 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/RSR/spatial/full_model_spatial_2.rds')
run_2 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/RSR/spatial/full_model_spatial_3.rds')
run_3 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/RSR/spatial/full_model_spatial_4.rds')


beta = abind(run_1$beta, 
             run_2$beta,
             run_3$beta, along = 3)

phi = abind(run_1$phi, 
            run_2$phi,
            run_3$phi, along = 3)

omega = abind(run_1$omega, 
              run_2$omega, 
              run_3$omega, along = 2)

sigma_species = abind(run_1$sigma_species, 
                      run_2$sigma_species,
                      run_3$sigma_species, along = 3)

rm(run_1, run_2, run_3)

# pars <- create_pars(fish_dat, mean_covs, mean_covs_log, mean_covs_logit, catch_covs) # non-spatial
pars <- create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs) # spatial
b_names = colnames(pars$X[[1]])
phi_names = colnames(pars$Z[[1]])

b_hat = apply(beta, c(1,2), mean)
phi_hat = apply(phi, c(1,2), mean)
omega_hat = matrix(apply(omega, 1, mean), nrow = pars$n_lakes)
species_hat = apply(sigma_species, c(1,2), mean)


# calculate DIC spatial -----------------------------------------------------------

# # data
# Y = pars$Y
# X = pars$X
# Z = pars$Z
# effort = pars$effort
# beta = b_hat
# phi = phi_hat
# omega = omega_hat
# Sigma_species = species_hat
# lake_id = pars$lake_id
# lake_index = pars$lake_index
# n_lakes = pars$n_lakes

# # parameters
# n = pars$n
# K = pars$K

ll_calc <- function(Y, effort, X, beta, Z, phi, omega, K, lake_id, lake_index){
  
  ind_array = tibble(id = lake_id, as_tibble(omega, .name_repair = ~LETTERS[1:K]))
  lake_array = tibble(id = lake_index)
  OMEGA = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))
  
  ll = 0
  for(k in 1:K){
    
    ll = ll + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + Z[[k]] %*% phi[k,] + OMEGA[,k]), log = T))
  }
  return(ll)
}

D = ll_calc(pars$Y, pars$effort, pars$X, b_hat, pars$Z, phi_hat, omega_hat, pars$K, pars$lake_id, pars$lake_index)

n_post_samps = dim(beta)[3]
D_post <- rep(NA, n_post_samps)

for(i in 1:n_post_samps){
  D_post[i] = ll_calc(pars$Y, pars$effort, pars$X, beta[,,1], pars$Z, phi[,,i], 
                      matrix(omega[,i], nrow = pars$n_lakes), pars$K, pars$lake_id, pars$lake_index)
  
  if(i %% seq(0, n_post_samps, by=1000)){
    print(i)
  }
  
}

plot(D_post, type = 'l')

pD = 0.5 * var(D_post[-c(1:5000)])

DIC_spatial = D + 2*pD

# calculate DIC non-spatial  ------------------------------------------------------

source('R/lewis_code/lewis_model_no_spatial.R')

# non spatial covs
mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:15,17, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]


run_1 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/RSR/non_spatial/')
run_2 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/RSR/spatial/full_model_spatial_3.rds')
run_3 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/RSR/spatial/full_model_spatial_4.rds')


beta_no = abind(run_1$beta, 
             run_2$beta,
             run_3$beta, along = 3)

phi_no = abind(run_1$phi, 
            run_2$phi,
            run_3$phi, along = 3)

omega_no = abind(run_1$omega, 
              run_2$omega, 
              run_3$omega, along = 2)

sigma_species_no = abind(run_1$sigma_species, 
                      run_2$sigma_species,
                      run_3$sigma_species, along = 3)

rm(run_1, run_2, run_3)

# pars <- create_pars(fish_dat, mean_covs, mean_covs_log, mean_covs_logit, catch_covs) # non-spatial
pars_no <- create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs) # spatial
b_names_no = colnames(pars$X[[1]])
phi_names_no = colnames(pars$Z[[1]])

b_hat_no = apply(beta_no, c(1,2), mean)
phi_hat_no = apply(phi_no, c(1,2), mean)
omega_hat_no = matrix(apply(omega_no, 1, mean), nrow = pars$n_lakes)
species_hat_no = apply(sigma_species_no, c(1,2), mean)


ll_calc_no <- function(Y, effort, X, beta, Z, phi, K){
  
  ll = 0
  for(k in 1:K){
    
    ll = ll + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + Z[[k]] %*% phi[k,]), log = T))
  }
  return(ll)
}

D = ll_calc_no(pars_no$Y, pars_no$effort, pars_no$X, b_hat_no, pars_no$Z, phi_hat_no, pars_no$K)

n_post_samps_no = dim(beta_no)[3]
D_post_no <- rep(NA, n_post_samps_no)

for(i in 1:n_post_samps_no){
  D_post[i] = ll_calc_no(pars_no$Y, par_nos$effort, pars_no$X, beta_no[,,1], pars_no$Z, phi_no[,,i], pars_no$K)
  
  if(i %% seq(0, n_post_samps_no, by=1000)){
    print(i)
  }
  
}

plot(D_post_no, type = 'l')

pD_no = 0.5 * var(D_post_no[-c(1:5000)])

DIC_no_spatial = D_no + 2*pD_no