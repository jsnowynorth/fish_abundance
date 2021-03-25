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





# load in data ------------------------------------------------------------

static = read_csv('data/Static_MN_Data_JE.csv') %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
effort = read_csv('data/MN_fish_catch_effort_all_lakes_years_surveys_JE.csv') %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
temp = read_rds('data/daily_degree_days_MN_lakes_glm2.rds') %>% ungroup()

temp = temp %>%
  select(date, temp_0, MNDOW_ID, DD5, DOY) %>% # C5 temperature
  rename(SURVEYDATE = date) %>%
  mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>%
  select(-MNDOW_ID)

GDD = temp %>% 
  mutate(year = year(SURVEYDATE)) %>% 
  group_by(year, DOW) %>% 
  summarise(DD5 = max(DD5)) %>% 
  ungroup()


secchi_year = read_csv('data/annual_median_remote_sensing_secchi.csv') %>% 
  mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0"))) %>% 
  select(year, annual.median.rs, DOW) %>% 
  rename('secchi' = 'annual.median.rs')

# secchi_year = read_csv('data/annual_median_secchi_KV_model.csv') %>% 
#   mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0"))) %>% 
#   select(DOW:median.secchi.m) %>% 
#   rename('secchi' = 'median.secchi.m',
#          'year' = 'Year')
# 
# secchi_year_2 = read_csv('data/Secchi_annual_CDOM_lakedepth.csv') %>% 
#   mutate(DOW = (str_pad(DOWLKNUM, 8, side = "left", pad = "0"))) %>% 
#   select(DOW, Year, avgSecchi.m) %>% 
#   rename('secchi' = 'avgSecchi.m',
#          'year' = 'Year')
# 
# secchi = secchi_year %>% 
#   full_join(secchi_year_2, by = c('DOW', 'year', 'secchi')) %>% 
#   group_by(year) %>% 
#   distinct(DOW, .keep_all = T) %>% 
#   ungroup() %>% 
#   arrange(DOW, year)
# 
# rm(secchi_year, secchi_year_2)
  

land = read_xlsx('data/MN_lake_landuse.xlsx')

land = land %>%
  select(MNDOW_ID, year:grass, lake_elevation_m) %>% 
  mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>% 
  select(-MNDOW_ID) %>% 
  filter(complete.cases(.)) %>% 
  mutate(DOW = as.factor(DOW)) %>% 
  distinct(DOW, .keep_all = T)


# data cleaning -----------------------------------------------------------

# combine survey data with some static variables
fish_dat = effort %>% 
  left_join(static, by = 'DOW') %>%
  select(DOW, COMMON_NAME, TOTAL_CATCH, EFFORT, GEAR, CPUE, SURVEYDATE, MAX_DEPTH_FEET, mean.gdd, LAKE_AREA_GIS_ACRES, 
         LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>%
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME),
         GEAR = as.factor(GEAR),
         SURVEYDATE = mdy(SURVEYDATE)) %>%
  filter(complete.cases(.)) %>% 
  filter(COMMON_NAME != 'white sucker',
         COMMON_NAME != 'smallmouth bass') %>%
  filter(SURVEYDATE >= '1993-01-01') %>% 
  filter(GEAR == 'GN' | GEAR == 'TN') %>%
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW),
         GEAR = droplevels(GEAR)) %>% 
  arrange(DOW) %>% 
  mutate(year = year(SURVEYDATE))

# add in land use
fish_dat = fish_dat %>% 
  left_join(land %>% select(-year), by = c('DOW')) %>% 
  filter(complete.cases(.))

# all <- fish_dat %>% 
#   group_by(SURVEYDATE, DOW) %>% 
#   tidyr::expand(COMMON_NAME, SURVEYDATE, GEAR) %>% 
#   ungroup() %>% 
#   arrange(DOW)

# all combos of fish, survey date, and gear
all <- fish_dat %>% 
  group_by(DOW) %>% 
  tidyr::expand(COMMON_NAME, SURVEYDATE, GEAR) %>% 
  ungroup() %>% 
  arrange(DOW)


fish_dat <- fish_dat %>% 
  right_join(all) %>%
  mutate(EFFORT = coalesce(EFFORT, 0L),
         TOTAL_CATCH = coalesce(TOTAL_CATCH, 0L),
         CPUE = coalesce(CPUE, 0L)) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:lake_elevation_m), list(~coalesce(., 0L))) %>% 
  group_by(GEAR, SURVEYDATE, DOW) %>% 
  mutate(EFFORT = max(EFFORT)) %>%
  ungroup() %>% 
  group_by(DOW, SURVEYDATE) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:lake_elevation_m), ~ mean(.[!is.na(.) & . != 0])) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:lake_elevation_m), ~ ifelse(is.na(.), 0, .)) %>% 
  ungroup() %>% 
  mutate(TN = ifelse(GEAR == 'TN', 1, 0),
         GN = ifelse(GEAR == 'GN', 1, 0)) %>% 
  select(-GEAR) %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW)) %>%
  mutate(CPUE = TOTAL_CATCH/EFFORT,
         CPUE = ifelse(is.na(CPUE), 0, CPUE)) %>% 
  arrange(DOW)


# add temperature 

fish_dat = fish_dat %>% 
  inner_join(GDD %>% 
               mutate(DD5 = (DD5 - mean(DD5))/sd(DD5))) %>% 
  inner_join(temp %>% 
               select(SURVEYDATE, temp_0, DOW) %>% 
               mutate(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0))) %>% 
  inner_join(secchi_year, by = c('DOW', 'year'))


tmp = fish_dat %>% 
  inner_join(GDD %>% 
               mutate(DD5 = (DD5 - mean(DD5))/sd(DD5))) %>% 
  inner_join(temp %>% 
               select(SURVEYDATE, temp_0, DOW)) %>% 
  inner_join(secchi_year, by = c('DOW', 'year'))

# fish_dat = fish_dat %>% 
#   inner_join(GDD) %>% 
#   mutate(DD5 = (DD5 - mean(DD5))/sd(DD5))

# add in DOY info

# fish_dat = fish_dat %>% 
#   mutate(DOY = yday(SURVEYDATE),
#          DOY_sin = sin(DOY/365 * 2*pi),
#          DOY_cos = cos(DOY/365 * 2*pi),
#          DOY_sin_semi = sin(DOY/365 * 4*pi),
#          DOY_cos_semi = cos(DOY/365 * 4*pi))

fish_dat = fish_dat %>% 
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/365 * 4*pi),
         DOY_cos_semi = cos(DOY/365 * 4*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0)


colnames(fish_dat)
mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:19, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9, 13:19)]
catch_covs = colnames(fish_dat)[c(24, 27, 28)] # no temp doy interaction
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]
# need columns: DOW, COMMON_NAME, EFFORT, SURVEYDATE

# ggplot(fish_dat, aes(x = MAX_DEPTH_FEET, y = secchi)) +
#   geom_point()

# fish_dat %>%
#   group_by(year) %>%
#   distinct(DOW, .keep_all = T) %>%
#   select(all_of(mean_covs)) %>%
#   ungroup() %>%
#   mutate_at(vars(mean_covs_log), ~ ifelse(. == 0, . + 0.001, .)) %>%
#   mutate_at(vars(mean_covs_log), ~ log(.)) %>%
#   ggpairs()
# 
# 
# fish_dat %>%
#   group_by(year) %>%
#   distinct(DOW, .keep_all = T) %>%
#   select(all_of(mean_covs), all_of(catch_covs)) %>%
#   ungroup() %>%
#   mutate_at(vars(all_of(mean_covs_log)), ~ ifelse(. == 0, . + 0.001, .)) %>%
#   mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>%
#   ggpairs()


# model -------------------------------------------------------------------

create_pars <- function(fish_dat, mean_covs, mean_covs_log, catch_covs){
  
  K = nlevels(fish_dat$COMMON_NAME)
  levs = levels(fish_dat$COMMON_NAME)
  lake_index = (fish_dat %>% filter(COMMON_NAME == levs[1]))$DOW # lake_index == lake_id[2]
  lake_id = levels(fish_dat$DOW)
  n_lakes = length(levels(fish_dat$DOW))
  
  pars = list()
  
  # data
  pars$Y = list()
  pars$X = list()
  pars$Z = list()
  pars$effort = list()
  
  for(k in 1:K){
    pars$Y[[k]] = (fish_dat %>% filter(COMMON_NAME == levs[k]) %>% select(TOTAL_CATCH))$TOTAL_CATCH
  }
  
  for(k in 1:K){
    X = fish_dat %>% 
      filter(COMMON_NAME == levs[k]) %>% 
      select(all_of(mean_covs)) %>% 
      mutate_at(vars(all_of(mean_covs_log)), ~ ifelse(. == 0, . + 0.001, .)) %>% 
      mutate_at(vars(all_of(mean_covs_log)), ~ log(.))
    
    Z = fish_dat %>% 
      filter(COMMON_NAME == levs[k]) %>%
      select(all_of(catch_covs), GN) %>% 
      mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN))
    
    pars$X[[k]] = as.matrix(X)
    pars$Z[[k]] = as.matrix(Z)
  }
  
  for(k in 1:K){
    pars$effort[[k]] = (fish_dat %>% filter(COMMON_NAME == levs[k]) %>% select(EFFORT))$EFFORT
  }
  
  
  # parameters
  pars$n = unlist(lapply(pars$Y, length))
  pars$p_beta = ncol(pars$X[[1]])
  pars$p_phi = ncol(pars$Z[[1]])
  pars$K = K
  
  pars$beta_0 = 0
  pars$beta_0_accept = 0
  pars$beta_0_prior_var = 100
  pars$beta = array(0, dim = c(K, pars$p_beta))
  pars$beta_accept =  array(0, dim = c(K, pars$p_beta))
  pars$beta_prior_var = 100
  pars$phi = array(0, dim = c(K, pars$p_phi))
  pars$phi_accept =  array(0, dim = c(K, pars$p_phi))
  pars$phi_prior_var = 100
  
  # hyperpriors
  pars$Sigma_species = diag(K)
  pars$nu_species = K + 10
  pars$Psi_species = diag(K)
  
  pars$omega = matrix(rep(0, K), ncol = K)
  
  # Proposal variances
  pars$sig_prop_beta = array(2, dim = c(K, pars$p_beta))
  pars$sig_prop_phi = array(2, dim = c(K, pars$p_phi))
  pars$sig_prop_beta_0 = 2
  
  # indexing
  pars$lake_index = lake_index
  pars$lake_id = lake_id
  pars$n_lakes = n_lakes
  pars$fish_names = levels(fish_dat$COMMON_NAME)
  
  return(pars)
  
}

update_beta <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  Z = pars$Z
  effort = pars$effort
  omega = pars$omega
  phi = pars$phi
  
  # parameters
  n = pars$n
  p_beta = pars$p_beta
  p_phi = pars$p_phi
  K = pars$K
  
  # beta monitor values
  beta_accept = array(0, dim = c(K, p_beta))
  beta_curr = pars$beta
  sig_prop_beta = pars$sig_prop_beta
  beta_prior_var = pars$beta_prior_var
  
  # beta_0 monitor
  beta_0 = pars$beta_0
  beta_0_accept = 0
  beta_0_prior_var = pars$beta_0_prior_var
  sig_prop_beta_0 = pars$sig_prop_beta_0
  
  # beta_0
  
  b_prop = rnorm(1, beta_0, sig_prop_beta_0)
  
  
  like_curr = like_prop = 0
  for(k in 1:K){
    like_curr = like_curr + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta_curr[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T)) + dnorm(beta_0, 0, sig_prop_beta_0, log = T)
    like_prop = like_prop + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(b_prop + X[[k]] %*% beta_curr[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T)) + dnorm(b_prop, 0, sig_prop_beta_0, log = T)
  }
  
  like_curr = like_prop = 0
  for(k in 1:K){
    like_curr = like_curr + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta_curr[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T)) + dnorm(beta_0, 0, sig_prop_beta_0, log = T)
    like_prop = like_prop + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(b_prop + X[[k]] %*% beta_curr[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T)) + dnorm(b_prop, 0, sig_prop_beta_0, log = T)
  }
  # if(is.na(like_prop) | is.na(like_curr)){print(b_prop, beta_0, beta_curr)}
  if((like_prop - like_curr) > log(runif(1))){
    beta_0 = b_prop
    beta_0_accept = 1
  }
  
  # all other betas
  for(i in 1:p_beta){
    for(k in 1:K){
      
      beta_prop = beta_curr
      
      b_prop = beta_prop[k,i] = rnorm(1, beta_curr[k,i], sig_prop_beta[k,i])
      b_curr = beta_curr[k,i]
      
      like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta_curr[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T)) + dnorm(b_curr, 0, beta_prior_var, log = T)
      like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta_prop[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T)) + dnorm(b_prop, 0, beta_prior_var, log = T)
      
      # if(is.na(like_prop) | is.na(like_curr)){print(b_prop, beta_0, beta_curr)}
      if((like_prop - like_curr) > log(runif(1))){
        beta_curr[k,i] = b_prop
        beta_accept[k,i] = 1
      }
      
    }
    
  }
  
  pars$beta_0 = beta_0
  pars$beta_0_accept = beta_0_accept
  pars$beta = beta_curr
  pars$beta_accept = beta_accept
  
  return(pars)
  
  
}

update_phi <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  Z = pars$Z
  effort = pars$effort
  omega = pars$omega
  beta = pars$beta
  beta_0 = pars$beta_0
  
  # parameters
  n = pars$n
  p_beta = pars$p_beta
  p_phi = pars$p_phi
  K = pars$K
  
  # beta monitor values
  phi_accept = array(0, dim = c(K, p_phi))
  phi_curr = pars$phi
  sig_prop_phi = pars$sig_prop_phi
  phi_prior_var = pars$phi_prior_var
  
  
  # update_phi
  for(i in 1:p_phi){
    for(k in 1:K){
      
      phi_prop = phi_curr
      
      p_prop = phi_prop[k,i] = rnorm(1, phi_curr[k,i], sig_prop_phi[k,i])
      p_curr = phi_curr[k,i]
      
      like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta[k,] + Z[[k]] %*% phi_curr[k,] + omega[k]), log = T)) + dnorm(p_curr, 0, phi_prior_var, log = T)
      like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta[k,] + Z[[k]] %*% phi_prop[k,] + omega[k]), log = T)) + dnorm(p_prop, 0, phi_prior_var, log = T)
      
      # if(is.na(like_prop) | is.na(like_curr)){print(b_prop, beta_0, beta_curr)}
      if((like_prop - like_curr) > log(runif(1))){
        phi_curr[k,i] = p_prop
        phi_accept[k,i] = 1
      }
      
    }
    
  }
  
  pars$phi = phi_curr
  pars$phi_accept = phi_accept
  
  return(pars)
  
  
}

ll_calc <- function(Y, effort, beta_0, X, beta, Z, phi, omega, K){
  ll = 0
  for(k in 1:K){
    ll = ll + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T))
  }
  return(ll)
}

update_omega <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  Z = pars$Z
  effort = pars$effort
  beta_0 = pars$beta_0
  beta = pars$beta
  phi = pars$phi
  omega = pars$omega
  Sigma_species = pars$Sigma_species
  
  # parameters
  n = pars$n
  p = pars$p
  K = pars$K
  
  # choose ellipse
  nu = rmvnorm(1, rep(0, K), Sigma_species)
  nu = nu - mean(nu)
  
  # ll threshold
  logy = ll_calc(Y, effort, beta_0, X, beta, Z, phi, omega, K) + log(runif(1))
  
  # draw initial proposal
  theta = runif(1, 0, 2*pi)
  theta_min = theta - 2*pi
  theta_max = theta
  
  f_proposal = omega * cos(theta) + nu * sin(theta)
  logf = ll_calc(Y, effort, beta_0, X, beta, Z, phi, f_proposal, K)
  
  if(logf > logy){
    keeper = f_proposal
  }else{
    
    while(logf < logy){
      if(theta < 0){
        theta_min = theta
      }else{
        theta_max = theta
      }
      theta = runif(1, theta_min, theta_max)
      f_proposal = omega * cos(theta) + nu * sin(theta)
      logf = ll_calc(Y, effort, beta_0, X, beta, Z, phi, f_proposal, K)
    }
    
    keeper = f_proposal
    
  }
  
  
  pars$omega = keeper
  
  return(pars)
  
  
}

update_sigma_species <- function(pars){
  
  # data
  Sigma_species = pars$Sigma_species
  omega = pars$omega
  
  nu_species = pars$nu_species
  Psi_species = pars$Psi_species
  
  # parameters
  K = pars$K
  
  nu_hat = nu_species + K
  psi_hat = Psi_species + t(omega) %*% omega
  
  pars$Sigma_species = MCMCpack::riwish(nu_hat, psi_hat)
  
  return(pars)
  
}

update_proposal_var_beta <- function(pars, beta_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop_beta
  
  bp = beta_accept_post[,,(i-check_num+1):i]
  accept_rate = apply(bp, c(1,2), mean)
  
  sig_prop = ifelse(accept_rate < 0.2, sig_prop*0.9, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.9, sig_prop)
  
  pars$sig_prop_beta = sig_prop
  
  return(pars)
  
  
}

update_proposal_var_phi <- function(pars, phi_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop_phi
  
  bp = phi_accept_post[,,(i-check_num+1):i]
  accept_rate = apply(bp, c(1,2), mean)
  
  sig_prop = ifelse(accept_rate < 0.2, sig_prop*0.9, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.9, sig_prop)
  
  pars$sig_prop_phi = sig_prop
  
  return(pars)
  
  
}

update_proposal_var_beta_0 <- function(pars, beta_0_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop_beta_0
  
  accept_rate = mean(beta_0_accept_post[(i-check_num+1):i])
  
  if(accept_rate < 0.2){
    sig_prop = sig_prop*0.9
  }else if(accept_rate > 0.45){
    sig_prop = sig_prop/0.9
  }else{
    sig_prop = sig_prop
  }
  
  pars$sig_prop_beta_0 = sig_prop
  
  return(pars)
  
  
}

sampler <- function(fish_dat, mean_covs, mean_covs_log, catch_covs, nits, check_num = 200){
  
  pars <- create_pars(fish_dat, mean_covs, mean_covs_log, catch_covs)
  
  p_beta = pars$p_beta
  p_phi = pars$p_phi
  K = pars$K
  
  beta_post = array(NA, dim = c(K, p_beta, nits))
  beta_accept_post = array(NA, dim = c(K, p_beta, nits))
  phi_post = array(NA, dim = c(K, p_phi, nits))
  phi_accept_post = array(NA, dim = c(K, p_phi, nits))
  beta_0_post = array(NA, dim = c(nits))
  beta_0_accept_post = array(NA, dim = c(nits))
  omega_post = array(NA, dim = c(length(pars$omega), nits))
  sigma_species_post = array(NA, dim = c(dim(pars$Sigma_species), nits))
  
  pb <- progress_bar$new(
    format = "  Running [:bar] :percent eta: :eta",
    total = nits, clear = FALSE, width = 60)
  
  for(i in seq(1, nits)){
    
    pars <- update_beta(pars)
    pars <- update_phi(pars)
    pars <- update_omega(pars)
    pars <- update_sigma_species(pars)
    
    beta_post[,,i] = pars$beta
    beta_accept_post[,,i] = pars$beta_accept
    phi_post[,,i] = pars$phi
    phi_accept_post[,,i] = pars$phi_accept
    beta_0_post[i] = pars$beta_0
    beta_0_accept_post[i] = pars$beta_0_accept
    omega_post[,i] = pars$omega
    sigma_species_post[,,i] = pars$Sigma_species
    
    
    if(i %in% seq(0, nits-1, by = check_num)){
      pars <- update_proposal_var_beta(pars, beta_accept_post, i, check_num)
      pars <- update_proposal_var_phi(pars, phi_accept_post, i, check_num)
      pars <- update_proposal_var_beta_0(pars, beta_0_accept_post, i, check_num)
    }
    
    
    pb$tick()
    
  }
  
  return(list(beta = beta_post,
              beta_accept = beta_accept_post,
              phi = phi_post,
              phi_accept = phi_accept_post,
              sig_prop_beta = pars$sig_prop_beta,
              sig_prop_phi = pars$sig_prop_phi,
              beta_0 = beta_0_post,
              beta_0_accept = beta_0_accept_post,
              sig_prop_beta_0 = pars$sig_prop_beta_0,
              omega = omega_post,
              sig_prop_omega = pars$sig_prop_omega,
              sigma_species = sigma_species_post))
  
}

run = sampler(fish_dat, mean_covs, mean_covs_log, catch_covs, nits = 50000, check_num = 100)

# saveRDS(run, file = '/Users/joshuanorth/Desktop/full_model_w_doy_temp_int.rds')
# saveRDS(run, file = '/Users/joshuanorth/Desktop/full_model_w_doy_temp_int_cov.rds')
# saveRDS(run, file = '/Users/joshuanorth/Desktop/full_model_w_doy_temp_int_cov_corrected.rds')
# run = readRDS('/Users/joshuanorth/Desktop/full_model_w_doy_temp_int.rds')
# run = readRDS('/Users/joshuanorth/Desktop/full_model_w_doy_temp_int_cov.rds')
run = readRDS('/Users/joshuanorth/Desktop/full_model_w_doy_temp_int_cov_corrected.rds')


# nits = 10000
burnin = 1:25000

apply(run$beta[,,-burnin], c(1,2), mean)
apply(run$beta_accept[,,-burnin], c(1,2), mean)
run$sig_prop_beta

apply(run$phi[,,-burnin], c(1,2), mean)
apply(run$phi_accept[,,-burnin], c(1,2), mean)
run$sig_prop_phi

plot(run$beta_0[-burnin], type = 'l')

pars <- create_pars(fish_dat, mean_covs, mean_covs_log, catch_covs)
b_names = colnames(pars$X[[1]])
phi_names = colnames(pars$Z[[1]])

apply(run$sigma_species[,,-burnin], c(1,2), mean)
cmat = cov2cor(apply(run$sigma_species[,,-burnin], c(1,2), mean))
colnames(cmat) = pars$fish_names
rownames(cmat) = pars$fish_names
cmat

par(mfrow = c(3,4))
for(i in 1:11){
  plot(run$beta[1,i,-burnin], type = 'l', main = b_names[i])
}


par(mfrow = c(3,4))
for(i in 1:11){
  plot(run$phi[1,i,-burnin], type = 'l', main = phi_names[i])
}


par(mfrow = c(2,3))
for(i in 1:6){
  plot(run$omega[i,-burnin], type = 'l', main = pars$fish_names[i])
}

# par(mfrow = c(2,3))
# for(i in 1:6){
#   plot(run$omega[i,], type = 'l', main = pars$fish_names[i])
# }


par(mfrow = c(2,3))
for(i in 1:6){
  plot(run$sigma_species[i,1,], type = 'l', main = pars$fish_names[i])
}



# xtable output -----------------------------------------------------------

pars <- create_pars(fish_dat, mean_covs, mean_covs_log, catch_covs)
b_names = colnames(pars$X[[1]])
phi_names = colnames(pars$Z[[1]])
fnames = pars$fish_names


# relative abundance

b_mean = apply(run$beta[,,-burnin], c(1,2), mean)
b_lower = apply(run$beta[,,-burnin], c(1,2), quantile, probs = 0.025)
b_upper = apply(run$beta[,,-burnin], c(1,2), quantile, probs = 0.975)

b_sig = array(0, dim = dim(b_mean))
b_sig = ifelse((b_upper > 0) & (b_lower < 0), 1, 0)

colnames(b_mean) = b_names
rownames(b_mean) = fnames

# b_col = ifelse(b_sig, paste0("\\cellcolor{red}", round(b_mean, 3)), paste0("\\cellcolor{white}", round(b_mean, 3)))
# colnames(b_col) = b_names
# rownames(b_col) = fnames
# print(xtable(b_col), sanitize.text.function = identity)
# print(xtable(t(b_col)), sanitize.text.function = identity)

b_mean_sig = round(b_mean, 3)
cint = round((b_upper - b_lower)/2, 3)

b_mean_sig = paste0(b_mean_sig, " (", cint, ")")
b_col = ifelse(b_sig, paste0("\\cellcolor{red}", b_mean_sig), paste0("\\cellcolor{white}", b_mean_sig))
colnames(b_col) = b_names
rownames(b_col) = fnames
print(xtable(t(b_col)), sanitize.text.function = identity)


# catchability

phi_mean = apply(run$phi[,,-burnin], c(1,2), mean)
phi_lower = apply(run$phi[,,-burnin], c(1,2), quantile, probs = 0.025)
phi_upper = apply(run$phi[,,-burnin], c(1,2), quantile, probs = 0.975)

p_sig = array(0, dim = dim(phi_mean))
p_sig = ifelse((phi_upper > 0) & (phi_lower < 0), 1, 0)

p_mean_sig = round(phi_mean, 3)
cint = round((phi_upper - phi_lower)/2, 3)

p_mean_sig = paste0(p_mean_sig, " (", cint, ")")
p_col = ifelse(p_sig, paste0("\\cellcolor{red}", p_mean_sig), paste0("\\cellcolor{white}", p_mean_sig))
colnames(p_col) = phi_names
rownames(p_col) = fnames
print(xtable(t(p_col)), sanitize.text.function = identity)



# species dependence

sig_mean = apply(run$sigma_species[,,-burnin], c(1,2), mean)
sig_lower = apply(run$sigma_species[,,-burnin], c(1,2), quantile, probs = 0.025)
sig_upper = apply(run$sigma_species[,,-burnin], c(1,2), quantile, probs = 0.975)

cor_struc = array(NA, dim = dim(run$sigma_species[,,-burnin]))
strt = length(burnin)
for(i in 1:(dim(run$sigma_species)[3] - burnin)){
  cor_struc[,,i] = cov2cor(run$sigma_species[,, strt + i])
}

sig_mean = apply(cor_struc, c(1,2), mean)
sig_lower = apply(cor_struc, c(1,2), quantile, probs = 0.025)
sig_upper = apply(cor_struc, c(1,2), quantile, probs = 0.975)


sig_sig = array(0, dim = dim(sig_mean))
sig_sig = ifelse((sig_upper > 0) & (sig_lower < 0), 1, 0)

sig_mean_sig = round(sig_mean, 3)
cint = round((sig_upper - sig_lower)/2, 3)

sig_mean_sig = paste0(sig_mean_sig, " (", cint, ")")
sig_col = ifelse(sig_sig, paste0("\\cellcolor{red}", sig_mean_sig), paste0("\\cellcolor{white}", sig_mean_sig))
colnames(sig_col) = fnames
rownames(sig_col) = fnames
print(xtable(t(sig_col)), sanitize.text.function = identity)


apply(run$sigma_species[,,-burnin], c(1,2), mean)
cmat = cov2cor(apply(run$sigma_species[,,-burnin], c(1,2), mean))
colnames(cmat) = pars$fish_names
rownames(cmat) = pars$fish_names
cmat


# relative abundance plot -------------------------------------------------

# run = read_rds('/Users/joshuanorth/Desktop/full_model.rds')

apply(run$omega[,-burnin], 1, mean)

apply(run$sigma_species[,,-burnin], c(1,2), mean)
apply(run$sigma_species[,,-burnin], c(1,2), quantile, probs = 0.025)
apply(run$sigma_species[,,-burnin], c(1,2), quantile, probs = 0.975)
cov2cor(apply(run$sigma_species[,,-burnin], c(1,2), mean))

pars <- create_pars(fish_dat, mean_covs, mean_covs_log, catch_covs)
b_names = colnames(pars$X[[1]])
phi_names = colnames(pars$Z[[1]])

b_hat = apply(run$beta[,,-c(burnin)], c(1,2), mean)
colnames(b_hat) = b_names

lam_hat = list()
for(k in 1:pars$K){
  lam_hat[[k]] = exp(pars$X[[k]] %*% b_hat[k,] + omega_hat[k])
}


rel_abun = tibble(lake_index = pars$lake_index, 
                  as_tibble(matrix(unlist(lam_hat), ncol = pars$K)),
                  .name_repair = 'unique')
colnames(rel_abun) = c('DOW', pars$fish_names)

rel_abun = rel_abun %>% 
  rename_all(~str_replace_all(., " ", "_")) %>% 
  pivot_longer(cols = -DOW, names_to = "Fish", values_to = "Abundance") %>% 
  group_by(Fish) %>% 
  distinct(DOW, .keep_all = T) %>% 
  ungroup() %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW')

# plot(rel_abun$bluegill ~ rel_abun$`largemouth bass`)

lats = range(rel_abun$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(rel_abun$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = rel_abun, 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = log(Abundance)), size = 3) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~Fish, 
             labeller = labeller(Fish = c('black_crappie' = 'Black Crappie',
                                          'bluegill' = 'Bluegill',
                                          'northern_pike' = 'Northern Pike',
                                          'yellow_perch' = 'Yellow Perch',
                                          'largemouth_bass' = 'Largemouth Bass',
                                          'walleye' = 'Walleye'))) +
  ggtitle("Relative Log Abundance") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/non_spatial_results/relative_log_abund.png', width = 12, height = 12)

# relative abundance by CPUE  ---------------------------------------------

rel_abun = tibble(lake_index = pars$lake_index, 
                  as_tibble(matrix(unlist(lam_hat), ncol = pars$K)),
                  .name_repair = 'unique')
colnames(rel_abun) = c('DOW', pars$fish_names)

rel_abun = rel_abun %>% 
  pivot_longer(cols = -DOW, names_to = "COMMON_NAME", values_to = "Abundance") %>% 
  mutate(COMMON_NAME = as.factor(COMMON_NAME)) %>% 
  group_by(COMMON_NAME) %>% 
  distinct(DOW, .keep_all = T) %>% 
  ungroup() %>% 
  left_join(fish_dat, by = c('DOW', 'COMMON_NAME')) %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW')



lats = range(rel_abun$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(rel_abun$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))


cpue_rel = rel_abun %>% 
  # filter(year(SURVEYDATE) == 2014) %>% 
  filter(CPUE < 80) %>%
  mutate(yr = year(SURVEYDATE)) %>% 
  group_by(DOW, COMMON_NAME, yr) %>% 
  mutate(CPUE = sum(CPUE, na.rm = T)) %>% 
  distinct(DOW, .keep_all = T) %>% 
  ungroup() %>% 
  select(COMMON_NAME, Abundance, CPUE, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5)

ggplot(cpue_rel, aes(x = CPUE, y = Abundance)) +
  geom_point(size = 3) +
  facet_wrap(~COMMON_NAME, 
             scales = 'free', 
             labeller = labeller(COMMON_NAME = c('black crappie' = 'Black Crappie',
                                                 'bluegill' = 'Bluegill',
                                                 'largemouth bass' = 'Largemouth Bass',
                                                 'northern pike' = 'Northern Pike',
                                                 'walleye' = 'Walleye',
                                                 'yellow perch' = 'Yellow Perch'))) +
  ylab('Relative Abuncance') +
  xlab('Catch Per Unit Effort') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/non_spatial_results/cpue_rel.png', width = 12, height = 10)


cpue_rel_spat = rel_abun %>% 
  # filter(year(SURVEYDATE) == 2014) %>% 
  filter(CPUE < 80) %>%
  mutate(yr = year(SURVEYDATE)) %>% 
  group_by(DOW, COMMON_NAME, yr) %>% 
  mutate(CPUE = sum(CPUE, na.rm = T)) %>% 
  distinct(DOW, .keep_all = T) %>% 
  ungroup() %>% 
  select(COMMON_NAME, Abundance, CPUE, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
  pivot_longer(c(Abundance, CPUE), names_to = 'Est', values_to = 'Abun')



p1 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'black crappie'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Black Crappie - CPUE") +
  guides(color = F)

p2 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'black crappie'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Black Crappie - Abundance") +
  guides(color = F)

p3 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'bluegill'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Bluegill - CPUE") +
  guides(color = F)

p4 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'bluegill'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Bluegill - Abundance") +
  guides(color = F)

p5 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'largemouth bass'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Largemouth Bass - CPUE") +
  guides(color = F)

p6 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'largemouth bass'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Largemouth Bass - Abundance") +
  guides(color = F)

p7 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'northern pike'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Northern Pike - CPUE") +
  guides(color = F)

p8 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'northern pike'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Northern Pike - Abundance") +
  guides(color = F)

p9 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'walleye'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Walleye - CPUE") +
  guides(color = F)

p10 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'walleye'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Walleye - Abundance") +
  guides(color = F)

p11 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'yellow perch'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Yellow Perch - CPUE") +
  guides(color = F)

p12 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'yellow perch'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Yellow Perch - Abundance") +
  guides(color = F)

p1_grid = cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, byrow = F)
p2_grid = cowplot::plot_grid(p7, p8, p9, p10, p11, p12, nrow = 2, byrow = F)
ggsave('results/non_spatial_results/cpue_rel_spat1.png',p1_grid, width = 15, height = 10)
ggsave('results/non_spatial_results/cpue_rel_spat2.png',p2_grid, width = 15, height = 10)


# seasonal plots ----------------------------------------------------------

pars = create_pars(fish_dat, mean_covs, catch_covs)
phi_names = colnames(pars$Z[[1]])


alpha_b = apply(run$phi[,,-burnin], c(1,2), mean)
alpha_b_lower = apply(run$phi[,,-burnin], c(1,2), quantile, probs = c(0.025))
alpha_b_upper = apply(run$phi[,,-burnin], c(1,2), quantile, probs = c(0.975))


fnames = c('crappie', 'bluegill', 'bass', 'pike', 'walleye', 'perch')
# lake = '38025600'
# # 06015200 78002500 34007900 46010900
# lake = '06015200'
# 
# temp_raw = read_rds('data/daily_degree_days_MN_lakes.rds') %>% ungroup()
# 
# temp_raw = temp_raw %>%
#   select(date, temp_0, MNDOW_ID) %>% # C5 temperature
#   rename(SURVEYDATE = date) %>%
#   mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>%
#   select(-MNDOW_ID) %>% 
#   mutate(DOY = yday(SURVEYDATE)) %>% 
#   group_by(SURVEYDATE) %>% 
#   mutate(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0)) %>% 
#   ungroup()
# 
# temp %>% 
#   filter(year(SURVEYDATE) == 2005) %>% 
#   group_by(SURVEYDATE) %>% 
#   summarise(mean(temp_0)) %>% 
#   ungroup()
# filter(DOW == '06015200')


d_plot = function(alpha_b, alpha_b_upper, alpha_b_lower, fish_names, yr, fish_dat){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  # tC = GDD %>% 
  #   filter(year == yr)
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == yr) %>% 
    select(catch_covs, GN, SURVEYDATE) %>% 
    mutate_at(vars(catch_covs), .funs = list(GN = ~.*GN))
  
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:TN) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == yr) %>% 
    pivot_longer(GN:TN, names_to = "Gear", values_to = "Ind") %>% 
    mutate(Gear = factor(Gear)) %>% 
    filter(Ind == 1) %>% 
    group_by(Gear, SURVEYDATE) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    ungroup() %>% 
    select(-EFFORT) %>% 
    spread(key = Gear, value = Ind, drop = F, fill = 0)
  
  TN_df = tibble(DD5 = 0, DOY = 1:365) %>% 
    mutate(DOY_sin = sin(DOY/365 * 2*pi),
           DOY_cos = cos(DOY/365 * 2*pi),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>% 
    select(-DOY) %>% 
    mutate(GN = 0) %>% 
    mutate_at(vars(catch_covs), .funs = list(GN = ~.*GN))
  
  TN_df_lower = tibble(DD5 = min(year_select$DD5), DOY = 1:365) %>% 
    mutate(DOY_sin = sin(DOY/365 * 2*pi),
           DOY_cos = cos(DOY/365 * 2*pi),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>% 
    select(-DOY) %>% 
    mutate(GN = 0) %>% 
    mutate_at(vars(catch_covs), .funs = list(GN = ~.*GN))
  
  TN_df_upper = tibble(DD5 = max(year_select$DD5), DOY = 1:365) %>% 
    mutate(DOY_sin = sin(DOY/365 * 2*pi),
           DOY_cos = cos(DOY/365 * 2*pi),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>% 
    select(-DOY) %>% 
    mutate(GN = 0) %>% 
    mutate_at(vars(catch_covs), .funs = list(GN = ~.*GN))
  
  GN_df = tibble(DD5 = 0, DOY = 1:365) %>% 
    mutate(DOY_sin = sin(DOY/365 * 2*pi),
           DOY_cos = cos(DOY/365 * 2*pi),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>% 
    select(-DOY) %>% 
    mutate(GN = 1) %>% 
    mutate_at(vars(catch_covs), .funs = list(GN = ~.*GN))
  
  GN_df_lower = tibble(DD5 = min(year_select$DD5), DOY = 1:365) %>% 
    mutate(DOY_sin = sin(DOY/365 * 2*pi),
           DOY_cos = cos(DOY/365 * 2*pi),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>% 
    select(-DOY) %>% 
    mutate(GN = 1) %>% 
    mutate_at(vars(catch_covs), .funs = list(GN = ~.*GN))
  
  GN_df_upper = tibble(DD5 = max(year_select$DD5), DOY = 1:365) %>% 
    mutate(DOY_sin = sin(DOY/365 * 2*pi),
           DOY_cos = cos(DOY/365 * 2*pi),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>% 
    select(-DOY) %>% 
    mutate(GN = 1) %>% 
    mutate_at(vars(catch_covs), .funs = list(GN = ~.*GN))
  
  
  
  TN = as.matrix(TN_df) %*% t(alpha_b)
  GN = as.matrix(GN_df) %*% t(alpha_b)
  
  # TN_lower = as.matrix(TN_df_lower) %*% t(alpha_b_lower)
  # GN_lower = as.matrix(GN_df_lower) %*% t(alpha_b_lower)
  # 
  # TN_upper = as.matrix(TN_df_upper) %*% t(alpha_b_upper)
  # GN_upper = as.matrix(GN_df_upper) %*% t(alpha_b_upper)
  
  TN_lower = as.matrix(TN_df_lower) %*% t(alpha_b)
  GN_lower = as.matrix(GN_df_lower) %*% t(alpha_b)
  
  TN_upper = as.matrix(TN_df_upper) %*% t(alpha_b)
  GN_upper = as.matrix(GN_df_upper) %*% t(alpha_b)
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  colnames(TN_lower) = paste0("TN_",fish_names)
  colnames(GN_lower) = paste0("GN_",fish_names)
  
  colnames(TN_upper) = paste0("TN_",fish_names)
  colnames(GN_upper) = paste0("GN_",fish_names)
  
  cnames = c(paste0("TN_",fish_names),
             paste0("GN_",fish_names))
  
  all_days = tibble(SURVEYDATE = seq(ymd(paste0(yr, '-01-01')),ymd(paste0(yr, '-12-31')), by = '1 day'))
  
  date_select = year_select %>% 
    distinct(SURVEYDATE) %>% 
    mutate(sample_day = SURVEYDATE) %>% 
    right_join(all_days, by = c('SURVEYDATE')) %>% 
    left_join(gear_ind, by = c('SURVEYDATE')) %>% 
    mutate_at(vars(GN:TN), ~ if_else(. == 1, SURVEYDATE, NaN)) %>% 
    select(-sample_day) %>% 
    arrange(SURVEYDATE) %>% 
    rename_at(vars(-SURVEYDATE), ~ paste0(., '_survey'))
  
  mean = tibble(as_tibble(TN), as_tibble(GN), date_select)
  lower = tibble(as_tibble(TN_lower), as_tibble(GN_lower), date_select)
  upper = tibble(as_tibble(TN_upper), as_tibble(GN_upper), date_select)
  
  return(list(mean = mean,
              lower = lower,
              upper = upper))
}


plt_dat = d_plot(alpha_b, alpha_b_upper, alpha_b_lower, fish_names, yr =2010, fish_dat)

plt_dat$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness") %>% 
  left_join(plt_dat$upper %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_upper"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat$lower %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_lower"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>%
  left_join(plt_dat$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-c(sample_day)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(size = 1) +
  facet_wrap(~ Gear) +
  geom_ribbon(aes(ymin = Effectiveness_lower, ymax = Effectiveness_upper), alpha = 0.2) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-03-01')), ymd(paste0(yr, '-12-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
  ggtitle(paste0('Gear Type Effectiveness for ', yr)) +
  xlab('Date') +
  ylab('Relative Effectiveness')



# posterior effectiveness -------------------------------------------------

d_plot = function(phi, fish_names, yr, fish_dat, DD5_val){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == yr) %>% 
    select(catch_covs, GN, SURVEYDATE) %>% 
    mutate_at(vars(catch_covs), .funs = list(GN = ~.*GN))
  
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:TN) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == yr) %>% 
    pivot_longer(GN:TN, names_to = "Gear", values_to = "Ind") %>% 
    mutate(Gear = factor(Gear)) %>% 
    filter(Ind == 1) %>% 
    group_by(Gear, SURVEYDATE) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    ungroup() %>% 
    select(-EFFORT) %>% 
    spread(key = Gear, value = Ind, drop = F, fill = 0)
  
  TN_df = tibble(DD5 = DD5_val, DOY = 1:365) %>% 
    mutate(DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>% 
    select(-DOY) %>% 
    mutate(GN = 0) %>% 
    mutate_at(vars(catch_covs), .funs = list(GN = ~.*GN))
  
  GN_df = tibble(DD5 = DD5_val, DOY = 1:365) %>% 
    mutate(DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>% 
    select(-DOY) %>% 
    mutate(GN = 1) %>% 
    mutate_at(vars(catch_covs), .funs = list(GN = ~.*GN))
  
  TN = as.matrix(TN_df)
  GN = as.matrix(GN_df)
  
  n_runs = dim(phi)[3]
  
  post_stores_TN = post_stores_GN = array(NA, dim = c(nrow(TN), length(fnames), n_runs))
  
  for(i in 1:n_runs){
    post_stores_TN[,,i] = TN %*% t(phi[,,i])
    post_stores_GN[,,i] = GN %*% t(phi[,,i])
  }
  
  TN = apply(post_stores_TN, c(1,2), mean)
  TN_lower = apply(post_stores_TN, c(1,2), quantile, probs = 0.025)
  TN_upper = apply(post_stores_TN, c(1,2), quantile, probs = 0.975)
  
  GN = apply(post_stores_GN, c(1,2), mean)
  GN_lower = apply(post_stores_GN, c(1,2), quantile, probs = 0.025)
  GN_upper = apply(post_stores_GN, c(1,2), quantile, probs = 0.975)
  

  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  colnames(TN_lower) = paste0("TN_",fish_names)
  colnames(GN_lower) = paste0("GN_",fish_names)
  
  colnames(TN_upper) = paste0("TN_",fish_names)
  colnames(GN_upper) = paste0("GN_",fish_names)
  
  cnames = c(paste0("TN_",fish_names),
             paste0("GN_",fish_names))
  
  all_days = tibble(SURVEYDATE = seq(ymd(paste0(yr, '-01-01')),ymd(paste0(yr, '-12-31')), by = '1 day')) %>% 
    mutate(DOY = yday(SURVEYDATE)) %>% 
    filter(DOY <= 365) %>% 
    select(-DOY)
  
  date_select = year_select %>% 
    distinct(SURVEYDATE) %>% 
    mutate(sample_day = SURVEYDATE) %>% 
    right_join(all_days, by = c('SURVEYDATE')) %>% 
    left_join(gear_ind, by = c('SURVEYDATE')) %>% 
    mutate_at(vars(GN:TN), ~ if_else(. == 1, SURVEYDATE, NaN)) %>% 
    select(-sample_day) %>% 
    arrange(SURVEYDATE) %>% 
    rename_at(vars(-SURVEYDATE), ~ paste0(., '_survey'))
  
  mean = tibble(as_tibble(TN), as_tibble(GN), date_select)
  lower = tibble(as_tibble(TN_lower), as_tibble(GN_lower), date_select)
  upper = tibble(as_tibble(TN_upper), as_tibble(GN_upper), date_select)
  
  return(list(mean = mean,
              lower = lower,
              upper = upper))
}

yr = 1994
dd5_range = fish_dat %>% 
  filter(year(SURVEYDATE) == yr) %>% 
  summarize(range = range(DD5),
            mean = mean(DD5))
fnames = c('crappie', 'bluegill', 'bass', 'pike', 'walleye', 'perch')
phi = run$phi[,,-burnin]
# range(fish_dat$DD5) # -3.518076  2.788096
plt_dat = d_plot(phi, fish_names = fnames, yr, fish_dat = fish_dat, DD5_val = dd5_range$mean[1])
plt_dat_low = d_plot(phi, fish_names = fnames, yr, fish_dat = fish_dat, DD5_val = dd5_range$range[1])
plt_dat_high = d_plot(phi, fish_names = fnames, yr, fish_dat = fish_dat, DD5_val = dd5_range$range[2])

## mean
p1 = plt_dat$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness") %>% 
  left_join(plt_dat$upper %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_upper"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat$lower %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_lower"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>%
  left_join(plt_dat$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-c(sample_day)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(size = 1) +
  facet_wrap(~ Gear) +
  geom_ribbon(aes(ymin = Effectiveness_lower, ymax = Effectiveness_upper), alpha = 0.2) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-04-01')), ymd(paste0(yr, '-10-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
  ggtitle(paste0('Gear Type Effectiveness for ', yr, ' Mean GDD')) +
  xlab('Date') +
  ylab('Relative Effectiveness') +
  ylim(c(-3, 3))


## lower
p2 = plt_dat_low$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness") %>% 
  left_join(plt_dat_low$upper %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_upper"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat_low$lower %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_lower"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>%
  left_join(plt_dat_low$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-c(sample_day)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(size = 1) +
  facet_wrap(~ Gear) +
  geom_ribbon(aes(ymin = Effectiveness_lower, ymax = Effectiveness_upper), alpha = 0.2) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-04-01')), ymd(paste0(yr, '-10-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
  ggtitle(paste0('Gear Type Effectiveness for ', yr, ' Lower Bound GDD')) +
  xlab('Date') +
  ylab('Relative Effectiveness') +
  ylim(c(-3, 3))


## upper
p3 = plt_dat_high$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness") %>% 
  left_join(plt_dat_high$upper %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_upper"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat_high$lower %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_lower"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>%
  left_join(plt_dat_high$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-c(sample_day)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(size = 1) +
  facet_wrap(~ Gear) +
  geom_ribbon(aes(ymin = Effectiveness_lower, ymax = Effectiveness_upper), alpha = 0.2) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-04-01')), ymd(paste0(yr, '-10-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
  ggtitle(paste0('Gear Type Effectiveness for ', yr, ' Upper Bound GDD')) +
  xlab('Date') +
  ylab('Relative Effectiveness') +
  ylim(c(-3, 3))

p1 = p1 + guides(color = F, fill = F)
p2 = p2 + xlab('') + ylab('') + guides(color = F, fill = F)
p3 = p3 + xlab('') + ylab('')

cowplot::plot_grid(p1, p2, p3, nrow = 1)
ggsave(paste0('results/non_spatial_results/effectiveness_', yr, '.png'), width = 20)



# north south comparison --------------------------------------------------

lake_dow_south = 53002800 # south
lake_dow_north = 16001900 # north
yr = 2015

GDD_center = GDD %>% 
  mutate(DD5  = (DD5 - mean(DD5))/sd(DD5))

lake_south = GDD_center %>% 
  filter(DOW == lake_dow_south) %>% 
  filter(year == yr)

lake_north = GDD_center %>% 
  filter(DOW == lake_dow_north) %>% 
  filter(year == yr)
  
fnames = c('crappie', 'bluegill', 'bass', 'pike', 'walleye', 'perch')
phi = run$phi[,,-burnin]


plt_dat_south = d_plot(phi, fish_names = fnames, yr, fish_dat = fish_dat, DD5_val = lake_south$DD5[1])
plt_dat_north = d_plot(phi, fish_names = fnames, yr, fish_dat = fish_dat, DD5_val = lake_north$DD5[1])

plt_dat_south$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_South") %>% 
  left_join(plt_dat_north$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_North"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  pivot_longer(cols = -c(SURVEYDATE:Fish), names_to = c("Type", "Lake"), names_sep = "_", values_to = "Effectiveness") %>% 
  select(-Type) %>% 
  left_join(plt_dat_south$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-sample_day) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(aes(linetype = Lake), size = 1.2) +
  facet_wrap(~ Gear) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-05-01')), ymd(paste0(yr, '-10-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
  ggtitle(paste0('Gear Type Effectiveness for a North and South lake')) +
  xlab('Date') +
  ylab('Relative Effectiveness')

# ggsave(paste0('results/non_spatial_results/catchability_lake_comp.png'), width = 10, height = 6)


# catchability - temperature ----------------------------------------

pars <- create_pars(fish_dat, mean_covs, mean_covs_log, catch_covs)
phi_names = colnames(pars$Z[[1]])


alpha_b = apply(run$phi[,,-burnin], c(1,2), mean)
alpha_b_lower = apply(run$phi[,,-burnin], c(1,2), quantile, probs = c(0.025))
alpha_b_upper = apply(run$phi[,,-burnin], c(1,2), quantile, probs = c(0.975))

fnames = c('crappie', 'bluegill', 'bass', 'pike', 'walleye', 'perch')

colnames(alpha_b) = phi_names
rownames(alpha_b) = fnames
alpha_b

lake_pos = fish_dat %>% 
  select(DOW, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>% 
  distinct(DOW, .keep_all = T) %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW')

# 16001900
lake_pos %>% 
  filter(LAKE_CENTER_LAT_DD5 > quantile(LAKE_CENTER_LAT_DD5, probs = 0.975)) %>% 
  filter(LAKE_CENTER_LONG_DD5 > quantile(LAKE_CENTER_LONG_DD5, probs = 0.975))

# 53002800
lake_pos %>% 
  filter(LAKE_CENTER_LAT_DD5 < quantile(LAKE_CENTER_LAT_DD5, probs = 0.025)) %>% 
  filter(LAKE_CENTER_LONG_DD5 < quantile(LAKE_CENTER_LONG_DD5, probs = 0.025))

# 18005000
lake_pos %>% 
  filter((LAKE_CENTER_LAT_DD5 < quantile(LAKE_CENTER_LAT_DD5, probs = 0.52)) & (LAKE_CENTER_LAT_DD5 > quantile(LAKE_CENTER_LAT_DD5, probs = 0.48))) %>% 
  filter((LAKE_CENTER_LONG_DD5 < quantile(LAKE_CENTER_LONG_DD5, probs = 0.52)) & (LAKE_CENTER_LONG_DD5 > quantile(LAKE_CENTER_LONG_DD5, probs = 0.48)))




d_plot_lake = function(phis, fish_names, yr, fish_dat, lake_dow){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  tC = temp %>%
    filter(year(SURVEYDATE) == yr) %>%
    filter(DOW == lake_dow) %>% 
    filter(DOY != 366) %>% 
    select(-c(DOW, DD5, DOY))

  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == yr) %>% 
    select(SURVEYDATE, all_of(catch_covs), GN) %>% 
    mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN))
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:TN) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == yr) %>% 
    pivot_longer(GN:TN, names_to = "Gear", values_to = "Ind") %>% 
    mutate(Gear = factor(Gear)) %>% 
    filter(Ind == 1) %>% 
    group_by(Gear, SURVEYDATE) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    ungroup() %>% 
    select(-EFFORT) %>% 
    spread(key = Gear, value = Ind, drop = F, fill = 0)
  
  
  TN_df = year_select %>% 
    filter(GN != 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 0, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>%
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  GN_df = year_select %>% 
    filter(GN == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 1, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>%
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  TN = as.matrix(TN_df)
  GN = as.matrix(GN_df)
  
  n_runs = dim(phis)[3]
  
  post_stores_TN = post_stores_GN = array(NA, dim = c(nrow(tC), length(fnames), n_runs))
  
  for(i in 1:n_runs){
    post_stores_TN[,,i] = TN %*% t(phis[,,i])
    post_stores_GN[,,i] = GN %*% t(phis[,,i])
  }
  
  TN = apply(post_stores_TN, c(1,2), mean)
  TN_lower = apply(post_stores_TN, c(1,2), quantile, probs = 0.025)
  TN_upper = apply(post_stores_TN, c(1,2), quantile, probs = 0.975)
  
  GN = apply(post_stores_GN, c(1,2), mean)
  GN_lower = apply(post_stores_GN, c(1,2), quantile, probs = 0.025)
  GN_upper = apply(post_stores_GN, c(1,2), quantile, probs = 0.975)
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  colnames(TN_lower) = paste0("TN_",fish_names)
  colnames(GN_lower) = paste0("GN_",fish_names)
  
  colnames(TN_upper) = paste0("TN_",fish_names)
  colnames(GN_upper) = paste0("GN_",fish_names)
  
  
  date_select = year_select %>% 
    distinct(SURVEYDATE) %>% 
    mutate(sample_day = SURVEYDATE) %>% 
    right_join(tC, by = c('SURVEYDATE')) %>% 
    left_join(gear_ind, by = c('SURVEYDATE')) %>% 
    select(-c(temp_0)) %>% 
    mutate_at(vars(GN:TN), ~ if_else(. == 1, SURVEYDATE, NaN)) %>% 
    select(-sample_day) %>% 
    arrange(SURVEYDATE) %>% 
    rename_at(vars(-SURVEYDATE), ~ paste0(., '_survey'))
  
  mean = tibble(as_tibble(TN), as_tibble(GN), date_select)
  lower = tibble(as_tibble(TN_lower), as_tibble(GN_lower), date_select)
  upper = tibble(as_tibble(TN_upper), as_tibble(GN_upper), date_select)
  
  return(list(mean = mean,
              lower = lower,
              upper = upper))
}

phis = run$phi[,,-burnin]
yr = 1994

lake_dow_south = 53002800 # south
lake_dow_north = 16001900 # north
lake_dow_center = 18005000 # center
plt_dat_south = d_plot_lake(phis, fnames, yr, fish_dat, lake_dow_south)
plt_dat_north = d_plot_lake(phis, fnames, yr, fish_dat, lake_dow_north)
plt_dat_center = d_plot_lake(phis, fnames, yr, fish_dat, lake_dow_center)


plt_dat_south$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_South") %>% 
  left_join(plt_dat_north$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_North"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat_center$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_Cetner"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  pivot_longer(cols = -c(SURVEYDATE:Fish), names_to = c("Type", "Lake"), names_sep = "_", values_to = "Effectiveness") %>% 
  select(-Type) %>% 
  left_join(plt_dat_south$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-sample_day) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(aes(linetype = Lake), size = 0.9) +
  facet_wrap(~ Gear) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-04-01')), ymd(paste0(yr, '-10-01')))) +
  scale_linetype_manual(values=c("solid", 'longdash', "dotted")) +
  ggtitle(paste0('Gear Type Effectiveness by Region, ', yr)) +
  xlab('') +
  ylab('Relative Effectiveness') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ylim(c(-25, 25))

# ggsave(paste0('results/non_spatial_results/catchability_lake_comp_', yr, '.png'), width = 10, height = 6)


plt_dat = plt_dat_south
plt_dat = plt_dat_north
plt_dat = plt_dat_center

plt_dat$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness") %>% 
  left_join(plt_dat$upper %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_upper"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat$lower %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_lower"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>%
  left_join(plt_dat$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-c(sample_day)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(size = 0.8) +
  facet_wrap(~ Gear) +
  geom_ribbon(aes(ymin = Effectiveness_lower, ymax = Effectiveness_upper), alpha = 0.2) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-04-01')), ymd(paste0(yr, '-10-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle(paste0('Gear Type Effectiveness for ', yr)) +
  xlab('') +
  ylab('Relative Effectiveness') +
  ylim(c(-30, 30))

# ggsave(paste0('results/non_spatial_results/catchability_south_lake_', yr, '.png'), width = 10, height = 6)
# ggsave(paste0('results/non_spatial_results/catchability_north_lake_', yr, '.png'), width = 10, height = 6)
# ggsave(paste0('results/non_spatial_results/catchability_cental_lake_', yr, '.png'), width = 10, height = 6)

