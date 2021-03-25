##----------------------------------------------------
## Name: Joshua North
##
## Date: 03/25/2020
##
## Project: Fish Abundance
##
## Objective: Simulate data to determine model fit
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

static = static %>% 
  select(DOW, LAKE_CENTER_LAT_DD5:LAKE_CENTER_UTM_NORTHING) %>% 
  mutate(DOW = as.factor(DOW))


# simulate data -----------------------------------------------------------




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

fish_dat = fish_dat %>% 
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/365 * 4*pi),
         DOY_cos_semi = cos(DOY/365 * 4*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0)

fish_dat = fish_dat %>% 
  mutate_at(vars(temp_0, DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN))


colnames(fish_dat)
mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:19, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9, 13:19)]
# catch_covs = colnames(fish_dat)[c(24, 27, 28)] # no temp doy interaction
catch_covs = colnames(fish_dat)[c(24, 27:35)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]
fish_names = levels(fish_dat$COMMON_NAME)


# simulate data -----------------------------------------------------------

set.seed(25032021)

mean_cov_coef = matrix(round(rnorm(length(mean_covs)*length(fish_names), 0.05, 0.15), 2),
                       ncol = length(mean_covs), nrow = length(fish_names))
colnames(mean_cov_coef) = mean_covs
rownames(mean_cov_coef) = fish_names

catch_cov_coef = matrix(round(rnorm(length(catch_covs)*length(fish_names), 0.2, 0.5), 2),
                       ncol = length(catch_covs), nrow = length(fish_names))
colnames(catch_cov_coef) = catch_covs
rownames(catch_cov_coef) = fish_names

mean_cov_coef
catch_cov_coef


cal_gamma <- function(df, f_name){
  
  f_name = f_name$COMMON_NAME
  
  dat = df %>% 
    ungroup() %>% 
    select(all_of(mean_covs)) %>% 
    as.matrix()
  
  covs = mean_cov_coef[rownames(mean_cov_coef) == f_name,]
  
  return(df %>% 
           ungroup() %>% 
           mutate(gamma = c(exp(dat %*% covs))))
}

cal_alpha <- function(df, f_name){
  
  f_name = f_name$COMMON_NAME
  
  dat = df %>% 
    ungroup() %>% 
    select(all_of(catch_covs)) %>% 
    as.matrix()
  
  covs = catch_cov_coef[rownames(catch_cov_coef) == f_name,]
  
  return(df %>% 
           ungroup() %>% 
           mutate(alpha = c(exp(dat %*% covs))))
}

Gamma = fish_dat %>% 
  select(all_of(mean_covs), COMMON_NAME, DOW, SURVEYDATE, TN, GN, EFFORT) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ ifelse(. == 0, . + 0.001, .)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  group_by(COMMON_NAME, .drop = F) %>% 
  group_modify(~ cal_gamma(df = .x, f_name = .y)) %>% 
  ungroup() %>% 
  arrange(DOW, SURVEYDATE)


Alpha = fish_dat %>% 
  select(all_of(catch_covs), COMMON_NAME, DOW, SURVEYDATE, TN, GN, EFFORT) %>% 
  group_by(COMMON_NAME, .drop = F) %>% 
  group_modify(~ cal_alpha(df = .x, f_name = .y)) %>% 
  ungroup() %>% 
  arrange(DOW, SURVEYDATE)


Lambda = Gamma %>% 
  select(COMMON_NAME, DOW:gamma) %>% 
  mutate(DOW = as.factor(DOW)) %>% 
  dplyr::left_join(Alpha %>% 
              select(COMMON_NAME, DOW:alpha) %>% 
              mutate(DOW = as.factor(DOW))) %>% 
  mutate(lambda = EFFORT * gamma * alpha) %>% 
  rowwise() %>% 
  mutate(total_catch = rpois(1, lambda)) %>% 
  ungroup() %>% 
  left_join(fish_dat) %>% 
  select(-c(TOTAL_CATCH, CPUE)) %>%
  rename(TOTAL_CATCH = total_catch) %>% 
  arrange(DOW, SURVEYDATE, TN)


colnames(Lambda)
mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:19, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9, 13:19)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]

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

run = sampler(Lambda, mean_covs, mean_covs_log, catch_covs, nits = 5000, check_num = 50)


# chains ------------------------------------------------------------------

mean_cov_coef
catch_cov_coef

burnin = 1:2500

round(apply(run$beta[,,-burnin], c(1,2), mean), 2)
apply(run$beta_accept[,,-burnin], c(1,2), mean)
run$sig_prop_beta

round(apply(run$phi[,,-burnin], c(1,2), mean), 2)
apply(run$phi_accept[,,-burnin], c(1,2), mean)
run$sig_prop_phi

plot(run$beta_0[-burnin], type = 'l')

pars <- create_pars(Lambda, mean_covs, mean_covs_log, catch_covs)
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


par(mfrow = c(2,3))
for(i in 1:6){
  plot(run$sigma_species[i,1,], type = 'l', main = pars$fish_names[i])
}

