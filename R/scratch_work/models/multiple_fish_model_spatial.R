##----------------------------------------------------
## Name: Joshua North
##
## Date: 11/30/2020
##
## Project: Fish Abundance
##
## Objective: Simple Poisson model for all species, no time, elliptical slice sampler
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
library(Matrix)



# load in data ------------------------------------------------------------

static = read_csv('data/Static_MN_Data_JE.csv')
effort = read_csv('data/MN_fish_catch_effort_all_lakes_years_surveys_JE.csv')
temp = read_rds('data/daily_degree_days_MN_lakes.rds') %>% ungroup()


# join data ---------------------------------------------------------------


static = static %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
effort = effort %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))



# data cleaning -----------------------------------------------------------


fish_dat = effort %>% 
  left_join(static, by = 'DOW') %>%
  select(DOW, COMMON_NAME, TOTAL_CATCH, EFFORT, GEAR, SURVEYDATE, MAX_DEPTH_FEET, mean.gdd, LAKE_AREA_GIS_ACRES, secchi.m, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>%
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME),
         GEAR = as.factor(GEAR),
         SURVEYDATE = mdy(SURVEYDATE),
         DOY = yday(SURVEYDATE),
         DOY = (DOY - mean(seq(1,365)))/sd(seq(1,365))) %>%
  filter(complete.cases(.)) %>% 
  filter(COMMON_NAME == 'yellow perch' | COMMON_NAME == 'northern pike' | COMMON_NAME == 'walleye' | COMMON_NAME == 'largemouth bass') %>% 
  filter(SURVEYDATE >= '1993-01-01') %>% 
  filter(GEAR != 'GSH') %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW),
         GEAR = droplevels(GEAR)) %>% 
  arrange(DOW)

all <- fish_dat %>% 
  group_by(SURVEYDATE, DOW) %>% 
  tidyr::expand(COMMON_NAME, SURVEYDATE, GEAR) %>% 
  ungroup() %>% 
  arrange(DOW)

fish_dat <- fish_dat %>% 
  right_join(all) %>%
  mutate(EFFORT = coalesce(EFFORT, 0L),
         TOTAL_CATCH = coalesce(TOTAL_CATCH, 0L)) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:DOY), list(~coalesce(., 0L))) %>% 
  arrange(SURVEYDATE) %>% 
  group_by(GEAR, SURVEYDATE, DOW) %>% 
  mutate(EFFORT = max(EFFORT)) %>%
  ungroup() %>% 
  group_by(DOW, SURVEYDATE) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:DOY), ~ mean(.[!is.na(.) & . != 0])) %>% 
  ungroup() %>% 
  mutate(row = row_number(),
         Gear_ind = as.integer(1)) %>% 
  pivot_wider(names_from = GEAR, values_from = Gear_ind, values_fill = 0) %>% 
  select(-c(row)) %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW)) %>%
  mutate(DOY = ifelse(is.na(DOY), 0, DOY)) %>% 
  arrange(DOW)




# fish_dat %>% filter(SURVEYDATE == '2009-08-03' & DOW == '01000100')



# fish_dat %>% filter(COMMON_NAME == 'northern pike')
# fish_dat %>% filter(COMMON_NAME == 'yellow perch')



# add temperature ---------------------------------------------------------


# fish_dat = fish_dat %>% filter(SURVEYDATE > '1980-01-01')
temp = temp %>%
  select(date, temp_0, MNDOW_ID) %>% # C5 temperature
  rename(SURVEYDATE = date) %>% 
  mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>% 
  select(-MNDOW_ID)


fish_dat = fish_dat %>% 
  inner_join(temp, by = c('SURVEYDATE', 'DOW')) %>% 
  mutate(DOW = factor(DOW))

# nrow(fish_dat %>% filter(TOTAL_CATCH == 0 & EFFORT != 0))/nrow(fish_dat)

rm(all)
# rm(temp, static, effort, all)


# reduce for now
nlevs = length(levels(fish_dat$DOW))
keepers = sample(levels(fish_dat$DOW), floor(nlevs*0.2))
fish_dat = fish_dat %>% 
  filter(DOW %in% keepers) %>% 
  mutate(DOW = droplevels(DOW))

# simple model ------------------------------------------------------------

create_pars <- function(fish_dat){
  
  K = nlevels(fish_dat$COMMON_NAME)
  levs = levels(fish_dat$COMMON_NAME)
  lake_index = (fish_dat %>% filter(COMMON_NAME == levs[1]))$DOW # lake_index == lake_id[2]
  lake_id = levels(fish_dat$DOW)
  n_lakes = length(levels(fish_dat$DOW))
  
  pars = list()
  
  # data
  pars$Y = list()
  pars$X = list()
  pars$effort = list()
  
  for(k in 1:K){
    pars$Y[[k]] = (fish_dat %>% filter(COMMON_NAME == levs[k]) %>% select(TOTAL_CATCH))$TOTAL_CATCH
  }
  
  for(k in 1:K){
    X = fish_dat %>% 
      filter(COMMON_NAME == levs[k]) %>% 
      mutate(Int = 1) %>% 
      select(Int, MAX_DEPTH_FEET:DOY, temp_0) %>% 
      mutate_at(vars(MAX_DEPTH_FEET:LAKE_AREA_GIS_ACRES), ~ ifelse(. == 0, . + 0.001, .)) %>% 
      mutate_at(vars(MAX_DEPTH_FEET:LAKE_AREA_GIS_ACRES), ~ log(.)) %>% 
      mutate(day_tmp = DOY*temp_0) %>% 
      select(-c(LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING, mean.gdd, MAX_DEPTH_FEET))
    Z = fish_dat %>% 
      filter(COMMON_NAME == levs[k]) %>%
      select(EF, EW, GN, SE)
    # select(EF, EW, GN, GSH, SE)
    pars$X[[k]] = as.matrix(tibble(X, Z) %>% 
                              mutate(DOY_EF = DOY * EF,
                                     DOY_EW = DOY * EW,
                                     DOY_GN = DOY * GN,
                                     DOY_SE = DOY * SE,
                                     temp_EF = temp_0 * EF,
                                     temp_EW = temp_0 * EW,
                                     temp_GN = temp_0 * GN,
                                     temp_SE = temp_0 * SE,
                                     day_tmp_EF = day_tmp * EF,
                                     day_tmp_EW = day_tmp * EW,
                                     day_tmp_GN = day_tmp * GN,
                                     day_tmp_SE = day_tmp * SE))
    # pars$X[[k]] = as.matrix(cbind(X, Z))
  }
  
  for(k in 1:K){
    pars$effort[[k]] = (fish_dat %>% filter(COMMON_NAME == levs[k]) %>% select(EFFORT))$EFFORT
  }
  
  
  # parameters
  pars$n = unlist(lapply(pars$Y, length))
  pars$p = ncol(pars$X[[1]])
  pars$K = K
  
  pars$beta = array(0, dim = c(K, pars$p))
  pars$beta_accept =  array(0, dim = c(K, pars$p))
  pars$beta_prior_var = 100
  
  # hyperpriors
  pars$Sigma_species = diag(K)
  pars$nu_species = K + 2
  pars$Psi_species = 1*diag(K)
  
  # spatial parameters
  spat_dat = fish_dat %>% 
    distinct(DOW, .keep_all = T) %>% 
    select(DOW, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>% 
    arrange(DOW)
  
  d = rdist(cbind(spat_dat$LAKE_CENTER_UTM_EASTING, spat_dat$LAKE_CENTER_UTM_NORTHING))/1000
  phi = max(d)/3
  pars$Sigma_spatial = exp(-(d^2)/phi)
  pars$spatial_var = 1
  pars$Sigma_spatial_inv = solve(pars$Sigma_spatial)
  
  # U = t(chol(kronecker(pars$Sigma_species, pars$Sigma_spatial)))
  # b = rnorm(pars$K * n_lakes)
  
  # pars$eta = matrix(U%*%b, nrow = n_lakes, ncol = pars$K)
  
  pars$eta = rep(0, n_lakes)
  
  
  # species parameters
  pars$W = rep(0, pars$K)
  pars$W_accept = rep(0, pars$K)
  pars$w_prior_var = 100
  
  # Proposal variances
  pars$sig_prop = array(0.5, dim = c(K, pars$p))
  pars$sig_prop_w = array(0.5, dim = c(K))
  
  
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
  effort = pars$effort
  eta = pars$eta
  
  # parameters
  n = pars$n
  p = pars$p
  K = pars$K
  
  # indexing
  lake_index = pars$lake_index
  lake_id = pars$lake_id
  n_lakes = pars$n_lakes
  
  # beta monitor values
  beta_accept = array(0, dim = c(K, p))
  beta_curr = pars$beta
  sig_prop = pars$sig_prop
  beta_prior_var = pars$beta_prior_var
  
  # set up spatial random effect
  ind_array = data.frame(id = lake_id, eta)
  lake_array = data.frame(id = lake_index)
  ETA = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))
  
  # species random effect
  W = pars$W
  
  # for(i in 1:p){
  #   for(k in 1:K){
  #     
  #     b_prop = rnorm(1, beta_curr[k, i], sig_prop[k, i])
  #     beta_prop = beta_curr
  #     beta_prop[k,i] = b_prop
  #     
  #     like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta_curr[k,] + ETA + W[k]), log = T)) + dnorm(beta_curr[k, i], 0, beta_prior_var)
  #     like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta_prop[k,] + ETA + W[k]), log = T)) + dnorm(b_prop, 0, beta_prior_var)
  #     
  #     if((like_prop - like_curr) > log(runif(1))){
  #       beta_curr[k,i] = b_prop
  #       beta_accept[k,i] = 1
  #     }
  #     
  #   }
  #   
  # }
  
  # beta_prop = matrix(rnorm(K*p, c(beta_curr), c(sig_prop)), nrow = K, ncol = p)
  beta_prop = beta_curr
  
  for(i in 1:p){
    for(k in 1:K){
      
      beta_prop = beta_curr
      
      b_prop = beta_prop[k,i] = rnorm(1, beta_curr[k,i], sig_prop[k,i])
      b_curr = beta_curr[k,i]
      
      # like_curr = like_prop = 0
      # for(j in 1:K){
      #   like_curr = like_curr + sum(dpois(Y[[j]], lambda = effort[[j]]*exp(X[[k]] %*% beta_curr[j,] + eta[j]), log = T))
      #   like_prop = like_prop + sum(dpois(Y[[j]], lambda = effort[[j]]*exp(X[[k]] %*% beta_prop[j,] + eta[j]), log = T))
      # }
      # 
      # like_curr = like_curr + dnorm(b_curr, 0, beta_prior_var, log = T)
      # like_prop = like_prop + dnorm(b_prop, 0, beta_prior_var, log = T)
      
      
      
      like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta_curr[k,] + eta[k]), log = T)) + dnorm(b_curr, 0, beta_prior_var, log = T)
      like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta_prop[k,] + eta[k]), log = T)) + dnorm(b_prop, 0, beta_prior_var, log = T)
      
      
      
      if((like_prop - like_curr) > log(runif(1))){
        beta_curr[k,i] = b_prop
        beta_accept[k,i] = 1
      }
      
    }
    
  }
  
  pars$beta = beta_curr
  pars$beta_accept = beta_accept
  
  return(pars)
  
  
}

update_eta <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  W = pars$W
  effort = pars$effort
  beta = pars$beta
  eta = pars$eta
  Sigma_species = pars$Sigma_species
  Sigma_spatial = pars$Sigma_spatial
  spatial_var = pars$spatial_var
  
  C = Matrix(spatial_var * Sigma_spatial)
  
  # parameters
  n = pars$n
  p = pars$p
  K = pars$K
  
  # indexing
  lake_index = pars$lake_index
  lake_id = pars$lake_id
  n_lakes = pars$n_lakes
  
  # choose ellipse v
  U = t(chol(C))
  b = rnorm(n_lakes)
  v = as.matrix(U%*%b)
  
  # set up species random effect
  ind_array = data.frame(id = lake_id, eta)
  lake_array = data.frame(id = lake_index)
  ETA = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))
  
  ll_calc <- function(eta, lake_id, lake_index, beta, K, Y, X, effort, W){
    
    # set up species random effect
    ind_array = data.frame(id = lake_id, eta)
    lake_array = data.frame(id = lake_index)
    ETA = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))
    
    like = 0
    for(k in 1:K){
      like = like + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + ETA + W[k]), log = T))
    }
    
    return(like)
    
  }
  
  Ly = ll_calc(eta, lake_id, lake_index, beta, K, Y, X, effort, W) + log(runif(1))
  
  theta = runif(1, 0, 2*pi)
  theta_max = theta
  theta_min = theta_max - 2*pi
  
  eta_prop = eta*cos(theta) + v*sin(theta)
  curr_like = ll_calc(eta_prop, lake_id, lake_index, beta, K, Y, X, effort, W)
  
  while(curr_like < Ly){
    
    if(theta < 0){
      theta_min = theta
    }else{
      theta_max = theta
    }
    
    theta = runif(1, theta_min, theta_max)
    
    eta_prop = eta*cos(theta) + v*sin(theta)
    curr_like = ll_calc(eta_prop, lake_id, lake_index, beta, K, Y, X, effort, W)
    
  }
  
  pars$eta = eta_prop

  return(pars)
  
  
}

update_w <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  beta = pars$beta
  effort = pars$effort
  eta = pars$eta
  
  Sigma_species = pars$Sigma_species
  
  # parameters
  n = pars$n
  p = pars$p
  K = pars$K
  
  # indexing
  lake_index = pars$lake_index
  lake_id = pars$lake_id
  n_lakes = pars$n_lakes
  
  # W monitor values
  w_accept = array(0, dim = c(K))
  w_curr = pars$W
  sig_prop_w = pars$sig_prop_w
  w_prior_var = pars$w_prior_var
  
  # set up spatial random effect
  ind_array = data.frame(id = lake_id, eta)
  lake_array = data.frame(id = lake_index)
  ETA = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))
  
  # propose values
  # w_prop = rmvnorm(1, w_curr , Sigma_species)
  # w_prop = rmvnorm(1, w_curr , diag(sig_prop_w))
  # 
  # for(k in 1:K){
  #   
  #   w_tmp = w_curr
  #   w_tmp[k] = w_prop[k]
  #   
  #   like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + ETA + w_curr[k]), log = T)) + 
  #     dmvnorm(w_curr, rep(0, K), Sigma_species, log = T)
  #   like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + ETA + w_prop[k]), log = T)) +
  #     dmvnorm(w_tmp, rep(0, K), Sigma_species, log = T)
  #   
  #   if((like_prop - like_curr) > log(runif(1))){
  #     w_curr[k] = w_prop[k]
  #     w_accept[k] = 1
  #   }
  #   
  # }
  
  # beta_prop = beta_curr
  # 
  # for(i in 1:p){
  #   for(k in 1:K){
  #     
  #     b_prop = beta_prop[k,i] = rnorm(1, beta_curr[k,i], sig_prop[k,i])
  #     b_curr = beta_curr[k,i]
  
  w_prop = w_curr
  
  for(k in 1:K){
    
    w_prop[k] = rnorm(1, w_curr[k], sig_prop_w[k])
    
    like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + ETA + w_curr[k]), log = T)) + dnorm(w_curr[k], 0, w_prior_var, log = T)
    like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + ETA + w_prop[k]), log = T)) + dnorm(w_prop[k], 0, w_prior_var, log = T)
      
    if((like_prop - like_curr) > log(runif(1))){
      w_curr[k] = w_prop[k]
      w_accept[k] = 1
    }
    
  }
  
  pars$W = w_curr
  pars$W_accept = w_accept
  
  return(pars)
  
  
}

update_sigma_species <- function(pars){
  
  # data
  X = pars$X
  beta = pars$beta
  Sigma_spatial_inv = pars$Sigma_spatial_inv
  spatial_var = pars$spatial_var
  eta = pars$eta
  W = pars$W
  K = pars$K

  nu_species = pars$nu_species
  Psi_species = pars$Psi_species
  
  nu_hat = nu_species + K
  psi_hat = Psi_species + W %*% t(W)
  
  pars$Sigma_species = MCMCpack::riwish(nu_hat, psi_hat)
  
  return(pars)
  
}

update_spatial_var <- function(pars){
  
  # data
  X = pars$X
  beta = pars$beta
  Sigma_spatial_inv = pars$Sigma_spatial_inv
  Sigma_species = pars$Sigma_species
  eta = pars$eta
  W = pars$W
  K = pars$K
  
  # parameters
  n_lakes = pars$n_lakes
  
  a_hat = 2 + n_lakes
  b_hat = 2 + c(t(eta) %*% Sigma_spatial_inv %*% eta)
  
  pars$spatial_var = 1/rgamma(1, a_hat, b_hat)
  
  return(pars)
  
}

update_proposal_var <- function(pars, beta_accept_post, i, check_num = 100){
  
  sig_prop = pars$sig_prop
  
  bp = beta_accept_post[,,(i-check_num):i]
  accept_rate = apply(bp, c(1,2), mean)
  
  sig_prop = ifelse(accept_rate < 0.25, sig_prop*0.7, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.7, sig_prop)
  
  pars$sig_prop = sig_prop
  
  return(pars)
  
  
}

update_proposal_var_w <- function(pars, W_accept_post, i, check_num = 100){
  
  sig_prop = pars$sig_prop_w
  
  bp = W_accept_post[,(i-check_num):i]
  accept_rate = apply(bp, c(1), mean)
  
  sig_prop = ifelse(accept_rate < 0.25, sig_prop*0.7, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.7, sig_prop)
  
  pars$sig_prop_w = sig_prop
  
  return(pars)
  
  
}

sampler <- function(fish_dat, nits, check_num = 100){
  
  pars = create_pars(fish_dat)
  
  p = pars$p
  K = pars$K
  n_lake = pars$n_lakes
  
  beta_post = array(NA, dim = c(K, p, nits))
  beta_accept_post = array(NA, dim = c(K, p, nits))
  W_post = array(NA, dim = c(K, nits))
  W_accept_post = array(NA, dim = c(K, nits))
  eta_post = array(NA, dim = c(n_lake, nits))
  sigma_species_post = array(NA, dim = c(dim(pars$Sigma_species), nits))
  spatial_var_post = array(NA, dim = c(dim(pars$spatial_var), nits))
  
  pb <- progress_bar$new(
    format = "  Running [:bar] :percent eta: :eta",
    total = nits, clear = FALSE, width = 60)
  
  for(i in seq(1, nits)){
    
    pars <- update_beta(pars)
    pars <- update_eta(pars)
    pars <- update_w(pars)
    pars <- update_sigma_species(pars)
    pars <- update_spatial_var(pars)
    
    beta_post[,,i] = pars$beta
    beta_accept_post[,,i] = pars$beta_accept
    W_post[,i] = pars$W
    W_accept_post[,i] = pars$W_accept
    eta_post[,i] = pars$eta
    sigma_species_post[,,i] = pars$Sigma_species
    spatial_var_post[i] = pars$spatial_var
    
    
    if(i %in% seq(0, nits-1, by = check_num)){
      pars <- update_proposal_var(pars, beta_accept_post, i, check_num)
      pars <- update_proposal_var_w(pars, W_accept_post, i, check_num)
    }
    
    
    pb$tick()
    
  }
  
  return(list(beta = beta_post,
              beta_accept = beta_accept_post,
              W = W_post,
              W_accept = W_accept_post,
              sig_prop = pars$sig_prop,
              eta = eta_post,
              sigma_species = sigma_species_post,
              spatial_var = spatial_var_post))
  
}

nits = 50000
burnin = 1:25000

run = sampler(fish_dat, nits, check_num = 500)


pars = create_pars(fish_dat)
nms = colnames(pars$X[[1]])

par(mfrow = c(2,2))
plot(run$beta[1, 1,-c(burnin)], type = 'l', main = nms[1])
plot(run$beta[2, 1,-c(burnin)], type = 'l', main = nms[1])
plot(run$beta[3, 1,-c(burnin)], type = 'l', main = nms[1])
plot(run$beta[4, 1,-c(burnin)], type = 'l', main = nms[1])

plot(run$beta[1, 2,-c(burnin)], type = 'l', main = nms[2])
plot(run$beta[2, 2,-c(burnin)], type = 'l', main = nms[2])

plot(run$beta[1, 3,-c(burnin)], type = 'l', main = nms[3])
plot(run$beta[2, 3,-c(burnin)], type = 'l', main = nms[3])

plot(run$beta[1, 4,-c(burnin)], type = 'l', main = nms[4])
plot(run$beta[2, 4,-c(burnin)], type = 'l', main = nms[4])

plot(run$beta[1, 5,-c(burnin)], type = 'l', main = nms[5])
plot(run$beta[2, 5,-c(burnin)], type = 'l', main = nms[5])

plot(run$beta[1, 6,-c(burnin)], type = 'l', main = nms[6])
plot(run$beta[2, 6,-c(burnin)], type = 'l', main = nms[6])

plot(run$beta[1, 7,-c(burnin)], type = 'l', main = nms[7])
plot(run$beta[2, 7,-c(burnin)], type = 'l', main = nms[7])

plot(run$beta[1, 8,-c(burnin)], type = 'l', main = nms[8])
plot(run$beta[2, 8,-c(burnin)], type = 'l', main = nms[8])

plot(run$beta[1, 9,-c(burnin)], type = 'l', main = nms[9])
plot(run$beta[2, 9,-c(burnin)], type = 'l', main = nms[9])

plot(run$beta[1, 10,-c(burnin)], type = 'l', main = nms[10])
plot(run$beta[2, 10,-c(burnin)], type = 'l', main = nms[10])

plot(run$beta[1, 11,-c(burnin)], type = 'l', main = nms[11])
plot(run$beta[2, 11,-c(burnin)], type = 'l', main = nms[11])

plot(run$beta[1, 12,-c(burnin)], type = 'l', main = nms[12])
plot(run$beta[2, 12,-c(burnin)], type = 'l', main = nms[12])

plot(run$beta[1, 13,-c(burnin)], type = 'l', main = nms[12])
plot(run$beta[2, 13,-c(burnin)], type = 'l', main = nms[12])


apply(run$beta_accept, c(1,2), mean)
round(run$sig_prop, 3)

cbind(nms, t(apply(run$beta, c(1,2), mean)))

apply(run$sigma_species, c(1,2), mean)
cov2cor(apply(run$sigma_species, c(1,2), mean))
plot(run$sigma_species[1, 1,-c(burnin)], type = 'l', main = nms[12])
plot(run$sigma_species[1, 2,-c(burnin)], type = 'l', main = nms[12])
plot(run$sigma_species[2, 1,-c(burnin)], type = 'l', main = nms[12])
plot(run$sigma_species[2, 2,-c(burnin)], type = 'l', main = nms[12])
sig = apply(run$sigma_species[,,-c(burnin)], c(1,2), mean)
cov = cov2cor(apply(run$sigma_species[,,-c(burnin)], c(1,2), mean))
colnames(sig) = c('Pike', 'Perch')
write.table(sig, "results/spatial_results/dependence_mat.txt")
write.table(cov, "results/spatial_results/correlation_mat.txt")

plot(run$eta[1,-c(burnin)], type = 'l')
plot(run$eta[2,-c(burnin)], type = 'l')
plot(run$eta[3,-c(burnin)], type = 'l')
plot(run$eta[4,-c(burnin)], type = 'l')


plot(run$eta[1,-c(burnin)]*run$beta[1, 1,-c(burnin)], type = 'l')
plot(run$eta[1,-c(burnin)]*run$beta[1, 2,-c(burnin)], type = 'l')
plot(run$eta[1,-c(burnin)]*run$beta[1, 2,-c(burnin)], type = 'l')
plot(run$eta[1,-c(burnin)]*run$beta[1, 2,-c(burnin)], type = 'l')

plot(run$eta[2,1,-c(burnin)]*run$beta[1, 1,-c(burnin)], type = 'l')
plot(run$eta[2,2,-c(burnin)]*run$beta[1, 2,-c(burnin)], type = 'l')

mean(run$spatial_var)

# relative abundance ------------------------------------------------------

pars = create_pars(fish_dat)
nms = colnames(pars$X[[1]])

b_hat = apply(run$beta[,,-c(burnin)], c(1,2), mean)
colnames(b_hat) = nms

int_inds = 2:3
int_b = b_hat[,int_inds]

eta_hat = apply(run$eta[,-c(burnin)], 1, mean)

w_hat = apply(run$W[,-c(burnin)], 1, mean)

ind_array = data.frame(id = pars$lake_id, eta_hat)
lake_array = data.frame(id = pars$lake_index)
ETA = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))

lam_hat = list()
for(k in 1:pars$K){
  lam_hat[[k]] = exp(pars$X[[k]][,int_inds] %*% int_b[k,] + ETA + w_hat[k])
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

lats = range(rel_abun$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(rel_abun$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_jitter(data = rel_abun, 
              aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance),
              width = 0.2, height = 0.2) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Relative Abundance") +
  facet_wrap(~Fish, 
             labeller = labeller(Fish = c('northern_pike' = 'Northern Pike',
                                          'yellow_perch' = 'Yellow Perch',
                                          'largemouth_bass' = 'Largemouth Bass',
                                          'walleye' = 'Walleye')))
# ggsave('results/spatial_results/relative_abun.png')


ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_jitter(data = rel_abun %>% filter(Fish == 'northern_pike'), 
              aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance),
              width = 0.2, height = 0.2) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Relative Abundance - Pike")
# ggsave('results/spatial_results/relative_abun_pike.png')

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_jitter(data = rel_abun %>% filter(Fish == 'yellow_perch'), 
              aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance),
              width = 0.2, height = 0.2) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Relative Abundance - Perch")
ggsave('results/spatial_results/relative_abun_perch.png')


ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_jitter(data = rel_abun %>% filter(Fish == 'walleye'), 
              aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance),
              width = 0.2, height = 0.2) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Relative Abundance - Walleye")
# ggsave('results/spatial_results/relative_abun_walleye.png')


ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_jitter(data = rel_abun %>% filter(Fish == 'largemouth_bass'), 
              aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance),
              width = 0.2, height = 0.2) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Relative Abundance - Bass")
# ggsave('results/spatial_results/relative_abun_bass.png')

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_jitter(data = rel_abun, 
              aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Spatial),
              width = 0.1, height = 0.1) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Spatial Random Effect") +
  facet_wrap(~Fish, labeller = labeller(Fish = c('northern_pike' = 'Northern Pike',
                                                 'yellow_perch' = 'Yellow Perch')))
# ggsave('results/spatial_results/spatial_re.png')


ggplot(rel_abun, aes(y = Spatial, x = LAKE_CENTER_LAT_DD5, col = Abun)) +
  geom_point()

ggplot(rel_abun, aes(y = Spatial, x = LAKE_CENTER_LONG_DD5, col = Abun)) +
  geom_point()


ggplot(fish_dat, aes(x = log(TOTAL_CATCH))) +
  geom_density() +
  facet_wrap(~COMMON_NAME)


# effectiveness of gear ---------------------------------------------------

pars = create_pars(fish_dat)
nms = colnames(pars$X[[1]])

b_hat = apply(run$beta[,,-c(burnin)], c(1,2), mean)
colnames(b_hat) = nms

# alpha_inds = c(1, 4:11)
alpha_inds = c(1, 4:10)
alpha_b = b_hat[,alpha_inds]

b_hat_lower = apply(run$beta[,,-c(burnin)], c(1,2), quantile, probs = c(0.025))
alpha_b_lower = b_hat_lower[,alpha_inds]

b_hat_upper = apply(run$beta[,,-c(burnin)], c(1,2), quantile, probs = c(0.975))
alpha_b_upper = b_hat_upper[,alpha_inds]


fish_dat_sd = fish_dat %>%
  select(SURVEYDATE) %>%
  mutate(DOY = yday(SURVEYDATE)) %>% 
  summarise_at(vars(DOY), list(mean = ~mean(.), 
                               sd = ~sd(.)))


d_plot = function(alpha_b, time, temp, fish_names){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  TN = cbind(1, time, temp, time*temp, 0, 0, 0, 0) %*% t(alpha_b)
  EF = cbind(1, time, temp, time*temp, 1, 0, 0, 0) %*% t(alpha_b)
  EW = cbind(1, time, temp, time*temp, 0, 1, 0, 0) %*% t(alpha_b)
  GN = cbind(1, time, temp, time*temp, 0, 0, 1, 0) %*% t(alpha_b)
  SE = cbind(1, time, temp, time*temp, 0, 0, 0, 1) %*% t(alpha_b)
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(EF) = paste0("EF_",fish_names)
  colnames(EW) = paste0("EW_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  colnames(SE) = paste0("SE_",fish_names)
  
  cnames = c(paste0("TN_",fish_names), 
             paste0("EF_",fish_names), 
             paste0("EW_",fish_names), 
             paste0("GN_",fish_names), 
             paste0("SE_",fish_names))
  
  
  return(tibble(as_tibble(TN), as_tibble(EF), as_tibble(EW), as_tibble(GN), as_tibble(SE)))
}
tC = temp %>% 
  filter(DOW == '71016700') %>% 
  filter(year(SURVEYDATE) == '2008') %>% 
  select(temp_0)


year = '2015'

tC = temp %>% 
  filter(year(SURVEYDATE) == year) %>% 
  mutate(SURVEYDATE = factor(yday(SURVEYDATE))) %>% 
  group_by(SURVEYDATE) %>% 
  summarize(temp_0 = mean(temp_0, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(temp_0)

tC = tC$temp_0
time_len = 1:length(tC)

fnames = c('bass', 'pike', 'walleye', 'perch')

d_plot(alpha_b, (time_len - fish_dat_sd$mean)/fish_dat_sd$sd, tC, fnames) %>% 
  mutate(Time = time_len) %>% 
  pivot_longer(cols = -Time, names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness") %>% 
  left_join(d_plot(alpha_b_upper, (time_len - fish_dat_sd$mean)/fish_dat_sd$sd, tC, fnames) %>% 
              mutate(Time = time_len) %>% 
              pivot_longer(cols = -Time, names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_upper"), by = c('Time', 'Gear', 'Fish')) %>% 
  left_join(d_plot(alpha_b_lower, (time_len - fish_dat_sd$mean)/fish_dat_sd$sd, tC, fnames) %>% 
              mutate(Time = time_len) %>% 
              pivot_longer(cols = -Time, names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_lower"), by = c('Time', 'Gear', 'Fish')) %>% 
  ggplot(., aes(x = Time, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(size = 1) +
  xlim(c(150, 350)) +
  facet_wrap(~ Gear) +
  geom_ribbon(aes(ymin = Effectiveness_lower, ymax = Effectiveness_upper), alpha = 0.2) +
  ggtitle(paste0('Relative Effectiveness, ', year))
ggsave(paste0('results/spatial_results/catchability_', year, '.png'))
  
  










