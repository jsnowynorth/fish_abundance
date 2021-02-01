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
  mutate(Gear_ind = as.integer(1)) %>% 
  pivot_wider(names_from = GEAR, values_from = Gear_ind, values_fill = 0) %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW)) %>%
  mutate(DOY = ifelse(is.na(DOY), 0, DOY)) %>% 
  arrange(DOW)

# add temperature 

temp = temp %>%
  select(date, temp_0, MNDOW_ID) %>% # C5 temperature
  rename(SURVEYDATE = date) %>% 
  mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>% 
  select(-MNDOW_ID)


fish_dat = fish_dat %>% 
  inner_join(temp, by = c('SURVEYDATE', 'DOW')) %>% 
  mutate(DOW = factor(DOW))


rm(all)


# # reduce for now
# nlevs = length(levels(fish_dat$DOW))
# keepers = sample(levels(fish_dat$DOW), floor(nlevs*0.5))
# fish_dat = fish_dat %>% 
#   filter(DOW %in% keepers) %>% 
#   mutate(DOW = droplevels(DOW))
# 


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
  pars$Psi_species = 10*diag(K)
  
  pars$eta = matrix(rep(0, K), ncol = K)
  pars$eta_accept =  matrix(rep(0, K), ncol = K)
  pars$eta_prior_var = 100
  
  
  # dt_tmp = tibble(x = spat_dat$LAKE_CENTER_UTM_EASTING/1000, y = spat_dat$LAKE_CENTER_UTM_NORTHING/1000, z = C[,50])
  # ggplot(dt_tmp, aes(x = x, y = y, color = z)) +
  #   geom_point() +
  #   scale_color_gradient(low = 'yellow', high = 'red')
  
  
  # Proposal variances
  pars$sig_prop = array(2, dim = c(K, pars$p))
  pars$sig_prop_eta = array(2, dim = c(K))
  
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
  
  # beta monitor values
  beta_accept = array(0, dim = c(K, p))
  beta_curr = pars$beta
  sig_prop = pars$sig_prop
  beta_prior_var = pars$beta_prior_var

  # for(i in 1:p){
  #   for(k in 1:K){
  #     
  #     b_prop = rnorm(1, beta_curr[k, i], sig_prop[k, i])
  #     beta_prop = beta_curr
  #     beta_prop[k,i] = b_prop
  #     
  #     like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta_curr[k,] + eta[k]), log = T)) + dnorm(beta_curr[k, i], 0, beta_prior_var)
  #     like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta_prop[k,] + eta[k]), log = T)) + dnorm(b_prop, 0, beta_prior_var)
  #     
  #     if((like_prop - like_curr) > log(runif(1))){
  #       beta_curr[k,i] = b_prop
  #       beta_accept[k,i] = 1
  #     }
  #     
  #   }
  #   
  # }
  
  
  # beta_prop = beta_curr
  
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
  effort = pars$effort
  beta = pars$beta
  eta = pars$eta
  sig_prop_eta = pars$sig_prop_eta
  Sigma_species = pars$Sigma_species
  eta_prior_var = pars$eta_prior_var
  
  # parameters
  n = pars$n
  p = pars$p
  K = pars$K
  
  eta_accept = rep(0, K)
  
  
  # propose eta
  eta_curr = eta
  eta_prop = rmvnorm(1, mean = eta_curr, sigma = diag(sig_prop_eta))
  
  for(k in 1:K){
    
    # like_curr = like_prop = 0
    # for(k in 1:K){
    #   like_prop = like_prop + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + eta_prop[k]), log = T))
    #   like_curr = like_curr + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + eta_curr[k]), log = T))
    # }
    # like_prop = like_prop + dnorm(eta_prop[j], 0, eta_prior_var)
    # like_curr = like_curr + dnorm(eta_curr[j], 0, eta_prior_var)
    # 
    # 
    
    like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + eta_prop[k]), log = T)) + dnorm(eta_prop[k], 0, eta_prior_var)
    like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + eta_curr[k]), log = T)) + dnorm(eta_curr[k], 0, eta_prior_var)
    
    
    
    
    if((like_prop - like_curr) > log(runif(1))){
      eta_curr[k] = eta_prop[k]
      eta_accept[k] = 1
    }
    
  }
  
  pars$eta = eta_curr
  pars$eta_accept = eta_accept
  
  return(pars)
  
  
}

update_sigma_species <- function(pars){
  
  # data
  Sigma_species = pars$Sigma_species
  eta = pars$eta
  
  nu_species = pars$nu_species
  Psi_species = pars$Psi_species
  
  # parameters
  K = pars$K
  
  nu_hat = nu_species + K
  psi_hat = Psi_species + t(eta) %*% eta
  
  pars$Sigma_species = MCMCpack::riwish(nu_hat, psi_hat)
  
  return(pars)
  
}

update_proposal_var <- function(pars, beta_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop
  
  bp = beta_accept_post[,,(i-check_num):i]
  accept_rate = apply(bp, c(1,2), mean)
  
  sig_prop = ifelse(accept_rate < 0.2, sig_prop*0.9, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.9, sig_prop)
  
  pars$sig_prop = sig_prop
  
  return(pars)
  
  
}

update_proposal_var_eta <- function(pars, eta_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop_eta
  
  bp = eta_accept_post[,(i-check_num):i]
  accept_rate = apply(bp, c(1), mean)
  
  sig_prop = ifelse(accept_rate < 0.2, sig_prop*0.9, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.9, sig_prop)
  
  pars$sig_prop_eta = sig_prop
  
  return(pars)
  
  
}

sampler <- function(fish_dat, nits, check_num = 200){
  
  pars = create_pars(fish_dat)
  
  p = pars$p
  K = pars$K
  
  beta_post = array(NA, dim = c(K, p, nits))
  beta_accept_post = array(NA, dim = c(K, p, nits))
  eta_post = array(NA, dim = c(length(pars$eta), nits))
  eta_accept_post = array(NA, dim = c(K, nits))
  sigma_species_post = array(NA, dim = c(dim(pars$Sigma_species), nits))
  
  pb <- progress_bar$new(
    format = "  Running [:bar] :percent eta: :eta",
    total = nits, clear = FALSE, width = 60)
  
  for(i in seq(1, nits)){
    
    pars <- update_beta(pars)
    pars <- update_eta(pars)
    pars <- update_sigma_species(pars)
    
    beta_post[,,i] = pars$beta
    beta_accept_post[,,i] = pars$beta_accept
    eta_post[,i] = pars$eta
    eta_accept_post[,i] = pars$eta_accept
    sigma_species_post[,,i] = pars$Sigma_species
    
    
    if(i %in% seq(0, nits-1, by = check_num)){
      pars <- update_proposal_var(pars, beta_accept_post, i, check_num)
      pars <- update_proposal_var_eta(pars, eta_accept_post, i, check_num)
    }
    
    
    pb$tick()
    
  }
  
  return(list(beta = beta_post,
              beta_accept = beta_accept_post,
              sig_prop = pars$sig_prop,
              eta = eta_post,
              eta_accept = eta_accept_post,
              sig_prop_eta = pars$sig_prop_eta,
              sigma_species = sigma_species_post))
  
}

nits = 50000
burnin = 1:25000

run = sampler(fish_dat, nits, check_num = 100)

pars = create_pars(fish_dat)
nms = colnames(pars$X[[1]])

par(mfrow = c(1,2))
plot(run$beta[1, 1,-c(burnin)], type = 'l', main = nms[1])
plot(run$beta[2, 1,-c(burnin)], type = 'l', main = nms[1])
plot(run$beta[3, 1,-c(burnin)], type = 'l', main = nms[1])
plot(run$beta[4, 1,-c(burnin)], type = 'l', main = nms[1])

plot(run$beta[1, 2,-c(burnin)], type = 'l', main = nms[2])
plot(run$beta[2, 2,-c(burnin)], type = 'l', main = nms[2])
plot(run$beta[3, 2,-c(burnin)], type = 'l', main = nms[2])
plot(run$beta[4, 2,-c(burnin)], type = 'l', main = nms[2])

plot(run$beta[1, 3,-c(burnin)], type = 'l', main = nms[3])
plot(run$beta[2, 3,-c(burnin)], type = 'l', main = nms[3])
plot(run$beta[3, 3,-c(burnin)], type = 'l', main = nms[3])
plot(run$beta[4, 3,-c(burnin)], type = 'l', main = nms[3])

plot(run$beta[1, 4,-c(burnin)], type = 'l', main = nms[4])
plot(run$beta[2, 4,-c(burnin)], type = 'l', main = nms[4])
plot(run$beta[3, 4,-c(burnin)], type = 'l', main = nms[4])
plot(run$beta[4, 4,-c(burnin)], type = 'l', main = nms[4])

plot(run$beta[1, 5,-c(burnin)], type = 'l', main = nms[5])
plot(run$beta[2, 5,-c(burnin)], type = 'l', main = nms[5])
plot(run$beta[3, 5,-c(burnin)], type = 'l', main = nms[5])
plot(run$beta[4, 5,-c(burnin)], type = 'l', main = nms[5])

plot(run$beta[1, 6,-c(burnin)], type = 'l', main = nms[6])
plot(run$beta[2, 6,-c(burnin)], type = 'l', main = nms[6])
plot(run$beta[3, 6,-c(burnin)], type = 'l', main = nms[6])
plot(run$beta[4, 6,-c(burnin)], type = 'l', main = nms[6])

plot(run$beta[1, 7,-c(burnin)], type = 'l', main = nms[7])
plot(run$beta[2, 7,-c(burnin)], type = 'l', main = nms[7])
plot(run$beta[3, 7,-c(burnin)], type = 'l', main = nms[7])
plot(run$beta[4, 7,-c(burnin)], type = 'l', main = nms[7])

plot(run$beta[1, 8,-c(burnin)], type = 'l', main = nms[8])
plot(run$beta[2, 8,-c(burnin)], type = 'l', main = nms[8])
plot(run$beta[3, 8,-c(burnin)], type = 'l', main = nms[8])
plot(run$beta[4, 8,-c(burnin)], type = 'l', main = nms[8])

plot(run$beta[1, 9,-c(burnin)], type = 'l', main = nms[9])
plot(run$beta[2, 9,-c(burnin)], type = 'l', main = nms[9])
plot(run$beta[3, 9,-c(burnin)], type = 'l', main = nms[9])
plot(run$beta[4, 9,-c(burnin)], type = 'l', main = nms[9])

plot(run$beta[1, 10,-c(burnin)], type = 'l', main = nms[10])
plot(run$beta[2, 10,-c(burnin)], type = 'l', main = nms[10])
plot(run$beta[3, 10,-c(burnin)], type = 'l', main = nms[10])
plot(run$beta[4, 10,-c(burnin)], type = 'l', main = nms[10])

plot(run$beta[1, 11,-c(burnin)], type = 'l', main = nms[11])
plot(run$beta[2, 11,-c(burnin)], type = 'l', main = nms[11])

plot(run$beta[1, 12,-c(burnin)], type = 'l', main = nms[12])
plot(run$beta[2, 12,-c(burnin)], type = 'l', main = nms[12])

plot(run$beta[1, 13,-c(burnin)], type = 'l', main = nms[12])
plot(run$beta[2, 13,-c(burnin)], type = 'l', main = nms[12])

plot(run$beta[1, 14,-c(burnin)], type = 'l', main = nms[14])
plot(run$beta[2, 14,-c(burnin)], type = 'l', main = nms[14])
plot(run$beta[3, 14,-c(burnin)], type = 'l', main = nms[14])
plot(run$beta[4, 14,-c(burnin)], type = 'l', main = nms[14])

plot(run$beta[1, 15,-c(burnin)], type = 'l', main = nms[15])
plot(run$beta[2, 15,-c(burnin)], type = 'l', main = nms[15])
plot(run$beta[3, 15,-c(burnin)], type = 'l', main = nms[15])
plot(run$beta[4, 15,-c(burnin)], type = 'l', main = nms[15])

plot(run$beta[1, 16,-c(burnin)], type = 'l', main = nms[16])
plot(run$beta[2, 16,-c(burnin)], type = 'l', main = nms[16])
plot(run$beta[3, 16,-c(burnin)], type = 'l', main = nms[16])
plot(run$beta[4, 16,-c(burnin)], type = 'l', main = nms[16])

plot(run$beta[1, 17,-c(burnin)], type = 'l', main = nms[17])
plot(run$beta[2, 17,-c(burnin)], type = 'l', main = nms[17])
plot(run$beta[3, 17,-c(burnin)], type = 'l', main = nms[17])
plot(run$beta[4, 17,-c(burnin)], type = 'l', main = nms[17])

plot(run$beta[1, 18,-c(burnin)], type = 'l', main = nms[18])
plot(run$beta[2, 18,-c(burnin)], type = 'l', main = nms[18])
plot(run$beta[3, 18,-c(burnin)], type = 'l', main = nms[18])
plot(run$beta[4, 18,-c(burnin)], type = 'l', main = nms[18])

plot(run$beta[1, 19,-c(burnin)], type = 'l', main = nms[19])
plot(run$beta[2, 19,-c(burnin)], type = 'l', main = nms[19])
plot(run$beta[3, 19,-c(burnin)], type = 'l', main = nms[19])
plot(run$beta[4, 19,-c(burnin)], type = 'l', main = nms[19])

plot(run$beta[1, 20,-c(burnin)], type = 'l', main = nms[20])
plot(run$beta[2, 20,-c(burnin)], type = 'l', main = nms[20])
plot(run$beta[3, 20,-c(burnin)], type = 'l', main = nms[20])
plot(run$beta[4, 20,-c(burnin)], type = 'l', main = nms[20])


apply(run$beta_accept, c(1,2), mean)
round(run$sig_prop, 3)

apply(run$eta_accept, 1, mean)
round(run$sig_prop_eta, 3)

cbind(nms, t(apply(run$beta, c(1,2), mean)))

apply(run$eta, 1, mean)

apply(run$sigma_species[,,-c(burnin)], c(1,2), mean)
cov2cor(apply(run$sigma_species[,,-c(burnin)], c(1,2), mean))
plot(run$sigma_species[1, 1,-c(burnin)], type = 'l', main = nms[12])
plot(run$sigma_species[1, 2,-c(burnin)], type = 'l', main = nms[12])
plot(run$sigma_species[2, 1,-c(burnin)], type = 'l', main = nms[12])
plot(run$sigma_species[2, 2,-c(burnin)], type = 'l', main = nms[12])
sig = apply(run$sigma_species, c(1,2), mean)
colnames(sig) = c('Pike', 'Perch')
# write.table(sig, "results/non_spatial_results//dependence_mat.txt")

plot(run$eta[1,-c(burnin)], type = 'l')
plot(run$eta[2,-c(burnin)], type = 'l')
plot(run$eta[3,-c(burnin)], type = 'l')
plot(run$eta[4,-c(burnin)], type = 'l')



# save values
keepers = c(12,3)

write.table(rbind(c('par', pars$fish_names), cbind(nms[keepers], t(apply(run$beta, c(1,2), mean))[keepers,])), "results/non_spatial_results/mean_betas.txt")
write.table(rbind(c('par', pars$fish_names), cbind(nms[keepers], t(apply(run$beta, c(1,2), quantile, probs = 0.025))[keepers,])), "results/non_spatial_results/lower_betas.txt")
write.table(rbind(c('par', pars$fish_names), cbind(nms[keepers], t(apply(run$beta, c(1,2), quantile, probs = 0.975))[keepers,])), "results/non_spatial_results/upper_betas.txt")





# relative abundance ------------------------------------------------------

fish_dat %>%
  select(COMMON_NAME, TOTAL_CATCH) %>%
  group_by(COMMON_NAME) %>%
  summarize(mean(TOTAL_CATCH), sd(TOTAL_CATCH), median(TOTAL_CATCH)) #range(TOTAL_CATCH), 



pars = create_pars(fish_dat)
nms = colnames(pars$X[[1]])

b_hat = apply(run$beta[,,-c(burnin)], c(1,2), mean)
colnames(b_hat) = nms

eta_hat = apply(run$eta[,-c(burnin)], c(1), mean)

int_inds = 2:3
int_b = b_hat[,int_inds]

lam_hat = list()
for(k in 1:pars$K){
  lam_hat[[k]] = exp(pars$X[[k]][,int_inds] %*% int_b[k,] + eta_hat[k])
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
ggsave('results/non_spatial_results/relative_abun.png')

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
# ggsave('results/non_spatial_results/relative_abun_pike.png')

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
# ggsave('results/non_spatial_results/relative_abun_perch.png')

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
# ggsave('results/non_spatial_results/relative_abun_bass.png')

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
# ggsave('results/non_spatial_results/relative_abun_walleye.png')



# effectiveness of gear ---------------------------------------------------

pars = create_pars(fish_dat)
nms = colnames(pars$X[[1]])

b_hat = apply(run$beta[,,-c(burnin)], c(1,2), mean)
colnames(b_hat) = nms

# alpha_inds = c(1, 4:11)
alpha_inds = c(1, 4:22)
alpha_b = b_hat[,alpha_inds]

b_hat_lower = apply(run$beta[,,-c(burnin)], c(1,2), quantile, probs = c(0.025))
alpha_b_lower = b_hat_lower[,alpha_inds]

b_hat_upper = apply(run$beta[,,-c(burnin)], c(1,2), quantile, probs = c(0.975))
alpha_b_upper = b_hat_upper[,alpha_inds]


fnames = c('bass', 'pike', 'walleye', 'perch')

d_plot = function(alpha_b, fish_names, year, fish_dat){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  tC = temp %>% 
    filter(year(SURVEYDATE) == year) %>% 
    group_by(SURVEYDATE) %>% 
    summarize(temp_0 = mean(temp_0, na.rm = TRUE)) %>% 
    ungroup() %>% 
    select(temp_0, SURVEYDATE) %>% 
    mutate(DOY = (seq(1, length(temp_0)) - mean(1:365))/sd(1:365))
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == year) %>% 
    mutate(Int = 1) %>% 
    select(SURVEYDATE, Int, EF, EW, GN, SE, DOY, temp_0)
  
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:SE) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == year) %>% 
    pivot_longer(GN:SE, names_to = "Gear", values_to = "Ind") %>% 
    mutate(Gear = factor(Gear)) %>% 
    filter(Ind == 1) %>% 
    group_by(Gear, SURVEYDATE) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    ungroup() %>% 
    select(-EFFORT) %>% 
    spread(key = Gear, value = Ind, drop = F, fill = 0) %>% # pivot_wider(names_from = Gear, values_from = Ind, values_fill = 0) %>% 
    relocate(GN, .after = SURVEYDATE) %>% 
    relocate(SE, .after = last_col())
  
  TN_df = year_select %>% 
    rowwise() %>% 
    filter(sum(c_across(Int:SE)) == 1) %>% 
    ungroup() %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0', 'DOY')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(Int = 1, EF = 0,  EW = 0, GN = 0, SE = 0, day_tmp = DOY*temp_0) %>%
    select(-SURVEYDATE) %>% 
    relocate(DOY:day_tmp, .after = Int) %>% 
    mutate(day_tmp = DOY*temp_0,
           DOY_EF = DOY * EF,
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
           day_tmp_SE = day_tmp * SE)
  
  EF_df = year_select %>% 
    filter(EF == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0', 'DOY')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(Int = 1, EF = 1,  EW = 0, GN = 0, SE = 0, day_tmp = DOY*temp_0) %>%
    select(-SURVEYDATE) %>% 
    relocate(DOY:day_tmp, .after = Int) %>% 
    mutate(day_tmp = DOY*temp_0,
           DOY_EF = DOY * EF,
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
           day_tmp_SE = day_tmp * SE)
  
  EW_df = year_select %>% 
    filter(EF == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0', 'DOY')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(Int = 1, EF = 0,  EW = 1, GN = 0, SE = 0, day_tmp = DOY*temp_0) %>%
    select(-SURVEYDATE) %>% 
    relocate(DOY:day_tmp, .after = Int) %>% 
    mutate(day_tmp = DOY*temp_0,
           DOY_EF = DOY * EF,
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
           day_tmp_SE = day_tmp * SE)
  
  GN_df = year_select %>% 
    filter(EF == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0', 'DOY')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(Int = 1, EF = 0, EW = 0, GN = 1, SE = 0, day_tmp = DOY*temp_0) %>%
    select(-SURVEYDATE) %>% 
    relocate(DOY:day_tmp, .after = Int) %>% 
    mutate(day_tmp = DOY*temp_0,
           DOY_EF = DOY * EF,
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
           day_tmp_SE = day_tmp * SE)
  
  SE_df = year_select %>% 
    filter(SE == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0', 'DOY')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(Int = 1, EF = 0,  EW = 0, GN = 0, SE = 1, day_tmp = DOY*temp_0) %>%
    select(-SURVEYDATE) %>% 
    relocate(DOY:day_tmp, .after = Int) %>% 
    mutate(day_tmp = DOY*temp_0,
           DOY_EF = DOY * EF,
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
           day_tmp_SE = day_tmp * SE)
  
  
  
  TN = as.matrix(TN_df) %*% t(alpha_b)
  EF = as.matrix(EF_df) %*% t(alpha_b)
  EW = as.matrix(EW_df) %*% t(alpha_b)
  GN = as.matrix(GN_df) %*% t(alpha_b)
  SE = as.matrix(SE_df) %*% t(alpha_b)
  
  
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
  
  # date_select = year_select %>% 
  #   distinct(SURVEYDATE) %>% 
  #   mutate(sample_day = SURVEYDATE) %>% 
  #   right_join(tC, by = c('SURVEYDATE')) %>% 
  #   select(-c(temp_0, DOY)) %>% 
  #   arrange(SURVEYDATE)
  
  date_select = year_select %>% 
    distinct(SURVEYDATE) %>% 
    mutate(sample_day = SURVEYDATE) %>% 
    right_join(tC, by = c('SURVEYDATE')) %>% 
    left_join(gear_ind, by = c('SURVEYDATE')) %>% 
    select(-c(temp_0, DOY)) %>% 
    mutate_at(vars(GN:SE), ~ if_else(. == 1, SURVEYDATE, NaN)) %>% 
    select(-sample_day) %>% 
    arrange(SURVEYDATE) %>% 
    rename_at(vars(-SURVEYDATE), ~ paste0(., '_survey'))
  
  return(tibble(as_tibble(TN), as_tibble(EF), as_tibble(EW), as_tibble(GN), as_tibble(SE), date_select))
}

year = 2010
mean = d_plot(alpha_b, fnames, year, fish_dat)
upper = d_plot(alpha_b_upper, fnames, year, fish_dat)
lower = d_plot(alpha_b_lower, fnames, year, fish_dat)


mean %>% select(-c(GN_survey:SE_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness") %>% 
  left_join(upper %>% select(-c(GN_survey:SE_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_upper"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(lower %>% select(-c(GN_survey:SE_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_lower"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>%
  left_join(mean %>%
              select(SURVEYDATE:SE_survey) %>% 
              pivot_longer(GN_survey:SE_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-c(sample_day)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(size = 1) +
  facet_wrap(~ Gear) +
  geom_ribbon(aes(ymin = Effectiveness_lower, ymax = Effectiveness_upper), alpha = 0.2) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(year, '-03-01')), ymd(paste0(year, '-12-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
  ggtitle(paste0('Gear Type Effectiveness for ', year)) +
  xlab('Date') +
  ylab('Relative Effectiveness')

ggsave(paste0('results/non_spatial_results/catchability_', year, '.png'), width = 10, height = 6)





# mean %>% 
#   pivot_longer(cols = -c(SURVEYDATE, sample_day), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness") %>% 
#   left_join(upper %>% 
#               pivot_longer(cols = -c(SURVEYDATE, sample_day), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_upper"), by = c('SURVEYDATE', 'Gear', 'Fish', 'sample_day')) %>% 
#   left_join(lower %>% 
#               pivot_longer(cols = -c(SURVEYDATE, sample_day), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_lower"), by = c('SURVEYDATE', 'Gear', 'Fish', 'sample_day')) %>% 
#   ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
#   geom_line(size = 1) +
#   facet_wrap(~ Gear) +
#   geom_ribbon(aes(ymin = Effectiveness_lower, ymax = Effectiveness_upper), alpha = 0.2) +
#   geom_rug(aes(x = sample_day), sides="b") +
#   scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(year, '-03-01')), ymd(paste0(year, '-12-01')))) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
#   ggtitle(paste0('Gear Type Effectiveness for ', year)) +
#   xlab('Date') +
#   ylab('Relative Effectiveness')
# 
# ggsave(paste0('results/non_spatial_results/catchability_', year, '.png'))


# heatmap of effectiveness ------------------------------------------------

pars = create_pars(fish_dat)
nms = colnames(pars$X[[1]])

b_hat = apply(run$beta[,,-c(burnin)], c(1,2), mean)
colnames(b_hat) = nms

alpha_inds = c(1, 4:22)
alpha_b = b_hat[,alpha_inds]

b_hat_lower = apply(run$beta[,,-c(burnin)], c(1,2), quantile, probs = c(0.025))
alpha_b_lower = b_hat_lower[,alpha_inds]

b_hat_upper = apply(run$beta[,,-c(burnin)], c(1,2), quantile, probs = c(0.975))
alpha_b_upper = b_hat_upper[,alpha_inds]

fnames = c('bass', 'pike', 'walleye', 'perch')

heat_plot = function(alpha_b, fish_names, year, fish_dat){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  tC = temp %>% 
    filter(year(SURVEYDATE) == year) %>% 
    group_by(SURVEYDATE) %>% 
    summarise_at(vars(temp_0), list(lower = ~range(., na.rm = T)[1],
                                    upper = ~range(., na.rm = T)[2])) %>% 
    ungroup() %>% 
    mutate(DOY = (seq(1, n()) - mean(1:365))/sd(1:365))
  
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == year) %>% 
    mutate(Int = 1) %>% 
    select(SURVEYDATE, Int, EF, EW, GN, SE, DOY, temp_0)
  
  TN_df = year_select %>% 
    rowwise() %>% 
    filter(sum(c_across(Int:SE)) == 1) %>% 
    ungroup() %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0', 'DOY')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(Int = 1, EF = 0,  EW = 0, GN = 0, SE = 0, day_tmp = DOY*temp_0) %>%
    select(-SURVEYDATE) %>% 
    relocate(DOY:day_tmp, .after = Int) %>% 
    mutate(day_tmp = DOY*temp_0,
           DOY_EF = DOY * EF,
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
           day_tmp_SE = day_tmp * SE)
  
  EF_df = year_select %>% 
    filter(EF == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0', 'DOY')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(Int = 1, EF = 1,  EW = 0, GN = 0, SE = 0, day_tmp = DOY*temp_0) %>%
    select(-SURVEYDATE) %>% 
    relocate(DOY:day_tmp, .after = Int) %>% 
    mutate(day_tmp = DOY*temp_0,
           DOY_EF = DOY * EF,
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
           day_tmp_SE = day_tmp * SE)
  
  EW_df = year_select %>% 
    filter(EF == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0', 'DOY')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(Int = 1, EF = 0,  EW = 1, GN = 0, SE = 0, day_tmp = DOY*temp_0) %>%
    select(-SURVEYDATE) %>% 
    relocate(DOY:day_tmp, .after = Int) %>% 
    mutate(day_tmp = DOY*temp_0,
           DOY_EF = DOY * EF,
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
           day_tmp_SE = day_tmp * SE)
  
  GN_df = year_select %>% 
    filter(EF == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0', 'DOY')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(Int = 1, EF = 0, EW = 0, GN = 1, SE = 0, day_tmp = DOY*temp_0) %>%
    select(-SURVEYDATE) %>% 
    relocate(DOY:day_tmp, .after = Int) %>% 
    mutate(day_tmp = DOY*temp_0,
           DOY_EF = DOY * EF,
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
           day_tmp_SE = day_tmp * SE)
  
  SE_df = year_select %>% 
    filter(SE == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0', 'DOY')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(Int = 1, EF = 0,  EW = 0, GN = 0, SE = 1, day_tmp = DOY*temp_0) %>%
    select(-SURVEYDATE) %>% 
    relocate(DOY:day_tmp, .after = Int) %>% 
    mutate(day_tmp = DOY*temp_0,
           DOY_EF = DOY * EF,
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
           day_tmp_SE = day_tmp * SE)
  
  
  
  TN = as.matrix(TN_df) %*% t(alpha_b)
  EF = as.matrix(EF_df) %*% t(alpha_b)
  EW = as.matrix(EW_df) %*% t(alpha_b)
  GN = as.matrix(GN_df) %*% t(alpha_b)
  SE = as.matrix(SE_df) %*% t(alpha_b)
  
  
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
  
  date_select = year_select %>% 
    distinct(SURVEYDATE) %>% 
    mutate(sample_day = SURVEYDATE) %>% 
    right_join(tC, by = c('SURVEYDATE')) %>% 
    select(-c(temp_0, DOY)) %>% 
    arrange(SURVEYDATE)
  
  return(tibble(as_tibble(TN), as_tibble(EF), as_tibble(EW), as_tibble(GN), as_tibble(SE), date_select))
}











# species by year by gear type --------------------------------------------

species_gear_count = effort %>% 
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
  inner_join(temp, by = c('SURVEYDATE', 'DOW')) %>% 
  arrange(DOW) %>% 
  mutate(year = year(SURVEYDATE)) %>% 
  select(year, COMMON_NAME, GEAR) %>% 
  group_by(year, COMMON_NAME, GEAR) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = GEAR, values_from = count, values_fill = NA)

fish_dat %>% 
  mutate(year = year(SURVEYDATE)) %>% 
  select(year, COMMON_NAME, GN:SE, EFFORT) %>% 
  filter(EFFORT != 0) %>% 
  select(-EFFORT) %>% 
  pivot_longer(GN:SE, names_to = 'GEAR', values_to = 'Gear_ind') %>% 
  filter(Gear_ind == 1) %>% 
  select(-Gear_ind) %>% 
  group_by(year, COMMON_NAME, GEAR) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = GEAR, values_from = count, values_fill = NA)

tmp = fish_dat %>% 
  mutate(year = year(SURVEYDATE)) %>% 
  select(year, COMMON_NAME, GN:SE, TOTAL_CATCH) %>% 
  filter(TOTAL_CATCH != 0) %>% 
  select(-TOTAL_CATCH) %>% 
  pivot_longer(GN:SE, names_to = 'GEAR', values_to = 'Gear_ind') %>% 
  filter(Gear_ind == 1) %>% 
  group_by(year, COMMON_NAME, GEAR) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = GEAR, values_from = count, values_fill = NA)


effort %>% 
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
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW),
         GEAR = droplevels(GEAR)) %>% 
  arrange(DOW) %>% 
  mutate(year = year(SURVEYDATE)) %>% 
  select(year, COMMON_NAME, GEAR) %>% 
  group_by(year, COMMON_NAME, GEAR) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = GEAR, values_from = count, values_fill = NA)

# write_csv(species_gear_count, 'results/non_spatial_results/species_gear_count.csv')