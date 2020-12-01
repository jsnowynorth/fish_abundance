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
         DOY = (DOY - mean(DOY))/sd(DOY),
         DOY_squared = DOY^2) %>%
  filter(complete.cases(.)) %>% 
  filter(COMMON_NAME == 'yellow perch' | COMMON_NAME == 'northern pike') %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW))


all <- fish_dat %>% group_by(SURVEYDATE, DOW) %>% expand(COMMON_NAME, SURVEYDATE, GEAR) %>% ungroup()
fish_dat <- fish_dat %>% 
  right_join(all) %>%
  mutate(EFFORT = coalesce(EFFORT, 0L),
         TOTAL_CATCH = coalesce(TOTAL_CATCH, 0L)) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:DOY_squared), list(~coalesce(., 0L))) %>% 
  arrange(SURVEYDATE) %>% 
  mutate(row = row_number(),
         Gear_ind = as.integer(1)) %>% 
  pivot_wider(names_from = GEAR, values_from = Gear_ind, values_fill = 0) %>% 
  select(-c(row)) %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW)) %>% 
  arrange(DOW)



# fish_dat %>% filter(COMMON_NAME == 'northern pike')
# fish_dat %>% filter(COMMON_NAME == 'yellow perch')



# add temperature ---------------------------------------------------------


fish_dat = fish_dat %>% filter(SURVEYDATE > '1980-01-01')
temp = temp %>%
  select(date, TempC5, MNDOW_ID) %>% # C5 temperature
  rename(SURVEYDATE = date) %>% 
  mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>% 
  select(-MNDOW_ID)


fish_dat = fish_dat %>% 
  inner_join(temp, by = c('SURVEYDATE', 'DOW')) %>% 
  mutate(DOW = factor(DOW))

rm(temp, static, effort, all)


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
      select(Int, MAX_DEPTH_FEET:DOY_squared, TempC5) %>% 
      mutate_at(vars(MAX_DEPTH_FEET:secchi.m), ~ ifelse(. == 0, . + 0.001, .)) %>% 
      mutate_at(vars(MAX_DEPTH_FEET:secchi.m), ~ log(.)) %>% 
      select(-c(LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING))
    Z = fish_dat %>% 
      filter(COMMON_NAME == levs[k]) %>%
      select(SE:GSH)
    pars$X[[k]] = as.matrix(cbind(X, Z))
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
  
  # hyperpriors
  pars$Sigma_species = diag(K)
  pars$nu_species = K + 2
  pars$Psi_species = 10*diag(K)
  
  # spatial parameters
  
  spat_dat = fish_dat %>% 
    distinct(DOW, .keep_all = T) %>% 
    select(DOW, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>% 
    arrange(DOW)
  
  d = rdist(cbind(spat_dat$LAKE_CENTER_UTM_EASTING, spat_dat$LAKE_CENTER_UTM_NORTHING))/1000
  phi = max(d)/3
  pars$Sigma_spatial = exp(-(d^2)/phi)
  pars$spatial_var = 1
  
  U = t(chol(kronecker(pars$Sigma_species, pars$Sigma_spatial)))
  b = rnorm(pars$K * n_lakes)
  
  pars$eta = matrix(U%*%b, nrow = n_lakes, ncol = pars$K)
  
  
  # dt_tmp = tibble(x = spat_dat$LAKE_CENTER_UTM_EASTING/1000, y = spat_dat$LAKE_CENTER_UTM_NORTHING/1000, z = C[,50])
  # ggplot(dt_tmp, aes(x = x, y = y, color = z)) +
  #   geom_point() +
  #   scale_color_gradient(low = 'yellow', high = 'red')
  
  
  # Proposal variances
  pars$sig_prop = array(0.5, dim = c(K, pars$p))
  
  
  # indexing
  pars$lake_index = lake_index
  pars$lake_id = lake_id
  pars$n_lakes = n_lakes
  
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
  
  # set up species random effect
  ind_array = data.frame(id = lake_id, eta)
  lake_array = data.frame(id = lake_index)
  ETA = lake_array %>% right_join(ind_array, by = 'id') %>% select(-id)
  
  
  for(i in 1:p){
    for(k in 1:K){
      
      b_prop = rnorm(1, beta_curr[k, i], sig_prop[k, i])
      beta_prop = beta_curr
      beta_prop[k,i] = b_prop
      
      like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta_curr[k,] + ETA[,k]), log = T))
      like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta_prop[k,] + ETA[,k]), log = T))
      
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
  Sigma_species = pars$Sigma_species
  Sigma_spatial = pars$Sigma_spatial
  spatial_var = pars$spatial_var
  
  C = kronecker(pars$Sigma_species, spatial_var * pars$Sigma_spatial)
  
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
  b = rnorm(pars$K * n_lakes)
  v = matrix(U%*%b, nrow = n_lakes, ncol = K)
  
  ll_calc <- function(eta, lake_id, lake_index, beta, K, Y, X, effort){
    
    # set up species random effect
    ind_array = data.frame(id = lake_id, eta)
    lake_array = data.frame(id = lake_index)
    ETA = lake_array %>% right_join(ind_array, by = 'id') %>% select(-id)
    
    like = 0
    for(k in 1:K){
      like = like + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta[k,] + ETA[,k]), log = T))
    }
    
    return(like)
    
  }
  
  Ly = ll_calc(eta, lake_id, lake_index, beta, K, Y, X, effort) + log(runif(1))
  
  theta = runif(1, 0, 2*pi)
  theta_max = theta
  theta_min = theta_max - 2*pi
  
  eta_prop = eta*cos(theta) + v*sin(theta)
  curr_like = ll_calc(eta_prop, lake_id, lake_index, beta, K, Y, X, effort)
  
  while(curr_like < Ly){
    
    if(theta < 0){
      theta_min = theta
    }else{
      theta_max = theta
    }
    
    theta = runif(1, theta_min, theta_max)
    # c(theta, theta_min, theta_max)
    # theta_max = theta
    # theta_min = theta_max - 2*pi
    
    eta_prop = eta*cos(theta) + v*sin(theta)
    curr_like = ll_calc(eta_prop, lake_id, lake_index, beta, K, Y, X, effort)
    
  }
  
  pars$eta = eta_prop

  return(pars)
  
  
}

update_sigma_species <- function(pars){
  
  # data
  Sigma_species = pars$Sigma_species
  Sigma_spatial = pars$Sigma_spatial
  spatial_var = pars$spatial_var
  
  nu_species = pars$nu_species
  Psi_species = pars$Psi_species
  
  # parameters
  n_lakes = pars$n_lakes
  
  nu_hat = nu_species + n_lakes
  psi_hat = Psi_species + t(eta) %*% solve(spatial_var * Sigma_spatial) %*% eta
  
  pars$Sigma_species = MCMCpack::riwish(nu_hat, psi_hat)
  
  return(pars)
  
}

update_proposal_var <- function(pars, beta_accept_post, i){
  
  sig_prop = pars$sig_prop
  
  bp = beta_accept_post[,,(i-9):i]
  accept_rate = apply(bp, c(1,2), mean)
  
  sig_prop = ifelse(accept_rate < 0.25, sig_prop*0.7, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.7, sig_prop)
  
  pars$sig_prop = sig_prop
  
  return(pars)
  
  
}

sampler <- function(fish_dat, nits){
  
  pars = create_pars(fish_dat)
  
  p = pars$p
  K = pars$K
  
  beta_post = array(NA, dim = c(K, p, nits))
  beta_accept_post = array(NA, dim = c(K, p, nits))
  eta_post = array(NA, dim = c(dim(pars$eta), nits))
  sigma_species_post = array(NA, dim = c(dim(pars$Sigma_species), nits))
  spatial_var_post = array(NA, dim = c(dim(pars$spatial_var), nits))
  
  pb <- progress_bar$new(
    format = "  Running [:bar] :percent eta: :eta",
    total = nits, clear = FALSE, width = 60)
  
  for(i in seq(1, nits)){
    
    pars <- update_beta(pars)
    pars <- update_eta(pars)
    pars <- update_sigma_species(pars)
    
    beta_post[,,i] = pars$beta
    beta_accept_post[,,i] = pars$beta_accept
    eta_post[,,i] = pars$eta
    sigma_species_post[,,i] = pars$Sigma_species
    spatial_var_post[i] = pars$spatial_var
    
    
    if(i %in% seq(0, nits-1, by = 10)){
      pars <- update_proposal_var(pars, beta_accept_post, i)
    }
    
    
    pb$tick()
    
  }
  
  return(list(beta = beta_post,
              beta_accept = beta_accept_post,
              sig_prop = pars$sig_prop,
              eta = eta_post,
              sigma_species = sigma_species_post,
              spatial_var = spatial_var_post))
  
}

nits = 1000
burnin = 1:500

run = sampler(fish_dat, nits)

pars = create_pars(fish_dat)
nms = colnames(pars$X[[1]])

par(mfrow = c(1,2))
plot(run$beta[1, 1,-c(burnin)], type = 'l', main = nms[1])
plot(run$beta[2, 1,-c(burnin)], type = 'l', main = nms[1])

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

plot(run$eta[1,1,], type = 'l')
plot(run$eta[1,2,], type = 'l')

plot(run$eta[2,1,], type = 'l')
plot(run$eta[2,2,], type = 'l')

# relative abundance ------------------------------------------------------

pars = create_pars(fish_dat)

b_hat = apply(run$beta[,-c(burnin)], 1, mean)
lam_hat = exp(pars$X[,-c(6:11)] %*% b_hat[-c(6:11)])
lam_hat_full = pars$effort*exp(pars$X %*% b_hat)
sum(pars$Y - lam_hat_full)/length(pars$Y)

plot(pars$Y - lam_hat)

res_dat = fish_dat %>% 
  mutate(rel_inten = lam_hat,
         rel_abun = lam_hat_full) %>% 
  left_join(static, by = 'DOW') %>% 
  select(c(DOW, rel_abun, rel_inten, TOTAL_CATCH, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5)) %>% 
  mutate(e_hat = scale(TOTAL_CATCH - rel_abun))

lats = range(res_dat$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(res_dat$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

p1 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = res_dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = rel_inten)) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Intensity")

# ggsave('results/model_plots/intensity.png', p1)

p2 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = res_dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = e_hat)) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red', limits = c(-3,3)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Residuals")

# ggsave('results/model_plots/residuals.png', p2)

p3 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = res_dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = log(rel_abun))) +
  scale_color_gradient(low = 'yellow', high = 'red', name = "Log Relative Abundance") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Log Relative Abundance")

# ggsave('results/model_plots/log_relative_abundance.png', p3)

p1
p2
p3


# effectiveness of gear ---------------------------------------------------

pars = create_pars(fish_dat)

b_hat = apply(run$beta[,-c(burnin)], 1, mean)
lam_hat = pars$effort * exp(pars$X %*% b_hat)
sum(pars$Y - lam_hat)/length(pars$Y)


gear_effect = as.matrix(fish_dat %>% select(GN:GSH))

# gear_ind = 
# fish_dat %>%
#   select(GN:GSH) %>% 
#   pivot_longer(GN:GSH, names_to = 'Gear', values_to = 'Ind') %>% 
#   filter(Ind == 1) %>% 
#   mutate(effect = c(pars$effort * exp(gear_effect %*% b_hat[15:20]))) %>% 
#   select(-Ind) %>% 
#   group_by(Gear) %>% 
#   summarize_all(mean)

fish_dat %>%
  select(GN:GSH) %>% 
  pivot_longer(GN:GSH, names_to = 'Gear', values_to = 'Ind') %>% 
  filter(Ind == 1) %>% 
  mutate(effect = c(exp(gear_effect %*% b_hat[6:11]))) %>% 
  select(-Ind) %>% 
  group_by(Gear) %>% 
  summarize_all(mean)



fish_dat %>%
  select(TOTAL_CATCH, EFFORT, GN:GSH) %>% 
  mutate(CPUE = TOTAL_CATCH/EFFORT) %>% 
  select(-c(TOTAL_CATCH, EFFORT)) %>% 
  pivot_longer(GN:GSH, names_to = 'Gear', values_to = 'Ind') %>% 
  filter(Ind == 1) %>% 
  select(-Ind) %>% 
  group_by(Gear) %>% 
  summarise_all(mean)











# test functions ----------------------------------------------------------

update_w <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  effort = pars$effort
  beta = pars$beta
  w = pars$w
  mu_w = pars$mu_w
  sigma_w = pars$sigma_w
  sig_prop_w = pars$sig_prop_w
  
  # parameters
  n = pars$n
  p = pars$p
  K = pars$K
  
  # indexing
  lake_index = pars$lake_index
  lake_id = pars$lake_id
  n_lakes = pars$n_lakes
  
  
  # w_prop = t(apply(w, 1, function(X) rmvnorm(1, X, sigma_w)))
  w_accept = rep(0, n_lakes)
  
  W_dat_X = W_dat_Y = W_dat_effort = list()
  # set up species random effect
  for(k in 1:K){
    
    W_dat_X[[k]] = data.frame(id = lake_index, X[[k]])
    W_dat_Y[[k]] = data.frame(id = lake_index, Y = Y[[k]])
    W_dat_effort[[k]] = data.frame(id = lake_index, effort = effort[[k]])
    
  }
  
  
  for(i in 1:n_lakes){
    
    lid = lake_id[i]
    like_curr = like_prop = 0
    w_prop_lake = rmvnorm(1, w[i,], sig_prop_w[i]*diag(K))
    
    for(k in 1:K){
      
      X_lake = as.matrix(W_dat_X[[k]] %>% filter(id == lid) %>% select(-id))
      Y_lake = as.matrix(W_dat_Y[[k]] %>% filter(id == lid) %>% select(-id))
      effort_lake = as.matrix(W_dat_effort[[k]] %>% filter(id == lid) %>% select(-id))
      w_curr_lake = w[i,k]
      # w_prop_lake = w_prop[i,k]
      # w_prop_lake = rmvnorm(1, w[i,], sig_prop_w[i]*diag(K))
      
      
      like_curr = like_curr + sum(dpois(Y_lake, lambda = effort_lake*exp(X_lake %*% beta[k,] + w_curr_lake), log = T))
      like_prop = like_prop + sum(dpois(Y_lake, lambda = effort_lake*exp(X_lake %*% beta[k,] + w_prop_lake[k]), log = T))
    }
    
    if((like_prop - like_curr) > log(runif(1))){
      # w[i,] = w_prop[i,]
      w[i,] = w_prop
      w_accept[i] = 1
    }
    
  }
  
  pars$w = w
  pars$w_accept = w_accept
  
  return(pars)
  
  
}
