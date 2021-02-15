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
         DOY_sin = sin(DOY/365 * 2*pi),
         DOY_cos = cos(DOY/365 * 2*pi),
         DOY_sin_semi = sin(DOY/365 * 4*pi),
         DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
  filter(complete.cases(.)) %>% 
  filter(COMMON_NAME != 'white sucker',
         COMMON_NAME != 'smallmouth bass') %>%
  filter(SURVEYDATE >= '1993-01-01') %>% 
  filter(GEAR == 'GN' | GEAR == 'TN') %>%
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
  group_by(GEAR, SURVEYDATE, DOW) %>% 
  mutate(EFFORT = max(EFFORT)) %>%
  ungroup() %>% 
  group_by(DOW, SURVEYDATE) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:DOY), ~ mean(.[!is.na(.) & . != 0])) %>% 
  ungroup() %>% 
  mutate(TN = ifelse(GEAR == 'TN', 1, 0),
         GN = ifelse(GEAR == 'GN', 1, 0)) %>% 
  select(-GEAR) %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW)) %>%
  mutate(DOY = ifelse(is.na(DOY), 0, DOY),
         DOY_sin = ifelse(is.na(DOY_sin), sin(DOY/365 * 2*pi), DOY_sin),
         DOY_cos = ifelse(is.na(DOY_cos), cos(DOY/365 * 2*pi), DOY_cos),
         DOY_sin_semi = ifelse(is.na(DOY_sin_semi), sin(DOY/365 * 4*pi), DOY_sin_semi),
         DOY_cos_semi = ifelse(is.na(DOY_cos_semi), cos(DOY/365 * 4*pi), DOY_cos_semi)) %>% 
  arrange(DOW)


# add temperature 

temp = temp %>%
  select(date, temp_0, MNDOW_ID) %>% # C5 temperature
  rename(SURVEYDATE = date) %>%
  mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>%
  select(-MNDOW_ID) %>% 
  mutate(DOY = yday(SURVEYDATE)) %>% 
  group_by(SURVEYDATE) %>% 
  mutate(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0)) %>% 
  ungroup()

fish_dat = fish_dat %>% 
  inner_join(temp)

rm(all)




# proportion of gear type by species --------------------------------------

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
      select(Int, MAX_DEPTH_FEET:DOY_cos_semi, temp_0) %>% 
      mutate_at(vars(MAX_DEPTH_FEET:LAKE_AREA_GIS_ACRES), ~ ifelse(. == 0, . + 0.001, .)) %>% 
      mutate_at(vars(MAX_DEPTH_FEET:LAKE_AREA_GIS_ACRES), ~ log(.)) %>% 
      mutate_at(vars(DOY_sin:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
      select(-c(LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING, mean.gdd, MAX_DEPTH_FEET, DOY, Int))
    Z = fish_dat %>% 
      filter(COMMON_NAME == levs[k]) %>%
      select(GN)
    # select(GN, GO)
    pars$X[[k]] = as.matrix(tibble(X, Z) %>% 
                              mutate_at(vars(DOY_sin:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>%
                              select(-c(DOY_sin:DOY_cos, DOY_sin_temp:DOY_cos_temp,
                                        DOY_sin_GN:DOY_cos_GN, DOY_sin_temp_GN:DOY_cos_temp_GN)))
  }
  
  for(k in 1:K){
    pars$effort[[k]] = (fish_dat %>% filter(COMMON_NAME == levs[k]) %>% select(EFFORT))$EFFORT
  }
  
  
  # parameters
  pars$n = unlist(lapply(pars$Y, length))
  pars$p = ncol(pars$X[[1]])
  pars$K = K
  
  pars$beta_0 = 0
  pars$beta_0_accept = 0
  pars$beta_0_prior_var = 100
  pars$beta = array(0, dim = c(K, pars$p))
  pars$beta_accept =  array(0, dim = c(K, pars$p))
  pars$beta_prior_var = 100
  
  # hyperpriors
  pars$Sigma_species = diag(K)
  pars$nu_species = K + 2
  pars$Psi_species = 10*diag(K)
  
  pars$omega = matrix(rep(0, K), ncol = K)
  
  # Proposal variances
  pars$sig_prop = array(2, dim = c(K, pars$p))
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
  effort = pars$effort
  omega = pars$omega
  
  # parameters
  n = pars$n
  p = pars$p
  K = pars$K
  
  # beta monitor values
  beta_accept = array(0, dim = c(K, p))
  beta_curr = pars$beta
  sig_prop = pars$sig_prop
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
    like_curr = like_curr + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta_curr[k,] + omega[k]), log = T)) + dnorm(beta_0, 0, sig_prop_beta_0, log = T)
    like_prop = like_prop + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(b_prop + X[[k]] %*% beta_curr[k,] + omega[k]), log = T)) + dnorm(b_prop, 0, sig_prop_beta_0, log = T)
  }
  # if(is.na(like_prop) | is.na(like_curr)){print(b_prop, beta_0, beta_curr)}
  if((like_prop - like_curr) > log(runif(1))){
    beta_0 = b_prop
    beta_0_accept = 1
  }
  
  # all other betas
  for(i in 1:p){
    for(k in 1:K){
      
      beta_prop = beta_curr
      
      b_prop = beta_prop[k,i] = rnorm(1, beta_curr[k,i], sig_prop[k,i])
      b_curr = beta_curr[k,i]
      
      like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta_curr[k,] + omega[k]), log = T)) + dnorm(b_curr, 0, beta_prior_var, log = T)
      like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta_prop[k,] + omega[k]), log = T)) + dnorm(b_prop, 0, beta_prior_var, log = T)
      
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

ll_calc <- function(Y, effort, beta_0, X, beta, omega, K){
  ll = 0
  for(k in 1:K){
    ll = ll + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta[k,] + omega[k]), log = T))
  }
  return(ll)
}

update_omega <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  effort = pars$effort
  beta_0 = pars$beta_0
  beta = pars$beta
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
  logy = ll_calc(Y, effort, beta_0, X, beta, omega, K) + log(runif(1))
  
  # draw initial proposal
  theta = runif(1, 0, 2*pi)
  theta_min = theta - 2*pi
  theta_max = theta
  
  f_proposal = omega * cos(theta) + nu * sin(theta)
  logf = ll_calc(Y, effort, beta_0, X, beta, f_proposal, K)
  
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
      logf = ll_calc(Y, effort, beta_0, X, beta, f_proposal, K)
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

update_proposal_var <- function(pars, beta_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop
  
  bp = beta_accept_post[,,(i-check_num+1):i]
  accept_rate = apply(bp, c(1,2), mean)
  
  sig_prop = ifelse(accept_rate < 0.2, sig_prop*0.9, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.9, sig_prop)
  
  pars$sig_prop = sig_prop
  
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

sampler <- function(fish_dat, nits, check_num = 200){
  
  pars = create_pars(fish_dat)
  
  p = pars$p
  K = pars$K
  
  beta_post = array(NA, dim = c(K, p, nits))
  beta_accept_post = array(NA, dim = c(K, p, nits))
  beta_0_post = array(NA, dim = c(nits))
  beta_0_accept_post = array(NA, dim = c(nits))
  omega_post = array(NA, dim = c(length(pars$omega), nits))
  sigma_species_post = array(NA, dim = c(dim(pars$Sigma_species), nits))
  
  pb <- progress_bar$new(
    format = "  Running [:bar] :percent eta: :eta",
    total = nits, clear = FALSE, width = 60)
  
  for(i in seq(1, nits)){
    
    pars <- update_beta(pars)
    pars <- update_omega(pars)
    pars <- update_sigma_species(pars)
    
    beta_post[,,i] = pars$beta
    beta_accept_post[,,i] = pars$beta_accept
    beta_0_post[i] = pars$beta_0
    beta_0_accept_post[i] = pars$beta_0_accept
    omega_post[,i] = pars$omega
    sigma_species_post[,,i] = pars$Sigma_species
    
    
    if(i %in% seq(0, nits-1, by = check_num)){
      pars <- update_proposal_var(pars, beta_accept_post, i, check_num)
      pars <- update_proposal_var_beta_0(pars, beta_0_accept_post, i, check_num)
    }
    
    
    pb$tick()
    
  }
  
  return(list(beta = beta_post,
              beta_accept = beta_accept_post,
              sig_prop = pars$sig_prop,
              beta_0 = beta_0_post,
              beta_0_accept = beta_0_accept_post,
              sig_prop_beta_0 = pars$sig_prop_beta_0,
              omega = omega_post,
              sig_prop_omega = pars$sig_prop_omega,
              sigma_species = sigma_species_post))
  
}

nits = 20000
burnin = 1:10000
thin = seq(25000, 50000, by = 100)

run = sampler(fish_dat, nits, check_num = 100)

pars = create_pars(fish_dat)
nms = colnames(pars$X[[1]])

plot(run$beta_0[-c(burnin)], type = 'l')


par(mfrow = c(3,5))
for(i in 1:13){
  plot(run$beta[6, i,-c(burnin)], type = 'l', main = nms[i])
}


fish_names = str_replace(pars$fish_names, " ", "_")
for(j in 1:6){
  plt_name = paste0('/Users/joshuanorth/Desktop/fish_plots/', fish_names[j], '.png')
  png(plt_name, width = 1400, height = 800)
  par(mfrow = c(3,5))
  for(i in 1:13){
    plot(run$beta[j, i,-c(burnin)], type = 'l', main = nms[i])
  }
  dev.off()
}

apply(run$beta_accept[,,-c(burnin)], c(1,2), mean)
round(run$sig_prop, 3)

cbind(nms, t(apply(run$beta, c(1,2), mean)))

apply(run$omega, 1, mean)

plot(run$omega[1,-c(burnin)], type = 'l')
plot(run$omega[2,-c(burnin)], type = 'l')
plot(run$omega[3,-c(burnin)], type = 'l')
plot(run$omega[4,-c(burnin)], type = 'l')
plot(run$omega[5,-c(burnin)], type = 'l')
plot(run$omega[6,-c(burnin)], type = 'l')



plt_name = paste0('/Users/joshuanorth/Desktop/fish_plots/crappie_int.png')
png(plt_name, width = 1400, height = 800)
par(mfrow = c(2,2))
plot(run$beta[1, 1,-c(burnin)], type = 'l', main = 'Intercept')
plot(run$eta[1,-c(burnin)], type = 'l', main = 'Random Effect')
plot(run$beta[1, 1,-c(burnin)] + run$eta[1,-c(burnin)], type = 'l', main = 'Intercept + Random Effect')
dev.off()

fish = 6
par(mfrow = c(2,2))
plot(run$beta[fish, 1,-c(burnin)], type = 'l', main = 'Intercept')
plot(run$eta[fish,-c(burnin)], type = 'l', main = 'Random Effect')
plot(run$beta[fish, 1,-c(burnin)] + run$eta[fish,-c(burnin)], type = 'l', main = 'Intercept + Random Effect')



# relative abundance ------------------------------------------------------

fish_dat %>%
  select(COMMON_NAME, TOTAL_CATCH) %>%
  group_by(COMMON_NAME) %>%
  summarize(mean(TOTAL_CATCH), sd(TOTAL_CATCH), median(TOTAL_CATCH)) #range(TOTAL_CATCH), 



pars = create_pars(fish_dat)
nms = colnames(pars$X[[1]])

b_hat = apply(run$beta[,,-c(burnin)], c(1,2), mean)
colnames(b_hat) = nms

omega_hat = apply(run$omega[,-c(burnin)], c(1), mean)

int_inds = 1:2
int_b = b_hat[,int_inds]

lam_hat = list()
for(k in 1:pars$K){
  lam_hat[[k]] = exp(pars$X[[k]][,int_inds] %*% int_b[k,] + omega_hat[k])
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
              aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = log(Abundance)),
              width = 0.2, height = 0.2) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Relative Log Abundance") +
  facet_wrap(~Fish, 
             labeller = labeller(Fish = c('black_crappie' = 'Black Crappie',
                                          'bluegill' = 'Bluegill',
                                          'northern_pike' = 'Northern Pike',
                                          'yellow_perch' = 'Yellow Perch',
                                          'largemouth_bass' = 'Largemouth Bass',
                                          'walleye' = 'Walleye')))
# ggsave('results/non_spatial_results/relative_abun.png')



ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_jitter(data = rel_abun %>% filter(Fish == 'black_crappie'), 
              aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance),
              width = 0.2, height = 0.2) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Relative Abundance - Black Crappie")
# ggsave('results/non_spatial_results/relative_abun_crappie.png')



ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_jitter(data = rel_abun %>% filter(Fish == 'bluegill'), 
              aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance),
              width = 0.2, height = 0.2) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Relative Abundance - Bluegill")
# ggsave('results/non_spatial_results/relative_abun_bluegill.png')


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
b_0 = mean(run$beta_0[-c(burnin)])
colnames(b_hat) = nms

# alpha_inds = c(1, 4:11)
alpha_inds = c(3:13)
alpha_b = b_hat[,alpha_inds]

b_hat_lower = apply(run$beta[,,-c(burnin)], c(1,2), quantile, probs = c(0.025))
alpha_b_lower = b_hat_lower[,alpha_inds]

b_hat_upper = apply(run$beta[,,-c(burnin)], c(1,2), quantile, probs = c(0.975))
alpha_b_upper = b_hat_upper[,alpha_inds]


fnames = c('crappie', 'bluegill', 'bass', 'pike', 'walleye', 'perch')
lake = '38025600'
# 06015200 78002500 34007900 46010900
lake = '06015200'

temp_raw = read_rds('data/daily_degree_days_MN_lakes.rds') %>% ungroup()

temp_raw = temp_raw %>%
  select(date, temp_0, MNDOW_ID) %>% # C5 temperature
  rename(SURVEYDATE = date) %>%
  mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>%
  select(-MNDOW_ID) %>% 
  mutate(DOY = yday(SURVEYDATE)) %>% 
  group_by(SURVEYDATE) %>% 
  mutate(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0)) %>% 
  ungroup()

temp %>% 
  filter(year(SURVEYDATE) == 2005) %>% 
  group_by(SURVEYDATE) %>% 
  summarise(mean(temp_0)) %>% 
  ungroup()
  filter(DOW == '06015200')


d_plot = function(alpha_b, fish_names, year, fish_dat){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  # tC = temp %>%
  #   filter(year(SURVEYDATE) == year) %>%
  #   group_by(SURVEYDATE) %>%
  #   summarize(temp_0 = mean(temp_0, na.rm = TRUE)) %>%
  #   ungroup() %>%
  #   select(temp_0, SURVEYDATE) %>%
  #   mutate(DOY = seq(1, length(temp_0)),
  #          temp_0 = (temp_0 - mean(temp_0))/sd(temp_0))
  
  tC = temp %>%
    filter(year(SURVEYDATE) == year) %>%
    group_by(SURVEYDATE) %>%
    summarize(temp_0 = mean(temp_0, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(DOY = seq(1, length(temp_0)))
  
  # tC = temp_raw %>%
  #   filter(year(date) == year) %>% 
  #   select(date, temp_0, MNDOW_ID) %>% # C5 temperature
  #   rename(SURVEYDATE = date) %>%
  #   mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>%
  #   select(-MNDOW_ID) %>% 
  #   mutate(DOY = yday(SURVEYDATE)) %>% 
  #   group_by(SURVEYDATE) %>% 
  #   summarise(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0),
  #             DOY = mean(DOY, na.rm = T)) %>% 
  #   ungroup() %>% 
  #   filter(DOY != 366)
  # 
  # 
  # tC = temp %>%
  #   filter(DOW == lake) %>%
  #   filter(year(SURVEYDATE) == year) %>%
  #   filter(DOY != 366) %>% 
  #   select(-DOW)
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == year) %>% 
    select(SURVEYDATE, DOY_sin:DOY_cos_semi, GN, temp_0) %>% 
    mutate_at(vars(DOY_sin:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin:DOY_cos_semi, temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(DOY_sin:DOY_cos, DOY_sin_temp:DOY_cos_temp,
              DOY_sin_GN:DOY_cos_GN, DOY_sin_temp_GN:DOY_cos_temp_GN)) %>% 
    relocate(GN, .after = DOY_cos_semi_temp)
  
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:TN) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == year) %>% 
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
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  GN_df = year_select %>% 
    filter(GN == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 1, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  TN = as.matrix(TN_df) %*% t(alpha_b)
  GN = as.matrix(GN_df) %*% t(alpha_b)
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  cnames = c(paste0("TN_",fish_names),
             paste0("GN_",fish_names))
  
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
    mutate_at(vars(GN:TN), ~ if_else(. == 1, SURVEYDATE, NaN)) %>% 
    select(-sample_day) %>% 
    arrange(SURVEYDATE) %>% 
    rename_at(vars(-SURVEYDATE), ~ paste0(., '_survey'))
  
  return(tibble(as_tibble(TN), as_tibble(GN), date_select))
}

year = 2000
mean = d_plot(alpha_b, fnames, year, fish_dat)
upper = d_plot(alpha_b_upper, fnames, year, fish_dat)
lower = d_plot(alpha_b_lower, fnames, year, fish_dat)


mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness") %>% 
  left_join(upper %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_upper"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(lower %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_lower"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>%
  left_join(mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
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

# ggsave(paste0('results/non_spatial_results/catchability_', year, '.png'), width = 10, height = 6)



d_plot = function(alpha_b, alpha_b_upper, alpha_b_lower, fish_names, year, fish_dat){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  tC = temp %>%
    filter(year(SURVEYDATE) == year) %>%
    group_by(SURVEYDATE) %>%
    summarize(temp_0 = mean(temp_0, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(DOY = seq(1, length(temp_0)))
  
  tC_lower = temp %>%
    filter(year(SURVEYDATE) == year) %>%
    group_by(SURVEYDATE) %>%
    summarize(temp_0 = quantile(temp_0, probs = 0.025, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(DOY = seq(1, length(temp_0)))
  
  tC_upper = temp %>%
    filter(year(SURVEYDATE) == year) %>%
    group_by(SURVEYDATE) %>%
    summarize(temp_0 = quantile(temp_0, probs = 0.975, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(DOY = seq(1, length(temp_0)))
  
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == year) %>% 
    select(SURVEYDATE, DOY_sin:DOY_cos_semi, GN, temp_0) %>% 
    mutate_at(vars(DOY_sin:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin:DOY_cos_semi, temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(DOY_sin:DOY_cos, DOY_sin_temp:DOY_cos_temp,
              DOY_sin_GN:DOY_cos_GN, DOY_sin_temp_GN:DOY_cos_temp_GN)) %>% 
    relocate(GN, .after = DOY_cos_semi_temp)
  
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:TN) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == year) %>% 
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
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  GN_df = year_select %>% 
    filter(GN == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 1, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  TN_df_lower = year_select %>% 
    filter(GN != 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC_lower, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 0, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  GN_df_lower = year_select %>% 
    filter(GN == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC_lower, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 1, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  TN_df_upper = year_select %>% 
    filter(GN != 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC_upper, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 0, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  GN_df_upper = year_select %>% 
    filter(GN == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC_upper, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 1, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  TN = as.matrix(TN_df) %*% t(alpha_b)
  GN = as.matrix(GN_df) %*% t(alpha_b)
  
  TN_lower = as.matrix(TN_df_lower) %*% t(alpha_b_lower)
  GN_lower = as.matrix(GN_df_lower) %*% t(alpha_b_lower)
  
  TN_upper = as.matrix(TN_df_upper) %*% t(alpha_b_upper)
  GN_upper = as.matrix(GN_df_upper) %*% t(alpha_b_upper)
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  colnames(TN_lower) = paste0("TN_",fish_names)
  colnames(GN_lower) = paste0("GN_",fish_names)
  
  colnames(TN_upper) = paste0("TN_",fish_names)
  colnames(GN_upper) = paste0("GN_",fish_names)
  
  cnames = c(paste0("TN_",fish_names),
             paste0("GN_",fish_names))
  
  date_select = year_select %>% 
    distinct(SURVEYDATE) %>% 
    mutate(sample_day = SURVEYDATE) %>% 
    right_join(tC, by = c('SURVEYDATE')) %>% 
    left_join(gear_ind, by = c('SURVEYDATE')) %>% 
    select(-c(temp_0, DOY)) %>% 
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

year = 2000

plt_dat = d_plot(alpha_b, alpha_b_upper, alpha_b_lower, fnames, year, fish_dat)


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
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(year, '-03-01')), ymd(paste0(year, '-12-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
  ggtitle(paste0('Gear Type Effectiveness for ', year)) +
  xlab('Date') +
  ylab('Relative Effectiveness')



# posterior gear type effectiveness ---------------------------------------

fnames = c('crappie', 'bluegill', 'bass', 'pike', 'walleye', 'perch')
lake = '38025600'
# 06015200 78002500 34007900 46010900
lake = '06015200'

betas = run$beta[,-c(1,2),-burnin]

d_plot = function(betas, fish_names, year, fish_dat){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  tC = temp %>%
    filter(year(SURVEYDATE) == year) %>%
    group_by(SURVEYDATE) %>%
    summarize(temp_0 = mean(temp_0, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(DOY = seq(1, length(temp_0))) %>% 
    filter(DOY != 366)
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == year) %>% 
    select(SURVEYDATE, DOY_sin:DOY_cos_semi, GN, temp_0) %>% 
    mutate_at(vars(DOY_sin:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin:DOY_cos_semi, temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(DOY_sin:DOY_cos, DOY_sin_temp:DOY_cos_temp,
              DOY_sin_GN:DOY_cos_GN, DOY_sin_temp_GN:DOY_cos_temp_GN)) %>% 
    relocate(GN, .after = DOY_cos_semi_temp)
  
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:TN) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == year) %>% 
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
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  GN_df = year_select %>% 
    filter(GN == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 1, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  TN = as.matrix(TN_df)
  GN = as.matrix(GN_df)
  
  n_runs = dim(betas)[3]
  
  post_stores_TN = post_stores_GN = array(NA, dim = c(nrow(tC), length(fnames), n_runs))
  
  for(i in 1:n_runs){
    post_stores_TN[,,i] = TN %*% t(betas[,,i])
    post_stores_GN[,,i] = GN %*% t(betas[,,i])
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
    select(-c(temp_0, DOY)) %>% 
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

year = 2000
plt_dat = d_plot(betas, fnames, year, fish_dat)

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
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(year, '-03-01')), ymd(paste0(year, '-12-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
  ggtitle(paste0('Gear Type Effectiveness for ', year)) +
  xlab('Date') +
  ylab('Relative Effectiveness')

# ggsave(paste0('results/non_spatial_results/catchability_', year, '.png'), width = 10, height = 6)


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



d_plot_lake = function(betas, fish_names, year, fish_dat, lake_dow){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  tC = temp %>%
    filter(year(SURVEYDATE) == year) %>%
    filter(DOW == lake_dow) %>% 
    filter(DOY != 366) %>% 
    select(-DOW)
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == year) %>% 
    select(SURVEYDATE, DOY_sin:DOY_cos_semi, GN, temp_0) %>% 
    mutate_at(vars(DOY_sin:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin:DOY_cos_semi, temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(DOY_sin:DOY_cos, DOY_sin_temp:DOY_cos_temp,
              DOY_sin_GN:DOY_cos_GN, DOY_sin_temp_GN:DOY_cos_temp_GN)) %>% 
    relocate(GN, .after = DOY_cos_semi_temp)
  
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:TN) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == year) %>% 
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
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  GN_df = year_select %>% 
    filter(GN == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 1, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  TN = as.matrix(TN_df)
  GN = as.matrix(GN_df)
  
  n_runs = dim(betas)[3]
  
  post_stores_TN = post_stores_GN = array(NA, dim = c(nrow(tC), length(fnames), n_runs))
  
  for(i in 1:n_runs){
    post_stores_TN[,,i] = TN %*% t(betas[,,i])
    post_stores_GN[,,i] = GN %*% t(betas[,,i])
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
    select(-c(temp_0, DOY)) %>% 
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

betas = run$beta[,-c(1,2),-burnin]
year = 2015
lake_dow = 53002800 # south
# lake_dow = 16001900 # north
plt_dat = d_plot_lake(betas, fnames, year, fish_dat, lake_dow)


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
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(year, '-03-01')), ymd(paste0(year, '-12-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
  ggtitle(paste0('Gear Type Effectiveness for ', year)) +
  xlab('Date') +
  ylab('Relative Effectiveness')

# ggsave(paste0('results/non_spatial_results/catchability_north_lake', year, '.png'), width = 10, height = 6)
ggsave(paste0('results/non_spatial_results/catchability_south_lake', year, '.png'), width = 10, height = 6)


# compare north and south lakes -------------------------------------------


# 16001900
lake_pos %>% 
  filter(LAKE_CENTER_LAT_DD5 > quantile(LAKE_CENTER_LAT_DD5, probs = 0.975)) %>% 
  filter(LAKE_CENTER_LONG_DD5 > quantile(LAKE_CENTER_LONG_DD5, probs = 0.975))

# 53002800
lake_pos %>% 
  filter(LAKE_CENTER_LAT_DD5 < quantile(LAKE_CENTER_LAT_DD5, probs = 0.025)) %>% 
  filter(LAKE_CENTER_LONG_DD5 < quantile(LAKE_CENTER_LONG_DD5, probs = 0.025))


d_plot_lake_year = function(betas, fish_names, year, fish_dat, lake_dow){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  tC = temp %>%
    filter(DOW == lake_dow) %>% 
    filter(DOY != 366) %>% 
    group_by(DOY) %>% 
    summarize(temp_0 = mean(temp_0)) %>% 
    ungroup() %>% 
    mutate(SURVEYDATE = as.Date(DOY, origin = paste0(year,"-01-01")))
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == year) %>% 
    select(SURVEYDATE, DOY_sin:DOY_cos_semi, GN, temp_0) %>% 
    mutate_at(vars(DOY_sin:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin:DOY_cos_semi, temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(DOY_sin:DOY_cos, DOY_sin_temp:DOY_cos_temp,
              DOY_sin_GN:DOY_cos_GN, DOY_sin_temp_GN:DOY_cos_temp_GN)) %>% 
    relocate(GN, .after = DOY_cos_semi_temp)
  
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:TN) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == year) %>% 
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
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  GN_df = year_select %>% 
    filter(GN == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 1, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi)) %>%
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi), .funs = list(temp = ~.*temp_0)) %>% 
    mutate_at(vars(DOY_sin_semi:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  TN = as.matrix(TN_df)
  GN = as.matrix(GN_df)
  
  n_runs = dim(betas)[3]
  
  post_stores_TN = post_stores_GN = array(NA, dim = c(nrow(tC), length(fnames), n_runs))
  
  for(i in 1:n_runs){
    post_stores_TN[,,i] = TN %*% t(betas[,,i])
    post_stores_GN[,,i] = GN %*% t(betas[,,i])
  }
  
  TN = apply(post_stores_TN, c(1,2), mean)
  GN = apply(post_stores_GN, c(1,2), mean)
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  
  date_select = year_select %>% 
    distinct(SURVEYDATE) %>% 
    mutate(sample_day = SURVEYDATE) %>% 
    right_join(tC, by = c('SURVEYDATE')) %>% 
    left_join(gear_ind, by = c('SURVEYDATE')) %>% 
    select(-c(temp_0, DOY)) %>% 
    mutate_at(vars(GN:TN), ~ if_else(. == 1, SURVEYDATE, NaN)) %>% 
    select(-sample_day) %>% 
    arrange(SURVEYDATE) %>% 
    rename_at(vars(-SURVEYDATE), ~ paste0(., '_survey'))
  
  mean = tibble(as_tibble(TN), as_tibble(GN), date_select)
  
  return(mean)
}

betas = run$beta[,-c(1,2),-burnin]
year = 2014
lake_dow_south = 53002800 # south
lake_dow_north = 16001900 # north
plt_dat_south = d_plot_lake_year(betas, fnames, year, fish_dat, lake_dow_south)
plt_dat_north = d_plot_lake_year(betas, fnames, year, fish_dat, lake_dow_north)


plt_dat_south %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_South") %>% 
  left_join(plt_dat_north %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_North"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  pivot_longer(cols = -c(SURVEYDATE:Fish), names_to = c("Type", "Lake"), names_sep = "_", values_to = "Effectiveness") %>% 
  select(-Type) %>% 
  left_join(plt_dat_south %>%
            select(SURVEYDATE:TN_survey) %>% 
            pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-sample_day) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(aes(linetype = Lake), size = 1.2) +
  facet_wrap(~ Gear) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(year, '-05-01')), ymd(paste0(year, '-10-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 10)) +
  ggtitle(paste0('Gear Type Effectiveness for a North and South lake')) +
  xlab('Date') +
  ylab('Relative Effectiveness')

ggsave(paste0('results/non_spatial_results/catchability_lake_comp.png'), width = 10, height = 6)












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
  filter(COMMON_NAME != 'white sucker',
         COMMON_NAME != 'smallmouth bass') %>%
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