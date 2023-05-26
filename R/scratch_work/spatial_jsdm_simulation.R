

# load libraries ----------------------------------------------------------

library(tidyverse)
library(rstan)
library(fields)
library(Matrix)
library(lubridate)
library(stringr)
library(LaplacesDemon)



# load data ---------------------------------------------------------------


center_temp = read_csv('data/center_temp.csv')
center_temp = center_temp %>%
  mutate(DOW = as.factor(DOW)) %>%
  arrange(DOW, year)


fish_dat = read_csv("data/fish_dat.csv")
fish_dat = fish_dat %>%
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME)) %>%
  arrange(DOW, year, COMMON_NAME)


# reduce the size of the data
# set.seed(1)
# lselect = fish_dat %>% 
#   select(DOW) %>% 
#   distinct() %>% 
#   sample_n(200) %>% 
#   droplevels()


lselect = fish_dat %>% 
  mutate(cntr = 1) %>% 
  group_by(DOW) %>% 
  summarise(cnt = sum(cntr)) %>% 
  ungroup() %>% 
  arrange(desc(cnt)) %>% 
  slice_head(n = 300) %>% 
  select(DOW) %>% 
  droplevels()




fish_dat_sub = fish_dat %>% 
  right_join(lselect) %>% 
  droplevels()

center_temp_sub = center_temp %>% 
  right_join(lselect) %>% 
  droplevels()


fish_dat %>% 
  group_by(year) %>% 
  mutate(cnt = 1) %>% 
  summarise(cnt = sum(cnt)/12)



# select covariates from data
# mean_covs = colnames(fish_dat)[c(7, 9, 23, 24, 14:15)]
# temporal_covs = colnames(fish_dat)[c(23, 24)]
# mean_covs_log = colnames(fish_dat)[c(7, 9)]
# mean_covs_logit = colnames(fish_dat)[c(14:15)]
# catch_covs = colnames(fish_dat)[c(25, 27:30)] # temp doy interaction
# gear_types = colnames(fish_dat)[c(21, 22)]

mean_covs = colnames(fish_dat)[c(7, 9, 24)]
temporal_covs = colnames(fish_dat)[c(24)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = NULL
catch_covs = colnames(fish_dat)[c(25, 27:30)]
gear_types = colnames(fish_dat)[c(21, 22)]

create_pars <- function(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, center_temp){
  
  fish_dat = fish_dat %>% 
    arrange(DOW, SURVEYDATE, COMMON_NAME, GN)
  
  K = nlevels(fish_dat$COMMON_NAME)
  levs = levels(fish_dat$COMMON_NAME)
  lake_index = (fish_dat %>% filter(COMMON_NAME == levs[1]))$DOW
  lake_id = levels(fish_dat$DOW)
  n_lakes = length(levels(fish_dat$DOW))
  n_obs = (fish_dat %>% filter(COMMON_NAME == levs[1]) %>% select(TOTAL_CATCH))$TOTAL_CATCH
  
  pars = list()
  
  
  pars$Y = fish_dat %>% 
    arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
    select(DOW, SURVEYDATE, COMMON_NAME, GN, TOTAL_CATCH) %>% 
    pivot_wider(names_from = COMMON_NAME, values_from = TOTAL_CATCH) %>% 
    select(-c(DOW, SURVEYDATE, GN)) %>% 
    as.matrix()
  
  pars$effort = fish_dat %>% 
    arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
    select(DOW, SURVEYDATE, COMMON_NAME, GN, EFFORT) %>% 
    pivot_wider(names_from = COMMON_NAME, values_from = EFFORT) %>% 
    select(-c(DOW, SURVEYDATE, GN)) %>% 
    as.matrix()
  
  X = fish_dat %>% 
    arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
    filter(COMMON_NAME == levs[1]) %>% 
    select(all_of(mean_covs), GN) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.))
  # mutate(secchi = secchi - mean(secchi))
  
  Z = fish_dat %>% 
    arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
    filter(COMMON_NAME == levs[1]) %>%
    select(all_of(catch_covs), GN) %>% 
    mutate(TN = 1) %>% 
    relocate(TN) %>% 
    mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN))
  
  temp_obs_lakes = center_temp %>% 
    filter(DOW %in% unique(fish_dat$DOW)) %>% 
    group_by(DOW) %>% 
    mutate(dow = weekdays(SURVEYDATE)) %>% 
    filter(dow == "Monday") %>% 
    select(-dow) %>% 
    ungroup()
  
  TN = temp_obs_lakes %>% 
    mutate(TN = 1,
           DOY = yday(SURVEYDATE)) %>% 
    relocate(temp_0, .after = DOY) %>% 
    mutate(DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/121 * 2*pi),
           DOY_cos_semi = cos(DOY/121 * 2*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0,
           GN = 0) %>% 
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY, DOW, year))
  
  GN = temp_obs_lakes %>% 
    mutate(TN = 1,
           DOY = yday(SURVEYDATE)) %>% 
    relocate(temp_0, .after = DOY) %>% 
    mutate(DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/121 * 2*pi),
           DOY_cos_semi = cos(DOY/121 * 2*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0,
           GN = 1) %>% 
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY, DOW, year))
  
  pars$X = as.matrix(X)
  pars$Z = as.matrix(Z %>% select(-c(TN, GN)))
  # pars$Zstar = rbind(TN, GN) %>% select(-c(TN, GN)) %>% as.matrix()
  pars$Zstar = as.matrix(Z %>% select(-c(TN, GN)))
  pars$y_vec = c(pars$Y)
  
  pars$Nstar = dim(pars$Zstar)[1]
  
  # parameters
  pars$N = nrow(pars$Y)
  pars$p_beta = ncol(pars$X)
  pars$p_phi = ncol(pars$Z)
  pars$K = K
  
  # spatial parameters
  spat_dat = fish_dat %>% 
    distinct(DOW, .keep_all = T) %>% 
    select(DOW, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING)
  
  d = rdist(cbind(spat_dat$LAKE_CENTER_UTM_EASTING, spat_dat$LAKE_CENTER_UTM_NORTHING))/1000
  phi = 10
  # pars$Sigma_spatial = Matrix(exp(-d/phi))
  pars$Sigma_spatial = exp(-d/phi)
  pars$Sigma_spatial_inv = solve(pars$Sigma_spatial)
  pars$up_chol_spatial = t(chol(exp(-d/phi)))
  
  # indexing
  pars$lake_index = lake_index
  # pars$each_lake = table(lake_index)
  pars$lake_id = lake_id
  pars$n_lakes = n_lakes
  pars$fish_names = levels(fish_dat$COMMON_NAME)
  
  
  pars$each_lake = fish_dat %>% 
    filter(COMMON_NAME == levs[1]) %>% 
    select(DOW) %>% 
    mutate(ind = 1,
           id = 1:n()) %>% 
    pivot_wider(names_from = DOW, values_from = ind, values_fill = 0) %>% 
    dplyr::select(-id) %>% 
    as.matrix() %>% 
    unname()
  
  return(pars)
  
}

# create_pars <- function(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, center_temp){
#   
#   fish_dat = fish_dat %>% 
#     arrange(DOW, SURVEYDATE, COMMON_NAME, GN)
#   
#   K = nlevels(fish_dat$COMMON_NAME)
#   levs = levels(fish_dat$COMMON_NAME)
#   lake_index = (fish_dat %>% filter(COMMON_NAME == levs[1]))$DOW
#   lake_id = levels(fish_dat$DOW)
#   n_lakes = length(levels(fish_dat$DOW))
#   n_obs = (fish_dat %>% filter(COMMON_NAME == levs[1]) %>% select(TOTAL_CATCH))$TOTAL_CATCH
#   
#   pars = list()
#   
#   
#   pars$Y = fish_dat %>% 
#     arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
#     select(DOW, SURVEYDATE, COMMON_NAME, GN, TOTAL_CATCH) %>% 
#     pivot_wider(names_from = COMMON_NAME, values_from = TOTAL_CATCH) %>% 
#     select(-c(DOW, SURVEYDATE, GN)) %>% 
#     as.matrix()
#   
#   pars$effort = fish_dat %>% 
#     arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
#     select(DOW, SURVEYDATE, COMMON_NAME, GN, EFFORT) %>% 
#     pivot_wider(names_from = COMMON_NAME, values_from = EFFORT) %>% 
#     select(-c(DOW, SURVEYDATE, GN)) %>% 
#     as.matrix()
#   
#   X = fish_dat %>% 
#     arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
#     filter(COMMON_NAME == levs[1]) %>% 
#     select(all_of(mean_covs)) %>% 
#     mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
#     mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
#     mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
#     mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.))
#   # mutate(secchi = secchi - mean(secchi))
#   
#   Z = fish_dat %>% 
#     arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
#     filter(COMMON_NAME == levs[1]) %>%
#     select(all_of(catch_covs), GN) %>% 
#     mutate(TN = 1) %>% 
#     relocate(TN) %>% 
#     mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN)) %>% 
#     select(-TN)
#   
#   temp_obs_lakes = center_temp %>% 
#     filter(DOW %in% unique(fish_dat$DOW)) %>% 
#     group_by(DOW) %>% 
#     mutate(dow = weekdays(SURVEYDATE)) %>% 
#     filter(dow == "Monday") %>% 
#     select(-dow) %>% 
#     ungroup()
#   
#   TN = temp_obs_lakes %>% 
#     mutate(TN = 1,
#            DOY = yday(SURVEYDATE)) %>% 
#     relocate(temp_0, .after = DOY) %>% 
#     mutate(DOY = yday(SURVEYDATE),
#            DOY_sin_semi = sin(DOY/121 * 2*pi),
#            DOY_cos_semi = cos(DOY/121 * 2*pi),
#            DOY_sin_semi_temp = DOY_sin_semi * temp_0,
#            DOY_cos_semi_temp = DOY_cos_semi * temp_0,
#            GN = 0) %>% 
#     mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
#     select(-c(SURVEYDATE, DOY, DOW, year))
#   
#   GN = temp_obs_lakes %>% 
#     mutate(TN = 1,
#            DOY = yday(SURVEYDATE)) %>% 
#     relocate(temp_0, .after = DOY) %>% 
#     mutate(DOY = yday(SURVEYDATE),
#            DOY_sin_semi = sin(DOY/121 * 2*pi),
#            DOY_cos_semi = cos(DOY/121 * 2*pi),
#            DOY_sin_semi_temp = DOY_sin_semi * temp_0,
#            DOY_cos_semi_temp = DOY_cos_semi * temp_0,
#            GN = 1) %>% 
#     mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
#     select(-c(SURVEYDATE, DOY, DOW, year))
#   
#   pars$X = as.matrix(X)
#   pars$Z = as.matrix(Z)
#   pars$Zstar = rbind(TN, GN) %>% select(-TN) %>% as.matrix()
#   pars$y_vec = c(pars$Y)
#   
#   pars$Nstar = dim(pars$Zstar)[1]
#   
#   # parameters
#   pars$N = nrow(pars$Y)
#   pars$p_beta = ncol(pars$X)
#   pars$p_phi = ncol(pars$Z)
#   pars$K = K
#   
#   # spatial parameters
#   spat_dat = fish_dat %>% 
#     distinct(DOW, .keep_all = T) %>% 
#     select(DOW, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING)
#   
#   d = rdist(cbind(spat_dat$LAKE_CENTER_UTM_EASTING, spat_dat$LAKE_CENTER_UTM_NORTHING))/1000
#   phi = 10
#   # pars$Sigma_spatial = Matrix(exp(-d/phi))
#   pars$Sigma_spatial = exp(-d/phi)
#   pars$Sigma_spatial_inv = solve(pars$Sigma_spatial)
#   pars$up_chol_spatial = t(chol(exp(-d/phi)))
#   
#   # indexing
#   pars$lake_index = lake_index
#   # pars$each_lake = table(lake_index)
#   pars$lake_id = lake_id
#   pars$n_lakes = n_lakes
#   pars$fish_names = levels(fish_dat$COMMON_NAME)
#   
#   
#   pars$each_lake = fish_dat %>% 
#     filter(COMMON_NAME == levs[1]) %>% 
#     select(DOW) %>% 
#     mutate(ind = 1,
#            id = 1:n()) %>% 
#     pivot_wider(names_from = DOW, values_from = ind, values_fill = 0) %>% 
#     dplyr::select(-id) %>% 
#     as.matrix() %>% 
#     unname()
#   
#   return(pars)
#   
# }


obs_only_temp = center_temp_sub %>% 
  right_join(fish_dat_sub %>% 
               select(SURVEYDATE, temp_0, DOW, year))



dat = create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, center_temp)

# dat_train = create_pars(fish_dat_sub %>% 
#                           mutate(year = year(SURVEYDATE)) %>% 
#                           filter(year != 2008 | year != 2018) %>% 
#                           droplevels(), 
#                         mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, 
#                         center_temp_sub %>% 
#                           mutate(year = year(SURVEYDATE)) %>% 
#                           filter(year != 2008 | year != 2018) %>% 
#                           droplevels())
# dat_train = create_pars(fish_dat_sub %>% 
#                           mutate(year = year(SURVEYDATE)) %>% 
#                           filter(year != 2008 | year != 2018) %>% 
#                           droplevels(), 
#                         mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, 
#                         center_temp_sub)
# 
# dat_test = create_pars(fish_dat_sub %>% 
#                          mutate(year = year(SURVEYDATE)) %>% 
#                          filter(year == 2008 | year == 2018) %>% 
#                          droplevels(), 
#                        mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, 
#                        center_temp_sub %>% 
#                          mutate(year = year(SURVEYDATE)) %>% 
#                          filter(year == 2008 | year == 2018) %>% 
#                          droplevels())

dat_train = create_pars(fish_dat_sub %>% 
                          mutate(year = year(SURVEYDATE)) %>% 
                          filter(year != 2008 | year != 2018) %>% 
                          droplevels(), 
                        mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, 
                        obs_only_temp)

dat_test = create_pars(fish_dat_sub %>% 
                         mutate(year = year(SURVEYDATE)) %>% 
                         filter(year == 2008 | year == 2018) %>% 
                         droplevels(), 
                       mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, 
                       obs_only_temp %>% 
                         mutate(year = year(SURVEYDATE)) %>% 
                         filter(year == 2008 | year == 2018) %>% 
                         droplevels())

str(dat)
str(dat_train)
str(dat_test)


# simulate data ----------------------------------------------------------


gen_Y = function(beta_0, beta, phi, Sigma_species, tau, omega, dat){
  
  OMEGA = dat$each_lake %*% omega
  
  # relative abundance
  B0 = matrix(rep(beta_0, dat$N), dat$N, dat$K, byrow = T)
  gamma = B0 + dat$X %*% beta + OMEGA
  # gamma = dat$X %*% beta + OMEGA
  
  # catchability
  theta_star = dat$Z %*% phi
  theta_star_star = dat$Zstar %*% phi
  # theta = theta_star - mean(theta_star_star) - var(c(theta_star_star))/2
  mux = mean(exp(theta_star_star))^2
  sx = var(exp(c(theta_star_star)))
  mu = log(mux/sqrt(mux+sx))
  s = log(1 + sx/mux)
  
  theta = theta_star - mu - s/2
  # mean(exp(theta))
  
  
  lambda = exp(log(dat$effort) + theta + gamma)
  y = rpois(dat$N*dat$K, c(lambda))
  Y = matrix(y, dat$N, dat$K)
  
  return(list(y_vec = y, Y = Y))
  
}


out = read_rds("/Users/jsnow/Desktop/stan_lake_random_effect.rds")
chains = extract(out, permuted = T)

# beta_0 = apply(chains$beta_0, 2, mean)
# beta = apply(chains$beta, c(2,3), mean)
# # phi = apply(chains$phi, c(2,3), mean)
# phi = chains$phi[900,,]
# omega = apply(chains$omega, c(2,3), mean)
# Sigma_species = apply(chains$Sigma_species, c(2,3), mean)
# tau = apply(chains$tau, 2, mean)

set.seed(1)
beta_0 = apply(chains$beta_0, 2, mean)
# beta = apply(chains$beta, c(2,3), mean)[c(1,2,4),]
beta = rbind(apply(chains$beta, c(2,3), mean)[c(1,2,4),], runif(6))
# tn = runif(6)
# gn = 1-tn
# phi = rbind(tn, gn)
x = rnorm(dat$p_phi * dat$K)
phi = matrix(x - mean(x), nrow = dat$p_phi)
sum(phi)
omega = apply(chains$omega, c(2,3), mean)
Sigma_species = apply(chains$Sigma_species, c(2,3), mean)
tau = apply(chains$tau, 2, mean)



# n = 10000
# x = rnorm(n)
# 
# # mean(x - mean(x) - var(x)/2)
# # 
# # mean(exp(x) - exp(mean(x) + var(x)/2))
# # 
# # mean(exp(x - mean(x) - var(x)/2))
# # mean(x-mean(x))
# 
# 
# theta = exp(x)
# mux = mean(theta)^2
# sx = var(x)
# 
# mu = log(mux/sqrt(mux+sx))
# s = log(1 + sx/mux)
# 
# mean(exp(x - mu - s/2))


# join data with subsetted data
omega_train = omega %>% 
  as_tibble() %>% 
  mutate(DOW = dat$lake_id) %>% 
  right_join(lselect) %>% 
  droplevels() %>% 
  select(-DOW) %>% 
  as.matrix() %>% 
  unname()

omega_test = omega %>% 
  as_tibble() %>% 
  mutate(DOW = dat$lake_id) %>% 
  right_join(tibble(DOW = dat_test$lake_id)) %>% 
  droplevels() %>% 
  select(-DOW) %>% 
  as.matrix() %>% 
  unname()
  

set.seed(1)
trainY = gen_Y(beta_0, beta, phi, Sigma_species, tau, omega_train, dat_train)
testY = gen_Y(beta_0, beta, phi, Sigma_species, tau, omega_test, dat_test)

dat_train$y_vec = trainY$y_vec
dat_train$Y = trainY$Y

testY$y_vec = testY$y_vec
testY$Y = testY$Y


rm(out, chains, center_temp, center_temp_sub)


# set up stan parameters --------------------------------------------------

out = stan(file = 'stan/spatial_jsdm_effort_scaling.stan',
           data = dat_train, iter = 2000, warmup = 1000, chains = 3, cores = 3, refresh = 10) # our model

# saveRDS(out, "data/stan_jsdm_simulation.rds")
# out = read_rds("data/stan_jsdm_simulation.rds")

chains = extract(out, permuted = T)
names(chains)
lapply(chains, dim)


b_names = colnames(dat$X)
phi_names = colnames(dat$Z)

apply(chains$beta_0, 2, mean)
b_hat = t(apply(chains$beta, c(2,3), mean))
phi_hat = t(apply(chains$phi, c(2,3), mean))
fnames = dat$fish_names %>% str_replace(., ' ', '_')

round(b_hat, 3) %>%
  as_tibble(.name_repair = ~b_names) %>%
  mutate(Species = dat$fish_names) %>%
  relocate(Species)

round(phi_hat, 3) %>%
  as_tibble(.name_repair = ~phi_names) %>%
  mutate(Species = dat$fish_names) %>%
  relocate(Species)



par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$beta_0[,i], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$beta[,4,i], type = 'l', main = i)
}

par(mfrow = c(3,4))
for(i in 1:12){
  plot(chains$phi[,i,2], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$omega[,40,i], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$z[,3,i], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$Sigma_species[,i,5], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$tau[,i], type = 'l', main = i)
}


plot(chains$lp__, type = 'l')





# check coverage ----------------------------------------------------------

b0upper = apply(chains$beta_0, c(2), quantile, probs = 0.975)
b0lower = apply(chains$beta_0, c(2), quantile, probs = 0.025)

bupper = apply(chains$beta, c(2,3), quantile, probs = 0.975)
blower = apply(chains$beta, c(2,3), quantile, probs = 0.025)

pupper = apply(chains$phi, c(2,3), quantile, probs = 0.975)
plower = apply(chains$phi, c(2,3), quantile, probs = 0.025)

oupper = apply(chains$omega, c(2,3), quantile, probs = 0.975)
olower = apply(chains$omega, c(2,3), quantile, probs = 0.025)

supper = apply(chains$Sigma_species, c(2,3), quantile, probs = 0.975)
slower = apply(chains$Sigma_species, c(2,3), quantile, probs = 0.025)

tupper = apply(chains$tau, c(2), quantile, probs = 0.975)
tlower = apply(chains$tau, c(2), quantile, probs = 0.025)




(beta_0 > b0lower) & (beta_0 < b0upper)
sign(b0lower) == sign(b0upper)
covBeta0 = sum((beta_0 > b0lower) & (beta_0 < b0upper))/length(beta_0)

(beta > blower) & (beta < bupper)
sign(blower) == sign(bupper)
covBeta = sum((beta > blower) & (beta < bupper))/length(beta)

(phi > plower) & (phi < pupper)
sign(plower) == sign(pupper)
covPhi = sum((phi > plower) & (phi < pupper))/length(phi)


(omega_train > olower) & (omega_train < oupper)
sign(olower) == sign(oupper)
covOmega = sum((omega_train > olower) & (omega_train < oupper))/length(omega_train)
sum(sign(olower) == sign(oupper))/length(omega_train)



(Sigma_species > slower) & (Sigma_species < supper)
sign(slower) == sign(supper)
# don't include diagonal and is symmetric
covSigma = sum(((Sigma_species > slower) & (Sigma_species < supper))[upper.tri((Sigma_species > slower) & (Sigma_species < supper))])/length((((Sigma_species > slower) & (Sigma_species < supper))[upper.tri((Sigma_species > slower) & (Sigma_species < supper))]))


(tau > tlower) & (tau < tupper)
sign(tlower) == sign(tupper)
covTau = sum((tau > tlower) & (tau < tupper))/length(tau)

# coverage percentage
sum(covBeta0*length(beta_0),
    covBeta*length(beta),
    covPhi*length(phi),
    covOmega*length(omega_train),
    covSigma*15,
    covTau*length(tau))/sum(length(beta_0),length(beta),length(phi),length(omega_train),15,length(tau))

sum(covBeta0*length(beta_0),
    covBeta*length(beta),
    covPhi*length(phi),
    covSigma*15,
    covTau*length(tau))/sum(length(beta),length(phi),15,length(tau))




# predictability of lambda ------------------------------------------------


# chains = extract(out, permuted = T)

beta_0_hat = apply(chains$beta_0, 2, mean)
beta_hat = apply(chains$beta, c(2,3), mean)
phi_hat = apply(chains$phi, c(2,3), mean)
omega_hat = apply(chains$omega, c(2,3), mean)
Sigma_species_hat = apply(chains$Sigma_species, c(2,3), mean)
tau_hat = apply(chains$tau, 2, mean)

omega_test_hat = omega_hat %>% 
  as_tibble() %>% 
  mutate(DOW = dat_train$lake_id) %>% 
  right_join(tibble(DOW = dat_test$lake_id)) %>% 
  droplevels() %>% 
  select(-DOW) %>% 
  as.matrix() %>% 
  unname()


computeLambda = function(beta_0, beta, phi, Sigma_species, tau, omega, dat){
  
  OMEGA = dat$each_lake %*% omega
  
  # relative abundance
  B0 = matrix(rep(beta_0, dat$N), dat$N, dat$K, byrow = T)
  gamma = B0 + dat$X %*% beta + OMEGA
  
  # catchability
  theta_star = dat$Z %*% phi
  theta_star_star = dat$Zstar %*% phi
  # theta = theta_star - mean(theta_star_star) - var(c(theta_star_star))/2
  
  mux = mean(exp(theta_star_star))^2
  sx = var(exp(c(theta_star_star)))
  mu = log(mux/sqrt(mux+sx))
  s = log(1 + sx/mux)
  
  theta = theta_star - mu - s/2
  
  
  # lambda = exp(log(dat$effort) + theta + gamma)
  # return(lambda)
  # return(exp(gamma))
  return(gamma)
  
}



theta_star = dat$Z %*% phi
theta_star_star = dat$Zstar %*% phi

mux = mean(exp(theta_star_star))^2
sx = var(exp(c(theta_star_star)))
mu = log(mux/sqrt(mux+sx))
s = log(1 + sx/mux)

theta = theta_star - mu - s/2




# computeLambda(beta_0_hat, beta_hat, phi_hat, Sigma_species_hat, tau_hat, omega_train, dat_train)
lambdaEst = computeLambda(beta_0_hat, beta_hat, phi_hat, Sigma_species_hat, tau_hat, omega_test_hat, dat_test)
lambdaTrue = computeLambda(beta_0, beta, phi, Sigma_species, tau, omega_test, dat_test)


lamHat = array(NA, dim = c(310, 6, 3000))

for(i in 1:3000){

  omega_test_hat = chains$omega[i,,] %>% 
    as_tibble() %>% 
    mutate(DOW = dat_train$lake_id) %>% 
    right_join(tibble(DOW = dat_test$lake_id)) %>% 
    droplevels() %>% 
    select(-DOW) %>% 
    as.matrix() %>% 
    unname()
  
  lamHat[,,i] = computeLambda(chains$beta_0[i,], chains$beta[i,,], chains$phi[i,,], chains$Sigma_species[i,,],
                chains$tau[i,], omega_test_hat, dat_test)
  
}


Est2008 = lamHat[yearind == 2008,,]
Est2018 = lamHat[yearind == 2018,,]
True2008 = lambdaTrue[yearind == 2008,]
True2018 = lambdaTrue[yearind == 2018,]

dSqrd2008 = array(NA, dim = dim(Est2008)[2:3])
dSqrd2018 = array(NA, dim = dim(Est2018)[2:3])
for(i in 1:3000){
  dSqrd2008[,i] = apply((Est2008[,,i] - True2008)^2, 2, function(x){sqrt(mean(x))})
  dSqrd2018[,i] = apply((Est2018[,,i] - True2018)^2, 2, function(x){sqrt(mean(x))})
}



apply(dSqrd2008, 1, mean)
apply(dSqrd2018, 1, mean)

apply(dSqrd2008, 1, median)
apply(dSqrd2018, 1, median)



tmp = as_tibble(t(dSqrd2008), .name_repair = ~fnames) %>% 
  mutate(id = 1:n(),
         year = 2008) %>% 
  pivot_longer(-c(id,year))

tmp2 = as_tibble(t(dSqrd2018), .name_repair = ~fnames) %>% 
  mutate(id = 1:n(),
         year = 2018) %>% 
  pivot_longer(-c(id,year))

tmp = rbind(tmp, tmp2) %>% 
  mutate(year = factor(year))

ggplot(tmp, aes(x = value)) +
  geom_histogram(aes(color = year)) +
  facet_wrap(~name, scales = "free")

# 
# dSqrd = array(NA, dim = c(310, 6, 3000))
# for(i in 1:3000){
#   dSqrd[,,i] = (lamHat[,,i] - lambdaTrue)^2
# }
# 
# apply(apply(dSqrd, c(2,3), function(X){sqrt(sum(X)/300)}), 1, mean)
# 
# 
# lambdaEst = apply(lamHat, c(1,2), mean)
# 
# 
# sqrt(sum((lambdaEst - lambdaTrue)^2)/length(lambdaEst))
# 
# 
# plot(c(lambdaEst), c(lambdaTrue))
# 
# plot(log(c(lambdaEst)), log(c(lambdaTrue)))
# lines(log(1:2000), log(1:2000))
# lines(-5:10,-5:10)


tibble(Estimate = lambdaEst, Target = lambdaTrue)


fnames = dat$fish_names
fnames = fnames %>%
  str_replace_all(., "_", " ")


yearind = fish_dat_sub %>% 
  mutate(year = year(SURVEYDATE)) %>% 
  filter(year == 2008 | year == 2018) %>% 
  droplevels() %>% 
  select(DOW, year, GN) %>% 
  group_by(DOW, year) %>% 
  distinct() %>% 
  ungroup() %>% 
  select(year) %>% 
  pull


preds = as_tibble(lambdaEst, .name_repair = ~fnames) %>% 
  mutate(ID = 1:n(),
         year = yearind) %>% 
  pivot_longer(-c(ID,year), names_to = "fish", values_to = "Predicted") %>% 
  left_join(as_tibble(lambdaTrue, .name_repair = ~fnames) %>% 
              mutate(ID = 1:n(),
                     year = yearind) %>% 
              pivot_longer(-c(ID,year), names_to = "fish", values_to = "Observed"),
            by = c("ID", "fish", 'year'))


ggplot(preds, aes(x = log(Predicted), y = log(Observed), color = factor(year))) +
  geom_point() +
  scale_color_manual(values = c('purple', 'forest green'), labels = c('9', '18'), name = "") +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fish, scales = 'free') +
  xlab("Log Predicted Relative Abundance") +
  ylab("Log Observed Relative Abundance")


ggplot(preds, aes(x = Predicted, y = Observed, color = factor(year))) +
  geom_point() +
  scale_color_manual(values = c('purple', 'forest green'), labels = c('9', '18'), name = "") +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fish, scales = 'free') +
  xlab("Predicted Relative Abundance") +
  ylab("Observed Relative Abundance")


# ggsave('results/OOSpredictions.png', width = 9, height = 7)
# ggsave('results/eps_figures/OOSpredictions.eps', dpi = 600, width = 9, height = 7, device=grDevices::cairo_ps)
# ggsave('results/OOSpredictionsLog.png', width = 9, height = 7)
# ggsave('results/eps_figures/OOSpredictionsLog.eps', dpi = 600, width = 9, height = 7, device=grDevices::cairo_ps)



preds %>% 
  group_by(fish, year) %>% 
  summarize(RMSE = sqrt(mean((Predicted - Observed)^2))) %>% 
  ungroup() %>% 
  pivot_wider(names_from = year, values_from = RMSE) %>% 
  as.data.frame()



# oos coverage ----------------------------------------------------------


checkCoverage = function(est, target, i, j, l, u){
  
  vals = quantile(est[i,j,], probs = c(l,u))
  if((target[i,j] > vals[1]) && (target[i,j] < vals[2])){
    return(1)
  }else{
    return(0)
  }
  
}


covMat = matrix(0, 310, 6)
for(i in 1:310){
  for(j in 1:6){
    covMat[i,j] = checkCoverage(lamHat, lambdaTrue, i, j, 0.025, 0.975)
  }
}

as_tibble(covMat, .name_repair = ~fnames) %>% 
  mutate(ID = 1:n(),
         year = yearind) %>% 
  pivot_longer(-c(ID,year), names_to = "fish", values_to = "cInd") %>% 
  group_by(year, fish) %>% 
  summarize(coverage = mean(cInd)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = year, values_from = coverage) %>% 
  as.data.frame()

t(rbind(fnames, apply(covMat, 2, mean)))
fnames
