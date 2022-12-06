

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
  slice_head(n = 200) %>% 
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
mean_covs = colnames(fish_dat)[c(7, 9, 23, 24, 14:15)]
temporal_covs = colnames(fish_dat)[c(23, 24)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(14:15)]
catch_covs = colnames(fish_dat)[c(25, 27:30)] # temp doy interaction
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
    select(all_of(mean_covs)) %>% 
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
  pars$Z = as.matrix(Z)
  pars$Zstar = rbind(TN, GN) %>% as.matrix()
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

dat = create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, center_temp)

dat_train = create_pars(fish_dat_sub %>% 
                          mutate(year = year(SURVEYDATE)) %>% 
                          filter(year != 2008 | year != 2018) %>% 
                          droplevels(), 
                        mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, 
                        center_temp_sub %>% 
                          mutate(year = year(SURVEYDATE)) %>% 
                          filter(year != 2008 | year != 2018) %>% 
                          droplevels())
dat_test = create_pars(fish_dat_sub %>% 
                         mutate(year = year(SURVEYDATE)) %>% 
                         filter(year == 2008 | year == 2018) %>% 
                         droplevels(), 
                       mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, 
                       center_temp_sub %>% 
                         mutate(year = year(SURVEYDATE)) %>% 
                         filter(year == 2008 | year == 2018) %>% 
                         droplevels())

str(dat)
str(dat_train)
str(dat_test)


# simulate data ----------------------------------------------------------

# set.seed(1)
# 
# beta_0 = rnorm(dat$K)
# beta = matrix(rnorm(dat$K*dat$p_beta, 0, 0.3), dat$p_beta, dat$K)
# Sigma_species = cov2cor(rinvwishart(12, diag(rep(1, dat$K))))
# tau = c(1.14, 1.2, 0.5, 1.23, 0.8, 1.1)
# phi = matrix(rnorm(dat$K*dat$p_phi, 0, 0.3), dat$p_phi, dat$K)
# 
# z = matrix(rnorm(dat$K*dat$n_lakes), dat$n_lakes, dat$K)
# 
# 
# # intercept
# B0 = matrix(rep(beta_0, dat$N), dat$N, dat$K, byrow = T)
# 
# # Random effects - mean zero by species
# omega_star = z %*% diag(tau) %*% Sigma_species %*% diag(tau)
# omega_mean = apply(omega_star, 2, mean) # calculate mean of omega
# omega = omega_star - matrix(rep(omega_mean, dat$n_lakes), dat$n_lakes, dat$K, byrow = T)
# OMEGA = dat$each_lake %*% omega
# 
# # relative abundance
# gamma = B0 + dat$X %*% beta + OMEGA
# range(exp(gamma))
# 
# # catchability
# theta_star = dat$Z %*% phi
# theta_star_star = dat$Zstar %*% phi
# # theta = theta_star - mean(theta_star_star) - variance(theta_star_star)/2
# theta = theta_star - mean(theta_star_star) - var(c(theta_star_star))/2
# # theta = theta_star - apply(theta_star_star, 2, mean) - apply(theta_star_star, 2, var)/2
# 
# 
# lambda = exp(log(dat$effort) + theta + gamma)
# y = rpois(dat$N*dat$K, c(lambda))
# Y = matrix(y, dat$N, dat$K)
# 
# dat$y_vec = y
# dat$Y = Y


gen_Y = function(beta_0, beta, phi, Sigma_species, tau, omega, dat){
  
  OMEGA = dat$each_lake %*% omega
  
  # relative abundance
  B0 = matrix(rep(beta_0, dat$N), dat$N, dat$K, byrow = T)
  gamma = B0 + dat$X %*% beta + OMEGA
  
  # catchability
  theta_star = dat$Z %*% phi
  theta_star_star = dat$Zstar %*% phi
  theta = theta_star - mean(theta_star_star) - var(c(theta_star_star))/2
  
  
  lambda = exp(log(dat$effort) + theta + gamma)
  y = rpois(dat$N*dat$K, c(lambda))
  Y = matrix(y, dat$N, dat$K)
  
  return(list(y_vec = y, Y = Y))
  
}


out = read_rds("/Users/jsnow/Desktop/stan_lake_random_effect.rds")
chains = extract(out, permuted = T)

beta0 = apply(chains$beta_0, 2, mean)
beta = apply(chains$beta, c(2,3), mean)
phi = apply(chains$phi, c(2,3), mean)
omega = apply(chains$omega, c(2,3), mean)
Sigma_species = apply(chains$Sigma_species, c(2,3), mean)
tau = apply(chains$tau, 2, mean)



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




# set up stan parameters --------------------------------------------------

out = stan(file = 'stan/spatial_jsdm_effort_scaling.stan', data = dat_train, iter = 200, warmup = 100, chains = 3, cores = 3, refresh = 10) # our model

# saveRDS(out, "data/stan_jsdm_simulation.rds")
# out = read_rds("data/stan_jsdm_simulation.rds")

chains = extract(out, permuted = T)
names(chains)
lapply(chains, dim)



b_names = colnames(dat$X)
phi_names = colnames(dat$Z)

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

par(mfrow = c(2,4))
for(i in 1:8){
  plot(chains$beta[,i,1], type = 'l', main = i)
}

par(mfrow = c(3,4))
for(i in 1:12){
  plot(chains$phi[,i,5], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$omega[,3,i], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$Sigma_species[,i,1], type = 'l', main = i)
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
sum(sign(omega_train) == sign(omega_train))/length(omega_train)


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



# predictability of lambda ------------------------------------------------


chains = extract(out, permuted = T)

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
  theta = theta_star - mean(theta_star_star) - var(c(theta_star_star))/2
  
  
  lambda = exp(log(dat$effort) + theta + gamma)
  return(lambda)
  
}

# computeLambda(beta_0_hat, beta_hat, phi_hat, Sigma_species_hat, tau_hat, omega_train, dat_train)
lambdaEst = computeLambda(beta_0_hat, beta_hat, phi_hat, Sigma_species_hat, tau_hat, omega_test_hat, dat_test)
lambdaTrue = computeLambda(beta_0, beta, phi, Sigma_species, tau, omega_test, dat_test)


sqrt(sum((lambdaEst - lambdaTrue)^2))
plot(c(lambdaEst), c(lambdaTrue), xlim = c(0, 1000))










