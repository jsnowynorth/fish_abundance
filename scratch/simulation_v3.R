

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
  slice_head(n = 100) %>% 
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


fnames = fish_dat_sub %>% 
  distinct(COMMON_NAME) %>% 
  pull()


fish_dat %>% 
  select(DOW, secchi, SURVEYDATE) %>% 
  mutate(year = year(SURVEYDATE)) %>% 
  distinct(DOW, year, .keep_all = T) %>% 
  right_join(lselect) %>% 
  droplevels()

temp_sampled = center_temp %>% 
  right_join(lselect) %>% 
  group_by(year, DOW) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  arrange(DOW)
  


fish_sim = fish_dat_sub %>% 
  distinct(DOW, .keep_all = T) %>% 
  droplevels() %>% 
  select(DOW, TOTAL_CATCH, MAX_DEPTH_FEET, TOTAL_CATCH, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>% 
  mutate("black crappie" = 1,
         "bluegill" = 1,
         "largemouth bass" = 1,
         "northern pike" = 1,
         "walleye" = 1,
         "yellow perch" = 1,
         "TN" = 1,
         "GN" = 0) %>% 
  right_join(temp_sampled, by = c("DOW")) %>% 
  left_join(fish_dat %>% select(DOW, secchi, year) %>% distinct(DOW, year, .keep_all = T)) %>% 
  pivot_longer(`black crappie`:`yellow perch`, names_to = "COMMON_NAME", values_to = "val") %>% 
  pivot_longer(c(TN, GN), names_to = "gear", values_to = "TN") %>% 
  mutate("GN" = 1-TN) %>% 
  select(-c(val, gear)) %>% 
  mutate(COMMON_NAME = as.factor(COMMON_NAME),
         EFFORT = 1) %>% 
  mutate(DOY = yday(SURVEYDATE)) %>% 
  mutate(DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>% 
  droplevels() %>% 
  group_by(DOW) %>% 
  mutate(secchi = ifelse(is.na(secchi), mean(secchi, na.rm = T), secchi)) %>% 
  ungroup() %>% 
  relocate("DOW", "secchi", "MAX_DEPTH_FEET", "TOTAL_CATCH", "LAKE_CENTER_UTM_EASTING", "LAKE_CENTER_UTM_NORTHING", "SURVEYDATE" , "COMMON_NAME", "TN", "GN", "EFFORT", "temp_0", "year", "DOY" , "DOY_sin_semi", "DOY_cos_semi", "DOY_sin_semi_temp", "DOY_cos_semi_temp")


# nreps = 20
# 
# fish_sim = fish_dat_sub %>% 
#   distinct(DOW, .keep_all = T) %>% 
#   droplevels() %>% 
#   select(DOW, secchi, MAX_DEPTH_FEET, TOTAL_CATCH, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>% 
#   mutate("black crappie" = 1,
#          "bluegill" = 1,
#          "largemouth bass" = 1,
#          "northern pike" = 1,
#          "walleye" = 1,
#          "yellow perch" = 1,
#          "TN" = 1,
#          "GN" = 0) %>% 
#   mutate(count = nreps) %>% 
#   uncount(count) %>% 
#   mutate(SURVEYDATE = rep(as_date(mdy("07-01-16"):(mdy("07-01-16") + days(nreps-1))), 100)) %>% 
#   pivot_longer(`black crappie`:`yellow perch`, names_to = "COMMON_NAME", values_to = "val") %>% 
#   pivot_longer(c(TN, GN), names_to = "gear", values_to = "TN") %>% 
#   mutate("GN" = 1-TN) %>% 
#   select(-c(val, gear)) %>% 
#   mutate(COMMON_NAME = as.factor(COMMON_NAME),
#          EFFORT = 1) %>% 
#   left_join(center_temp, by = c("DOW", "SURVEYDATE")) %>%
#   mutate(DOY = yday(SURVEYDATE)) %>% 
#   mutate(DOY_sin_semi = sin(DOY/121 * 2*pi),
#          DOY_cos_semi = cos(DOY/121 * 2*pi),
#          DOY_sin_semi_temp = DOY_sin_semi * temp_0,
#          DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>% 
#   droplevels()


mean_covs = colnames(fish_sim)[c(2, 3)]
temporal_covs = NULL
mean_covs_log = colnames(fish_sim)[c(3)]
mean_covs_logit = NULL
# catch_covs = colnames(fish_dat)[c(25, 27:30)]
catch_covs = colnames(fish_dat)[c(25, 27:30)]
# catch_covs = NULL
gear_types = colnames(fish_sim)[c(9, 10)]

# dat = fish_sim


create_pars <- function(dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, center_temp){
  
  dat = dat %>% 
    arrange(DOW, SURVEYDATE, COMMON_NAME, GN)
  
  K = nlevels(dat$COMMON_NAME)
  levs = levels(dat$COMMON_NAME)
  lake_index = (dat %>% filter(COMMON_NAME == levs[1]))$DOW
  lake_id = levels(dat$DOW)
  n_lakes = length(levels(dat$DOW))
  n_obs = (dat %>% filter(COMMON_NAME == levs[1]) %>% select(TOTAL_CATCH))$TOTAL_CATCH
  
  pars = list()
  
  
  pars$Y = dat %>% 
    arrange(DOW, COMMON_NAME, GN) %>% 
    select(DOW, SURVEYDATE, COMMON_NAME, GN, TOTAL_CATCH) %>% 
    pivot_wider(names_from = COMMON_NAME, values_from = TOTAL_CATCH) %>% 
    select(-c(DOW, SURVEYDATE, GN)) %>% 
    as.matrix()
  
  pars$effort = dat %>% 
    arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
    select(DOW, SURVEYDATE, COMMON_NAME, GN, EFFORT) %>% 
    pivot_wider(names_from = COMMON_NAME, values_from = EFFORT) %>% 
    select(-c(DOW, SURVEYDATE, GN)) %>% 
    as.matrix()
  
  X = dat %>% 
    arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
    filter(COMMON_NAME == levs[1]) %>% 
    select(all_of(mean_covs)) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.))
  # mutate(secchi = secchi - mean(secchi))
  
  Z = dat %>%
    arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>%
    filter(COMMON_NAME == levs[1]) %>%
    select(all_of(catch_covs), GN) %>%
    mutate(TN = 1) %>%
    relocate(TN) %>%
    mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN))
  # %>%
    # rename("temp_0_TN" = "temp_0",
           # "temp_0_GN" = "GN")
  
  # Z = dat %>% 
  #   arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  #   filter(COMMON_NAME == levs[1]) %>%
  #   select(all_of(catch_covs), GN) %>% 
  #   mutate(TN = ifelse(GN == 1, 0, 1)) %>% 
  #   relocate(TN) %>% 
  #   mutate_at(vars(all_of(catch_covs)), .funs = list(TN = ~.*TN)) %>% 
  #   mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN)) %>% 
  #   select(-all_of(catch_covs))
  
  temp_obs_lakes = center_temp %>% 
    filter(DOW %in% unique(dat$DOW)) %>% 
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
  
  pars$X = cbind(rep(1, dim(X)[1]), as.matrix(X))
  pars$Z = as.matrix(Z %>% select(-c(TN, GN)))
  pars$Zstar = rbind(TN, GN) %>% select(-c(TN, GN)) %>% as.matrix()
  
  # pars$Z = as.matrix(Z %>% select(-c(TN)))
  # pars$Zstar = rbind(TN, GN) %>% select(c(temp_0, temp_0_GN)) %>% as.matrix()
  
  # pars$Z = as.matrix(Z)
  # pars$Zstar = rbind(TN, GN) %>% select(c(TN, GN, temp_0, temp_0_GN)) %>% as.matrix()
  
  # pars$Zstar = as.matrix(Z %>% select(-c(TN)))
  # pars$Z = as.matrix(Z)
  # pars$Zstar = rbind(TN, GN) %>% select(c(TN, GN)) %>% as.matrix()
  pars$y_vec = c(pars$Y)
  
  pars$Nstar = dim(pars$Zstar)[1]
  
  # parameters
  pars$N = nrow(pars$Y)
  pars$p_beta = ncol(pars$X)
  pars$p_phi = ncol(pars$Z)
  pars$K = K
  pars$ngears = 2
  
  # indexing
  pars$lake_index = lake_index
  # pars$each_lake = table(lake_index)
  pars$lake_id = lake_id
  pars$n_lakes = n_lakes
  pars$fish_names = levels(dat$COMMON_NAME)
  
  
  pars$each_lake = dat %>% 
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


# dat = create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, center_temp)

dat_train = create_pars(fish_sim,
                        mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs,
                        center_temp_sub)



str(dat)
str(dat_train)
str(dat_test)


# simulate data ----------------------------------------------------------

gen_Y = function(beta, theta, phi, omega, dat){
  
  
  # relative abundance
  OMEGA = dat$each_lake %*% omega
  gamma = dat$X %*% beta + OMEGA
  
  # catchability temporal
  phiz_star = dat$Z %*% phi
  phiz_star_star = dat$Zstar %*% phi
  mux = mean(exp(phiz_star_star))^2
  sx = var(exp(c(phiz_star_star)))
  mu = log(mux/sqrt(mux+sx))
  s = log(1 + sx/mux)

  phiz = phiz_star - mu - s/2
  
  # phiz = dat$Z %*% phi
  
  
  lambda = exp(log(dat$effort * theta) + phiz + gamma)
  # lambda = exp(log(dat$effort * theta) + gamma)
  # lambda = exp(log(dat$effort) + gamma)
  y = rpois(dat$N*dat$K, c(lambda))
  Y = matrix(y, dat$N, dat$K)
  summary(Y)
  
  return(list(y_vec = y, Y = Y))
  
}


out = read_rds("/Users/jsnow/Desktop/stan_lake_random_effect.rds")
chains = extract(out, permuted = T)


set.seed(1)
# beta
# beta_0 = apply(chains$beta_0, 2, mean)
# # beta_0[1] = 0
# beta = apply(chains$beta, c(2,3), mean)[c(1),]
# beta = rbind(beta_0, beta)

# beta_0 = c(2.5, 3.5, 0.5, 1.5, 0.9, 0.7)
# beta1 = c(-0.2, 0.2, 0.3, -0.1, -0.1, 0.4)
# beta2 = c(0.2, -0.2, -0.1, 0.1, 0.4, 0.2)
# beta = rbind(beta_0, beta1, beta2)

beta_0 = c(2.5, 3.5, 0.5, 1.5, 0.9, 0.7)
beta1 = c(-0.4, 0.3, 0.4, -0.3, -0.2, 0.5)
beta2 = c(0.3, -0.5, -0.2, 0.8, 1, 0.3)
beta = rbind(beta_0, beta1, beta2)


# theta
theta = matrix(1, 2, 6, byrow = T)
dfact = runif(6, min = -0.9, max = 0.9)
theta[1,] = theta[1,] - dfact
theta[2,] = theta[2,] + dfact

nreps = dim(dat_train$effort)[1]

# theta = matrix(rep(theta, nreps/2), nreps, dat_train$K, byrow = T)
theta = matrix(rep(t(theta), nreps/2), nreps, dat_train$K, byrow = T)

theta[1:2,]

# phi_int = solve(t(dat_train$Z[1:2,1:2]) %*% dat_train$Z[1:2,1:2]) %*% t(dat_train$Z[1:2,1:2]) %*% exp(theta[1:2,])


# phi_int_TN = c(4, 2, 3.5, 2, 2.5, 3)
# phi_int_GN = c(0.5, 1, 1.2, -0.2, 0.4, 0.3)
# phi = solve(t(dat_train$Z[1:2,]) %*% dat_train$Z[1:2,]) %*% t(dat_train$Z[1:2,]) %*% log(theta[1:2,])

# x = rnorm(dat_train$p_phi * dat_train$K)
# phi = matrix(x - mean(x), nrow = dat_train$p_phi)
# sum(phi)


phi = matrix(runif(dat_train$K * dat_train$p_phi, min = -1, max = 1), dat_train$p_phi, dat_train$K)

# phi = matrix(1, 2, 6, byrow = T)
# dfact = runif(6, min = -0.9, max = 0.9)
# phi[1,] = phi[1,] - dfact
# phi[2,] = phi[2,] + dfact

# phi = rbind(phi_int_TN, phi_int_GN, phi)


# random effect
# omega = apply(chains$omega, c(2,3), mean)
Sigma_species = apply(chains$Sigma_species, c(2,3), mean)
tau = apply(chains$tau, 2, mean)

Sig_tmp = diag(tau) %*% Sigma_species %*% diag(tau)


U = t(chol(Sig_tmp))
z = t(U %*% t(matrix(rnorm(dat_train$n_lakes*dat_train$K), dat_train$n_lakes, dat_train$K)))
omega = z - matrix(rep(t(solve(Sig_tmp) %*% rep(1, dat_train$K) %*% (1/(sum(solve(Sig_tmp))*dat_train$n_lakes)) %*% sum(z)), dat_train$n_lakes), ncol = dat_train$K, byrow = T)

# z = t(U %*% t(matrix(rnorm(10000 * 6), 10000, 6)))
# z = matrix(rnorm(10000 * 6), 10000, 6) %*% U
# omega = z - matrix(rep(t(solve(Sig_tmp) %*% rep(1, 6) %*% (1/(sum(solve(Sig_tmp))*10000)) %*% sum(z)), 10000), ncol = 6, byrow = T)
# var(omega)
# 
# 
# 
# (diag(tau) %*% t(chol(Sigma_species))) %*% t(diag(tau) %*% t(chol(Sigma_species)))


omega = matrix(0, dat_train$n_lakes, dat_train$K)


# set.seed(1)
trainY = gen_Y(beta, theta, phi, omega, dat_train)
# trainY = gen_Y(beta, phi, omega, dat_train)
summary(trainY$Y)
apply(trainY$Y, 2, quantile, probs = seq(0,1,0.01))

dat_train$y_vec = trainY$y_vec
dat_train$Y = trainY$Y


rm(out, chains, center_temp, center_temp_sub)



# set up stan parameters --------------------------------------------------

# out = stan(file = 'scratch/jdsm_effort_scaling.stan',
#            data = dat_train, iter = 500, warmup = 250, chains = 3, cores = 3, refresh = 10, init = 0) # our model

# dat_train$X = cbind(dat_train$X, dat_train$Z)
# # dat_train$X = cbind(dat_train$X, dat_train$Z[,c(1,6)])
# dat_train$p_beta = dim(dat_train$X)[2]
# str(dat_train)

out = stan(file = 'scratch/jdsm_simplex.stan',
           data = dat_train, iter = 500, warmup = 250, chains = 3, cores = 3, refresh = 10) # our model

# saveRDS(out, "C:/Users/jsnow/Desktop/stan_jsdm_simulation.rds")
# out = read_rds("C:/Users/jsnow/Desktop/stan_jsdm_simulation.rds")

# out_naive = stan(file = 'scratch/jdsm_simplex_naive.stan',
#            data = dat_train, iter = 1000, warmup = 500, chains = 3, cores = 3, refresh = 10) # naive model

# saveRDS(out_naive, "C:/Users/jsnow/Desktop/stan_jsdm_simulation_naive.rds")
out_naive = read_rds("C:/Users/jsnow/Desktop/stan_jsdm_simulation_naive.rds")

chains = extract(out, permuted = T)
chains_naive = extract(out_naive, permuted = T)

names(chains)
lapply(chains, dim)

s <- summary(out, probs = c(0.05, 0.95))
s$summary


b_names = colnames(dat_train$X)
b_names[1] = "Int"
phi_names = colnames(dat_train$Z)


b_hat = t(apply(chains$beta, c(2,3), mean))
phi_hat = t(apply(chains$phi, c(2,3), mean))
tau_hat = t(apply(chains$tau, c(2), mean))
sig_corr_hat = t(apply(chains$Sigma_species, c(2,3), mean))
sigma_hat = t(apply(chains$Sig, c(2,3), mean))
fnames = dat_train$fish_names %>% str_replace(., ' ', '_')


diag(c(tau_hat)) %*% sig_corr_hat %*% t(sig_corr_hat) %*% diag(c(tau_hat))

apply(chains$theta, c(2,3), mean) * 12
theta[1:2,]


theta[1:2,] - apply(chains$theta, c(2,3), mean) * 12

phi - t(phi_hat)

beta - t(b_hat)



round(b_hat, 3) %>%
  as_tibble(.name_repair = ~b_names) %>%
  mutate(Species = dat_train$fish_names) %>%
  relocate(Species)

round(phi_hat, 3) %>%
  as_tibble(.name_repair = ~phi_names) %>%
  mutate(Species = dat_train$fish_names) %>%
  relocate(Species)


j=3
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains_naive$beta[,j,i], type = 'l', main = i)
  abline(h = beta[j, i], col = "blue")
}


j = 1
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$beta[,j,i], type = 'l', main = i)
  abline(h = beta[j, i], col = "blue")
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(12*chains$theta[,j,i], type = 'l', main = i)
  abline(h = theta[j, i], col = "blue")
}


j = 1
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$phi[,j,i], type = 'l', main = i)
  abline(h = phi[j, i], col = "blue")
}





par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$omega[,4,i], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$z[,3,i], type = 'l', main = i)
}


Sig = array(NA, dim = dim(chains$Sigma_species))
for (i in 1:dim(Sig)[1]){
  Sig[i,,] = diag(chains$tau[i,]) %*% chains$Sigma_species[i,,] %*% t(chains$Sigma_species[i,,]) %*% diag(chains$tau[i,])
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(Sig[,i,1], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$tau[,i], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$Sig[,i,5], type = 'l', main = i)
}

plot(chains$lp__, type = 'l')

cov2cor(apply(chains$Sig, c(2,3), mean))




# check coverage ----------------------------------------------------------

bupper_naive = apply(chains_naive$beta, c(2,3), quantile, probs = 0.975)
blower_naive = apply(chains_naive$beta, c(2,3), quantile, probs = 0.025)

bupper = apply(chains$beta, c(2,3), quantile, probs = 0.975)
blower = apply(chains$beta, c(2,3), quantile, probs = 0.025)

pupper = apply(chains$phi, c(2,3), quantile, probs = 0.975)
plower = apply(chains$phi, c(2,3), quantile, probs = 0.025)

tupper = 12*apply(chains$theta, c(2,3), quantile, probs = 0.975)
tlower = 12*apply(chains$theta, c(2,3), quantile, probs = 0.025)

oupper = 12*apply(chains$omega, c(2,3), quantile, probs = 0.975)
olower = 12*apply(chains$omega, c(2,3), quantile, probs = 0.025)


(beta > blower_naive) & (beta < bupper_naive)
sign(blower_naive) == sign(bupper_naive)
covBetaNaive = sum((beta > blower_naive) & (beta < bupper_naive))/length(beta)


(beta > blower) & (beta < bupper)
sign(blower) == sign(bupper)
covBeta = sum((beta > blower) & (beta < bupper))/length(beta)

(phi > plower) & (phi < pupper)
sign(plower) == sign(pupper)
covPhi = sum((phi > plower) & (phi < pupper))/length(phi)


(theta[1:2,] > tlower) & (theta[1:2,] < tupper)
sign(tlower) == sign(tupper)
covTheta = sum((theta[1:2,] > tlower) & (theta[1:2,] < tupper))/length(theta[1:2,])
sum(sign(tlower) == sign(tupper))/length(theta[1:2,])

(omega > olower) & (omega < oupper)
sign(olower) == sign(oupper)
covOmega = sum((omega > olower) & (omega < oupper))/length(omega)




# coverage percentage
sum(covBeta*length(beta),
    covPhi*length(phi),
    covTheta*length(theta[1:2,]),
    covOmega*length(omega))/sum(length(beta),length(phi),length(theta[1:2,]),length(omega))


sum(covBeta*length(beta),
    covPhi*length(phi),
    covTheta*length(theta[1:2,]))/sum(length(beta),length(phi),length(theta[1:2,]))



# rmse in sample ----------------------------------------------------------


computeGamma = function(beta, dat){
  
  # relative abundance
  gamma = dat$X %*% beta
  
  return(exp(gamma))
  
}


gammaTrue = computeGamma(beta, dat_train)

nlks = dim(gammaTrue)[1]
nsmps = dim(chains$beta)[1]

gamHat = array(NA, dim = c(nlks, 6, nsmps))
gamHat_naive = array(NA, dim = c(nlks, 6, nsmps))

for(i in 1:nsmps){
  
  gamHat[,,i] = computeGamma(chains$beta[i,,], dat_train)
  gamHat_naive[,,i] = computeGamma(chains_naive$beta[i,,], dat_train)
  
}


dSqrd = array(NA, dim = dim(gamHat)[2:3])
dSqrd_naive = array(NA, dim = dim(gamHat_naive)[2:3])

for(i in 1:nsmps){
  dSqrd[,i] = apply((gamHat[,,i] - gammaTrue)^2, 2, function(x){sqrt(mean(x))})
  dSqrd_naive[,i] = apply((gamHat_naive[,,i] - gammaTrue)^2, 2, function(x){sqrt(mean(x))})
}



apply(dSqrd, 1, mean)
apply(dSqrd_naive, 1, mean)


fnames = dat$fish_names
fnames = fnames %>%
  str_replace_all(., "_", " ")




preds = as_tibble(apply(gamHat, c(1,2), median), .name_repair = ~fnames) %>% 
  mutate(ID = 1:n()) %>% 
  pivot_longer(-c(ID), names_to = "fish", values_to = "Predicted") %>% 
  left_join(as_tibble(gammaTrue, .name_repair = ~fnames) %>% 
              mutate(ID = 1:n()) %>% 
              pivot_longer(-c(ID), names_to = "fish", values_to = "Observed"),
            by = c("ID", "fish"))

preds_naive = as_tibble(apply(gamHat_naive, c(1,2), median), .name_repair = ~fnames) %>% 
  mutate(ID = 1:n()) %>% 
  pivot_longer(-c(ID), names_to = "fish", values_to = "Predicted") %>% 
  left_join(as_tibble(gammaTrue, .name_repair = ~fnames) %>% 
              mutate(ID = 1:n()) %>% 
              pivot_longer(-c(ID), names_to = "fish", values_to = "Observed"),
            by = c("ID", "fish"))


p1 = ggplot(preds, aes(x = Predicted, y = Observed)) +
  geom_point() +
  scale_color_manual(values = c('purple', 'forest green'), labels = c('9', '18'), name = "") +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fish, scales = 'free') +
  xlab("Log Predicted Relative Abundance") +
  ylab("Log Observed Relative Abundance")


p2 = ggplot(preds_naive, aes(x = Predicted, y = Observed)) +
  geom_point() +
  scale_color_manual(values = c('purple', 'forest green'), labels = c('9', '18'), name = "") +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fish, scales = 'free') +
  xlab("Log Predicted Relative Abundance") +
  ylab("Log Observed Relative Abundance")


cowplot::plot_grid(p1, p2)

preds %>% 
  group_by(fish) %>% 
  summarize(RMSE = sqrt(mean((Predicted - Observed)^2))) %>% 
  ungroup() %>% 
  as.data.frame()

preds_naive %>% 
  group_by(fish) %>% 
  summarize(RMSE = sqrt(mean((Predicted - Observed)^2))) %>% 
  ungroup() %>% 
  as.data.frame()



# prediction data set -----------------------------------------------------

pselect = fish_dat %>% 
  mutate(cntr = 1) %>% 
  group_by(DOW) %>% 
  summarise(cnt = sum(cntr)) %>% 
  ungroup() %>% 
  arrange(desc(cnt)) %>% 
  slice(101:200) %>% 
  select(DOW) %>% 
  droplevels()

fish_dat_pred = fish_dat %>% 
  right_join(pselect) %>% 
  droplevels()

center_temp_pred = center_temp %>% 
  right_join(pselect) %>% 
  droplevels()


# temp_sampled_pred = center_temp %>% 
#   right_join(pselect) %>% 
#   filter(year == 2018) %>% 
#   group_by(DOW) %>% 
#   sample_n(1) %>% 
#   ungroup() %>% 
#   arrange(DOW)

temp_sampled_pred = center_temp %>% 
  right_join(pselect) %>% 
  group_by(DOW, year) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  arrange(DOW)



fish_pred = fish_dat_pred %>% 
  distinct(DOW, .keep_all = T) %>% 
  droplevels() %>% 
  select(DOW, TOTAL_CATCH, MAX_DEPTH_FEET, TOTAL_CATCH, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>% 
  mutate("black crappie" = 1,
         "bluegill" = 1,
         "largemouth bass" = 1,
         "northern pike" = 1,
         "walleye" = 1,
         "yellow perch" = 1,
         "TN" = 1,
         "GN" = 0) %>% 
  right_join(temp_sampled_pred, by = c("DOW")) %>% 
  left_join(fish_dat %>% select(DOW, secchi, year) %>% distinct(DOW, year, .keep_all = T)) %>% 
  pivot_longer(`black crappie`:`yellow perch`, names_to = "COMMON_NAME", values_to = "val") %>% 
  pivot_longer(c(TN, GN), names_to = "gear", values_to = "TN") %>% 
  mutate("GN" = 1-TN) %>% 
  select(-c(val, gear)) %>% 
  mutate(COMMON_NAME = as.factor(COMMON_NAME),
         EFFORT = 1) %>% 
  mutate(DOY = yday(SURVEYDATE)) %>% 
  mutate(DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>% 
  droplevels() %>% 
  group_by(DOW) %>% 
  filter(!is.na(secchi)) %>%  
  ungroup() %>% 
  relocate("DOW", "secchi", "MAX_DEPTH_FEET", "TOTAL_CATCH", "LAKE_CENTER_UTM_EASTING", "LAKE_CENTER_UTM_NORTHING", "SURVEYDATE" , "COMMON_NAME", "TN", "GN", "EFFORT", "temp_0", "year", "DOY" , "DOY_sin_semi", "DOY_cos_semi", "DOY_sin_semi_temp", "DOY_cos_semi_temp")








# pselect = fish_dat %>% 
#   mutate(cntr = 1) %>% 
#   group_by(DOW) %>% 
#   summarise(cnt = sum(cntr)) %>% 
#   ungroup() %>% 
#   arrange(desc(cnt)) %>% 
#   slice(101:150) %>% 
#   select(DOW) %>% 
#   droplevels()
# 
# fish_dat_pred = fish_dat %>% 
#   right_join(pselect) %>% 
#   droplevels()
# 
# center_temp_pred = center_temp %>% 
#   right_join(pselect) %>% 
#   droplevels()
# 
# 
# nreps = 30
# 
# fish_pred = fish_dat_pred %>% 
#   distinct(DOW, .keep_all = T) %>% 
#   droplevels() %>% 
#   select(DOW, secchi, MAX_DEPTH_FEET, TOTAL_CATCH, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>% 
#   mutate("black crappie" = 1,
#          "bluegill" = 1,
#          "largemouth bass" = 1,
#          "northern pike" = 1,
#          "walleye" = 1,
#          "yellow perch" = 1,
#          "TN" = 1,
#          "GN" = 0) %>% 
#   mutate(count = nreps) %>% 
#   uncount(count) %>% 
#   mutate(SURVEYDATE = rep(as_date(mdy("07-01-16"):(mdy("07-01-16") + days(nreps-1))), 50)) %>% 
#   pivot_longer(`black crappie`:`yellow perch`, names_to = "COMMON_NAME", values_to = "val") %>% 
#   pivot_longer(c(TN, GN), names_to = "gear", values_to = "TN") %>% 
#   mutate("GN" = 1-TN) %>% 
#   select(-c(val, gear)) %>% 
#   mutate(COMMON_NAME = as.factor(COMMON_NAME),
#          EFFORT = 1) %>% 
#   left_join(center_temp, by = c("DOW", "SURVEYDATE")) %>%
#   mutate(DOY = yday(SURVEYDATE)) %>% 
#   mutate(DOY_sin_semi = sin(DOY/121 * 2*pi),
#          DOY_cos_semi = cos(DOY/121 * 2*pi),
#          DOY_sin_semi_temp = DOY_sin_semi * temp_0,
#          DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>% 
#   droplevels()




dat_test = create_pars(fish_pred,
                        mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs,
                       center_temp_pred)

# predictability of lambda ------------------------------------------------


computeGamma = function(beta, dat){
  
  # relative abundance
  gamma = dat$X %*% beta

  return(exp(gamma))
  
}


gammaTrue = computeGamma(beta, dat_test)

nlks = dim(gammaTrue)[1]
nsmps = dim(chains$beta)[1]

gamHat = array(NA, dim = c(nlks, 6, nsmps))
gamHat_naive = array(NA, dim = c(nlks, 6, nsmps))

for(i in 1:nsmps){
  
  gamHat[,,i] = computeGamma(chains$beta[i,,], dat_test)
  gamHat_naive[,,i] = computeGamma(chains_naive$beta[i,,], dat_test)
  
}


dSqrd = array(NA, dim = dim(gamHat)[2:3])
dSqrd_naive = array(NA, dim = dim(gamHat_naive)[2:3])

for(i in 1:nsmps){
  dSqrd[,i] = apply((gamHat[,,i] - gammaTrue)^2, 2, function(x){sqrt(mean(x))})
  dSqrd_naive[,i] = apply((gamHat_naive[,,i] - gammaTrue)^2, 2, function(x){sqrt(mean(x))})
}



apply(dSqrd, 1, mean)
apply(dSqrd_naive, 1, mean)


fnames = dat$fish_names
fnames = fnames %>%
  str_replace_all(., "_", " ")


preds = as_tibble(apply(gamHat, c(1,2), mean), .name_repair = ~fnames) %>% 
  mutate(ID = 1:n()) %>% 
  pivot_longer(-c(ID), names_to = "fish", values_to = "Predicted") %>% 
  left_join(as_tibble(gammaTrue, .name_repair = ~fnames) %>% 
              mutate(ID = 1:n()) %>% 
              pivot_longer(-c(ID), names_to = "fish", values_to = "Observed"),
            by = c("ID", "fish"))

preds_naive = as_tibble(apply(gamHat_naive, c(1,2), mean), .name_repair = ~fnames) %>% 
  mutate(ID = 1:n()) %>% 
  pivot_longer(-c(ID), names_to = "fish", values_to = "Predicted") %>% 
  left_join(as_tibble(gammaTrue, .name_repair = ~fnames) %>% 
              mutate(ID = 1:n()) %>% 
              pivot_longer(-c(ID), names_to = "fish", values_to = "Observed"),
            by = c("ID", "fish"))


ggplot(preds, aes(x = Predicted, y = Observed)) +
  geom_point() +
  scale_color_manual(values = c('purple', 'forest green'), labels = c('9', '18'), name = "") +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fish, scales = 'free') +
  xlab("Log Predicted Relative Abundance") +
  ylab("Log Observed Relative Abundance")


ggplot(preds_naive, aes(x = Predicted, y = Observed)) +
  geom_point() +
  scale_color_manual(values = c('purple', 'forest green'), labels = c('9', '18'), name = "") +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fish, scales = 'free') +
  xlab("Log Predicted Relative Abundance") +
  ylab("Log Observed Relative Abundance")


preds %>% 
  group_by(fish) %>% 
  summarize(RMSE = sqrt(mean((Predicted - Observed)^2))) %>% 
  ungroup() %>% 
  as.data.frame()

preds_naive %>% 
  group_by(fish) %>% 
  summarize(RMSE = sqrt(mean((Predicted - Observed)^2))) %>% 
  ungroup() %>% 
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


covMat = matrix(0, nlks, 6)
covMat_naive = matrix(0, nlks, 6)
for(i in 1:nlks){
  for(j in 1:6){
    covMat[i,j] = checkCoverage(gamHat, gammaTrue, i, j, 0.025, 0.975)
    covMat_naive[i,j] = checkCoverage(gamHat_naive, gammaTrue, i, j, 0.025, 0.975)
  }
}

as_tibble(covMat, .name_repair = ~fnames) %>% 
  mutate(ID = 1:n()) %>% 
  pivot_longer(-c(ID), names_to = "fish", values_to = "cInd") %>% 
  group_by(fish) %>% 
  summarize(coverage = mean(cInd)) %>% 
  ungroup() %>% 
  as.data.frame()

as_tibble(covMat_naive, .name_repair = ~fnames) %>% 
  mutate(ID = 1:n()) %>% 
  pivot_longer(-c(ID), names_to = "fish", values_to = "cInd") %>% 
  group_by(fish) %>% 
  summarize(coverage = mean(cInd)) %>% 
  ungroup() %>% 
  as.data.frame()

