

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




mean_covs = colnames(fish_sim)[c(2, 3)]
temporal_covs = NULL
mean_covs_log = colnames(fish_sim)[c(3)]
mean_covs_logit = NULL
catch_covs = colnames(fish_dat)[c(25, 27:30)]
gear_types = colnames(fish_sim)[c(9, 10)]



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
  # pars$Z = as.matrix(Z %>% select(-c(TN, GN)))
  # pars$Zstar = rbind(TN, GN) %>% select(-c(TN, GN)) %>% as.matrix()
  
  pars$Z = as.matrix(Z)
  pars$Zstar = rbind(TN, GN) %>% as.matrix()
  

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


dat_train = create_pars(fish_sim %>% 
                          filter(year != 2019),
                        mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs,
                        center_temp_sub)



str(dat)
str(dat_train)
str(dat_test)


# simulate data ----------------------------------------------------------

gen_Y = function(beta, theta, phi, omega, dat){
  
  OMEGA = dat$each_lake %*% omega
  
  # relative abundance
  gamma = dat$X %*% beta + OMEGA
  
  # catchability temporal
  # phiz_star = dat$Z %*% phi
  # phiz_star_star = dat$Zstar %*% phi
  # mux = mean(exp(phiz_star_star))^2
  # sx = var(exp(c(phiz_star_star)))
  # mu = log(mux/sqrt(mux+sx))
  # s = log(1 + sx/mux)
  
  phiz_star = dat$Z %*% phi
  phiz_star_star = dat$Zstar %*% phi
  mux = apply(exp(phiz_star_star), 2, function(x){mean(x)^2})
  sx = apply(exp(phiz_star_star), 2, var)
  mu = log(mux/sqrt(mux+sx))
  s = log(1 + sx/mux)
  
  phiz = t(t(phiz_star) - mu - s/2)
  
  # phiz = dat$Z %*% phi

  lambda = exp(log(dat$effort) + phiz + gamma)
  # lambda = exp(log(dat$effort * theta) + phiz + gamma)
  # lambda = exp(log(dat$effort * theta) + gamma)
  y = rpois(dat$N*dat$K, c(lambda))
  Y = matrix(y, dat$N, dat$K)
  summary(Y)
  
  return(list(y_vec = y, Y = Y))
  
}


out = read_rds("/Users/jsnow/Desktop/stan_lake_random_effect.rds")
chains = extract(out, permuted = T)


set.seed(1)

beta_0 = c(2.5, 3.5, 0.5, 1.5, 0.9, 0.7)
beta1 = c(-0.4, 0.3, 0.4, -0.3, -0.2, 0.5)
beta2 = c(0.3, -0.5, -0.2, 0.8, 1, 0.3)
beta = rbind(beta_0, beta1, beta2)


# theta
theta = matrix(1, 2, 6, byrow = T)
dfact = runif(6, min = -0.9, max = 0.9)
theta[1,] = theta[1,] - dfact
theta[2,] = theta[2,] + dfact


# nreps = dim(dat_train$effort)[1]
# theta = matrix(rep(t(theta), nreps/2), nreps, dat_train$K, byrow = T)
# theta[1:2,]



phi = matrix(runif(dat_train$K * dat_train$p_phi, min = -1, max = 1), dat_train$p_phi, dat_train$K)
phi[1,] = theta[1,]
phi[7,] = theta[2,]

# random effect

# Sigma_species = apply(chains$Sigma_species, c(2,3), mean)
# tau = apply(chains$tau, 2, mean)
# Sig_tmp = diag(tau) %*% Sigma_species %*% t(diag(tau) %*% Sigma_species)

S = rinvwishart(20, diag(rep(1, dat_train$K)))
C = diag(1 / sqrt(diag(S))) %*% S %*% diag(1 / sqrt(diag(S)))
tau = c(0.4, 0.7, 0.6, 1.2, 0.8, 0.7)

Cchol = t(chol(C))

# cov2cor(diag(tau) %*% Cchol %*% t(diag(tau) %*% Cchol))
Sigma = diag(tau) %*% Cchol %*% t(diag(tau) %*% Cchol)


U = t(chol(Sigma))
z = t(U %*% t(matrix(rnorm(dat_train$n_lakes*dat_train$K), dat_train$n_lakes, dat_train$K)))
omega = z - matrix(rep(t(solve(Sigma) %*% rep(1, dat_train$K) %*% (1/(sum(solve(Sigma))*dat_train$n_lakes)) %*% sum(z)), dat_train$n_lakes), ncol = dat_train$K, byrow = T)

# z = t(U %*% t(matrix(rnorm(10000 * 6), 10000, 6)))
# omega = z - matrix(rep(t(solve(Sig_tmp) %*% rep(1, 6) %*% (1/(sum(solve(Sig_tmp))*10000)) %*% sum(z)), 10000), ncol = 6, byrow = T)
# var(omega)
# 
# 
# 
# (diag(tau) %*% t(chol(Sigma_species))) %*% t(diag(tau) %*% t(chol(Sigma_species)))


# omega = matrix(0, dat_train$n_lakes, dat_train$K)


# set.seed(1)
trainY = gen_Y(beta, theta, phi, omega, dat_train)
# trainY = gen_Y(beta, phi, omega, dat_train)
summary(trainY$Y)
apply(trainY$Y, 2, quantile, probs = seq(0,1,0.01))

dat_train$y_vec = trainY$y_vec
dat_train$Y = trainY$Y


rm(out, chains, center_temp, center_temp_sub)



# set up stan parameters --------------------------------------------------


#### Proposed model
# out = stan(file = 'scratch/jdsm_original_full_model.stan',
#            data = dat_train, iter = 1500, warmup = 1000, chains = 3, cores = 3, refresh = 10) # our model

# saveRDS(out, "C:/Users/jsnow/Desktop/original_model_sim.rds")
out = read_rds("C:/Users/jsnow/Desktop/original_model_sim.rds")


#### Naive model
# out_naive = stan(file = 'scratch/jdsm_naive_model.stan',
           # data = dat_train, iter = 1500, warmup = 1000, chains = 3, cores = 3, refresh = 10) # naive model

# saveRDS(out_naive, "C:/Users/jsnow/Desktop/naive_model_sim.rds")
out_naive = read_rds("C:/Users/jsnow/Desktop/naive_model_sim.rds")


# out = stan(file = 'scratch/jdsm_simples_v4.stan',
#            data = dat_train, iter = 250, warmup = 150, chains = 3, cores = 3, refresh = 10) # our model

# saveRDS(out, "C:/Users/jsnow/Desktop/stan_jsdm_simulation_v4.rds")
# out = read_rds("C:/Users/jsnow/Desktop/stan_jsdm_simulation_v4.rds")


# out_naive = read_rds("C:/Users/jsnow/Desktop/stan_jsdm_simulation_naive.rds")

chains = extract(out, permuted = T)
chains_naive = extract(out_naive, permuted = T)

names(chains)
lapply(chains, dim)



# check convergence
library(coda)

mcmcList = As.mcmc.list(out, pars = list("beta", "phi", "tau"))
gelman.diag(mcmcList)

mcmcList = As.mcmc.list(out, pars = list("tau"))
gelman.diag(mcmcList)

mcmcList = As.mcmc.list(out, pars = list("omega"))
gelman.diag(mcmcList)





b_names = colnames(dat_train$X)
b_names[1] = "Int"
phi_names = colnames(dat_train$Z)


b_hat = t(apply(chains$beta, c(2,3), mean))
# t_hat = t(apply(chains$theta, c(2,3), mean))
p_hat = t(apply(chains$phi, c(2,3), mean))
o_hat = t(apply(chains$omega, c(2,3), mean))
fnames = dat_train$fish_names %>% str_replace(., ' ', '_')

b_hat
t_hat
p_hat

t(b_hat) - beta
t(t_hat) - theta[1:2,]
t(p_hat) - phi


apply(chains$theta_simplex, c(2,3), mean)

t(apply(chains$theta, c(2,3), mean))*2

apply(chains$Sig, c(2,3), mean)

t(t_hat * 2)
theta[1:2,]

round(b_hat, 3) %>%
  as_tibble(.name_repair = ~b_names) %>%
  mutate(Species = dat_train$fish_names) %>%
  relocate(Species)



for(j in 1:3){
png(paste0("C:/Users/jsnow/Desktop/originalSimulationChains/beta_", j, ".png"))
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$beta[,j,i], type = 'l', main = i)
  abline(h = beta[j,i], col = "blue")
}
dev.off()
}


J = c(1, 2, 6, 7, 12)
for(j in J){
png(paste0("C:/Users/jsnow/Desktop/originalSimulationChains/phi_", j, ".png"))
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$phi[,j,i], type = 'l', main = i)
  abline(h = phi[j, i], col = "blue")
}
dev.off()
}


J = c(1, 2, 50, 78, 100)
for(j in J){
png(paste0("C:/Users/jsnow/Desktop/originalSimulationChains/omega_", j, ".png"))
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$omega[,j,i], type = 'l', main = i)
  abline(h = omega[j,i], col = "blue")
}
dev.off()
}



j=2
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$theta[,j,i], type = 'l', main = i)
  abline(h = theta[j, i], col = "blue")
}


for(j in 1:6){
png(paste0("C:/Users/jsnow/Desktop/originalSimulationChains/sigma_", j, ".png"))
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$Sig[,j,i], type = 'l', main = i)
  abline(h = Sigma[j, i], col = "blue")
}
dev.off()
}


png(paste0("C:/Users/jsnow/Desktop/originalSimulationChains/tau.png"))
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$tau[,i], type = 'l', main = i)
  abline(h = tau[i], col = "blue")
}
dev.off()








gen_lambda = function(beta, theta, dat){
  
  
  # relative abundance
  gamma = dat$X %*% beta
  
  return(lambda)
  
}

LT = gen_lambda(beta, theta, dat_train)

LE = array(NA, dim = c(dim(chains$beta)[1], dim(LT)))

for(i in 1:dim(chains$beta)[1]){
  LE[i,,] = gen_lambda(chains$beta[i,,], chains$THETA[i,,], dat_train)
}
LE = apply(LE, c(2,3), mean)

plot(LT[,1], LE[,1])
head(LT)
head(LE)


LE = exp(log(apply(chains$Etilde, c(2,3), mean)) + dat_train$X %*% t(b_hat))

# check coverage ----------------------------------------------------------

# naive coverage
bupper_naive = apply(chains_naive$beta, c(2,3), quantile, probs = 0.975)
blower_naive = apply(chains_naive$beta, c(2,3), quantile, probs = 0.025)

oupper_naive = apply(chains_naive$omega, c(2,3), quantile, probs = 0.975)
olower_naive = apply(chains_naive$omega, c(2,3), quantile, probs = 0.025)

supper_naive = apply(chains_naive$Sig, c(2,3), quantile, probs = 0.975)
slower_naive = apply(chains_naive$Sig, c(2,3), quantile, probs = 0.025)

tupper_naive = apply(chains_naive$tau, c(2), quantile, probs = 0.975)
tlower_naive = apply(chains_naive$tau, c(2), quantile, probs = 0.025)


(beta > blower_naive) & (beta < bupper_naive)
sign(blower_naive) == sign(bupper_naive)
covBetaNaive = sum((beta > blower_naive) & (beta < bupper_naive))/length(beta)

(omega > olower_naive) & (omega < oupper_naive)
sign(olower_naive) == sign(oupper_naive)
covOmegaNaive = sum((omega > olower_naive) & (omega < oupper_naive))/length(omega)

(Sigma > slower_naive) & (Sigma < supper_naive)
sign(slower_naive) == sign(supper_naive)
# don't include diagonal and is symmetric
covSigmaNaive = sum(((Sigma > slower_naive) & (Sigma < supper_naive))[upper.tri((Sigma > slower_naive) & (Sigma < supper_naive))])/length((((Sigma > slower_naive) & (Sigma < supper_naive))[upper.tri((Sigma > slower_naive) & (Sigma < supper_naive))]))


(tau > tlower_naive) & (tau < tupper_naive)
sign(tlower_naive) == sign(tupper_naive)
covTauNaive = sum((tau > tlower_naive) & (tau < tupper_naive))/length(tau)


# our model coverage
bupper = apply(chains$beta, c(2,3), quantile, probs = 0.975)
blower = apply(chains$beta, c(2,3), quantile, probs = 0.025)

pupper = apply(chains$phi, c(2,3), quantile, probs = 0.975)
plower = apply(chains$phi, c(2,3), quantile, probs = 0.025)

oupper = apply(chains$omega, c(2,3), quantile, probs = 0.975)
olower = apply(chains$omega, c(2,3), quantile, probs = 0.025)

supper = apply(chains$Sig, c(2,3), quantile, probs = 0.975)
slower = apply(chains$Sig, c(2,3), quantile, probs = 0.025)

tupper = apply(chains$tau, c(2), quantile, probs = 0.975)
tlower = apply(chains$tau, c(2), quantile, probs = 0.025)



(beta > blower) & (beta < bupper)
sign(blower) == sign(bupper)
covBeta = sum((beta > blower) & (beta < bupper))/length(beta)

(phi > plower) & (phi < pupper)
sign(plower) == sign(pupper)
covPhi = sum((phi > plower) & (phi < pupper))/length(phi)

(omega > olower) & (omega < oupper)
sign(olower) == sign(oupper)
covOmega = sum((omega > olower) & (omega < oupper))/length(omega)


(Sigma > slower) & (Sigma < supper)
sign(slower) == sign(supper)
# don't include diagonal and is symmetric
covSigma = sum(((Sigma > slower) & (Sigma < supper))[upper.tri((Sigma > slower) & (Sigma < supper))])/length((((Sigma > slower) & (Sigma < supper))[upper.tri((Sigma > slower) & (Sigma < supper))]))


(tau > tlower) & (tau < tupper)
sign(tlower) == sign(tupper)
covTau = sum((tau > tlower) & (tau < tupper))/length(tau)



# coverage percentage
sum(covBeta*length(beta),
    covPhi*length(phi),
    covOmega*length(omega),
    covSigma*15,
    covTau*length(tau))/sum(length(beta),length(phi),length(omega),15,length(tau))

# 0.9578059

sum(covBetaNaive*length(beta),
    covOmegaNaive*length(omega),
    covSigmaNaive*15,
    covTauNaive*length(tau))/sum(length(beta),length(omega),15,length(tau))

# 0.6682316



# plot of parameter values compared to true -------------------------------


betaDiff = array(NA, dim = dim(chains$beta))
betaDiffNaive = array(NA, dim = dim(chains_naive$beta))
phiDiff = array(NA, dim = dim(chains$phi))
  
for(i in 1:1500){
  betaDiff[i,,] = chains$beta[i,,] - beta
  betaDiffNaive[i,,] = chains_naive$beta[i,,] - beta
  phiDiff[i,,] = chains$phi[i,,] - phi
}



bupper_naive = apply(betaDiffNaive, c(2,3), quantile, probs = 0.975)
blower_naive = apply(betaDiffNaive, c(2,3), quantile, probs = 0.025)
b_naive = t(apply(betaDiffNaive, c(2,3), mean))

bupper = apply(betaDiff, c(2,3), quantile, probs = 0.975)
blower = apply(betaDiff, c(2,3), quantile, probs = 0.025)
b_hat = t(apply(betaDiff, c(2,3), mean))

pupper = apply(phiDiff, c(2,3), quantile, probs = 0.975)
plower = apply(phiDiff, c(2,3), quantile, probs = 0.025)
p_hat = t(apply(phiDiff, c(2,3), mean))

b_names = colnames(dat_train$X)
b_names[1] = "Int"
phi_names = colnames(dat_train$Z)

b_names = c("intercept", "secchi disk depth", "maximum depth")
phi_names = c("TN", "temp", "sine", "cosine", "temp x sine", "temp x cosine", 
              "GN", "GN x temp", "GN x sine", "GN x cosine", "GN x temp x sine", "GN x temp x cosine")




betaDFvary = as_tibble(b_hat, .name_repair = ~b_names) %>% 
  mutate(fish = fnames) %>% 
  pivot_longer(-fish, names_to = "par", values_to = "mean") %>% 
  left_join(as_tibble(t(bupper), .name_repair = ~b_names) %>% 
              mutate(fish = fnames) %>% 
              pivot_longer(-fish, names_to = "par", values_to = "upper")) %>% 
  left_join(as_tibble(t(blower), .name_repair = ~b_names) %>% 
              mutate(fish = fnames) %>% 
              pivot_longer(-fish, names_to = "par", values_to = "lower")) %>% 
  mutate(model = "varying catchability")

betaDFconst = as_tibble(b_naive, .name_repair = ~b_names) %>% 
  mutate(fish = fnames) %>% 
  pivot_longer(-fish, names_to = "par", values_to = "mean") %>% 
  left_join(as_tibble(t(bupper_naive), .name_repair = ~b_names) %>% 
              mutate(fish = fnames) %>% 
              pivot_longer(-fish, names_to = "par", values_to = "upper")) %>% 
  left_join(as_tibble(t(blower_naive), .name_repair = ~b_names) %>% 
              mutate(fish = fnames) %>% 
              pivot_longer(-fish, names_to = "par", values_to = "lower")) %>% 
  mutate(model = "constant catchability")


betaDF = rbind(betaDFvary, betaDFconst)


phiDF = as_tibble(p_hat, .name_repair = ~phi_names) %>% 
  mutate(fish = fnames) %>% 
  pivot_longer(-fish, names_to = "par", values_to = "mean") %>% 
  left_join(as_tibble(t(pupper), .name_repair = ~phi_names) %>% 
              mutate(fish = fnames) %>% 
              pivot_longer(-fish, names_to = "par", values_to = "upper")) %>% 
  left_join(as_tibble(t(plower), .name_repair = ~phi_names) %>% 
              mutate(fish = fnames) %>% 
              pivot_longer(-fish, names_to = "par", values_to = "lower")) %>% 
  mutate(par = factor(par, levels = phi_names))




ggplot(betaDF, aes(x = mean, y = fish, color = model)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.3) +
  geom_vline(xintercept = 0) +
  facet_wrap(~par, scales = "free_x") +
  scale_color_manual(values = c("red", "blue"), name = "") +
  ylab("") +
  xlab("") +
  theme(legend.position="bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'),
        axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1),
        plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))

# ggsave('C:/Users/jsnow/Desktop/originalSimulationChains/betaEstimates.eps', dpi = 600, height = 4, width = 12, device=grDevices::cairo_ps)
# ggsave("C:/Users/jsnow/Desktop/originalSimulationChains/betaEstimates.png", width = 12, height = 4)



ggplot(phiDF, aes(x = mean, y = fish)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.3) +
  geom_vline(xintercept = 0) +
  facet_wrap(~par, scales = "free_x", nrow = 6, dir = "v") +
  ylab("") +
  xlab("") +
  theme(legend.position="bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'),
        axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1),
        panel.spacing = unit(1, "lines"),
        plot.margin = unit(c(1,1,1,1), "cm"))
  
# ggsave('C:/Users/jsnow/Desktop/originalSimulationChains/phiEstimates.eps', dpi = 600, width = 8, height = 14, device=grDevices::cairo_ps)
# ggsave("C:/Users/jsnow/Desktop/originalSimulationChains/phiEstimates.png", width = 8, height = 14)



# rmse in sample ----------------------------------------------------------


computeGamma = function(beta, omega, dat){
  
  OMEGA = dat$each_lake %*% omega
  
  # relative abundance
  gamma = dat$X %*% beta + OMEGA
  
  return(exp(gamma))
  
}


gammaTrue = computeGamma(beta, omega, dat_train)

nlks = dim(gammaTrue)[1]
nsmps = dim(chains$beta)[1]

gamHat = array(NA, dim = c(nlks, 6, nsmps))
gamHat_naive = array(NA, dim = c(nlks, 6, nsmps))

for(i in 1:nsmps){
  
  gamHat[,,i] = computeGamma(chains$beta[i,,], chains$omega[i,,], dat_train)
  gamHat_naive[,,i] = computeGamma(chains_naive$beta[i,,], chains$omega[i,,], dat_train)
  
}


dSqrd = array(NA, dim = dim(gamHat)[2:3])
dSqrd_naive = array(NA, dim = dim(gamHat_naive)[2:3])

for(i in 1:nsmps){
  dSqrd[,i] = apply((gamHat[,,i] - gammaTrue)^2, 2, function(x){sqrt(mean(x))})
  dSqrd_naive[,i] = apply((gamHat_naive[,,i] - gammaTrue)^2, 2, function(x){sqrt(mean(x))})
}



apply(dSqrd, 1, mean)
apply(dSqrd_naive, 1, mean)


fnames = dat_train$fish_names
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
  xlab("Predicted Relative Abundance") +
  ylab("Observed Relative Abundance") +
  ggtitle("Our Model")


p2 = ggplot(preds_naive, aes(x = Predicted, y = Observed)) +
  geom_point() +
  scale_color_manual(values = c('purple', 'forest green'), labels = c('9', '18'), name = "") +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fish, scales = 'free') +
  xlab("Predicted Relative Abundance") +
  ylab("Observed Relative Abundance") +
  ggtitle("Naive Model")

# ggsave("C:/Users/jsnow/Desktop/originalSimulationChains/insamplePredictions.png", cowplot::plot_grid(p1, p2), width = 12, height = 6)


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


# relative abundance check ------------------------------------------------

rep(1:5, each = 3)

tmp = as_tibble(apply(gamHat, c(1,2), median), .name_repair = ~fnames) %>% 
  mutate(ID = 1:n(),
         lakeID = rep(1:(n()/20), each = 20)) %>% 
  pivot_longer(-c(ID, lakeID), names_to = "fish", values_to = "Predicted") %>% 
  left_join(as_tibble(gammaTrue, .name_repair = ~fnames) %>% 
              mutate(ID = 1:n(),
                     lakeID = rep(1:(n()/20), each = 20)) %>% 
              pivot_longer(-c(ID, lakeID), names_to = "fish", values_to = "Observed"),
            by = c("ID", "fish", "lakeID")) %>% 
  group_by(fish) %>% 
  distinct(lakeID, .keep_all = T) %>% 
  select(-ID) %>% 
  pivot_longer(c(Predicted, Observed), names_to = "name", values_to = "val") %>% 
  pivot_wider(id_cols = c(lakeID, name), names_from = fish, values_from = val) 

write_csv(tmp %>% filter(lakeID %in% sample(1:200, size = 5)),"C:/Users/jsnow/Desktop/originalSimulationChains/withinLakeComparison.csv")


# prediction data set -----------------------------------------------------


fish_pred = fish_sim %>%
  filter(year == 2019)

dat_test = create_pars(fish_pred,
                       mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs,
                       center_temp_pred)

# predictability of lambda ------------------------------------------------


computeGamma = function(beta, omega, dat){
  
  OMEGA = dat$each_lake %*% omega
  
  # relative abundance
  gamma = dat$X %*% beta + OMEGA
  
  return(exp(gamma))
  
}


gammaTrue = computeGamma(beta, omega, dat_test)

nlks = dim(gammaTrue)[1]
nsmps = dim(chains$beta)[1]

gamHat = array(NA, dim = c(nlks, 6, nsmps))
gamHat_naive = array(NA, dim = c(nlks, 6, nsmps))

for(i in 1:nsmps){
  
  gamHat[,,i] = computeGamma(chains$beta[i,,], chains$omega[i,,], dat_test)
  gamHat_naive[,,i] = computeGamma(chains_naive$beta[i,,], chains_naive$omega[i,,], dat_test)
  
}


dSqrd = array(NA, dim = dim(gamHat)[2:3])
dSqrd_naive = array(NA, dim = dim(gamHat_naive)[2:3])

for(i in 1:nsmps){
  dSqrd[,i] = apply((gamHat[,,i] - gammaTrue)^2, 2, function(x){sqrt(mean(x))})
  dSqrd_naive[,i] = apply((gamHat_naive[,,i] - gammaTrue)^2, 2, function(x){sqrt(mean(x))})
}



apply(dSqrd, 1, mean)
apply(dSqrd_naive, 1, mean)


fnames = dat_train$fish_names
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


p1 = ggplot(preds, aes(x = Predicted, y = Observed)) +
  geom_point() +
  scale_color_manual(values = c('purple', 'forest green'), labels = c('9', '18'), name = "") +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fish, scales = 'free') +
  xlab("Predicted Relative Abundance") +
  ylab("Observed Relative Abundance") +
  ggtitle("Our Model")


p2 = ggplot(preds_naive, aes(x = Predicted, y = Observed)) +
  geom_point() +
  scale_color_manual(values = c('purple', 'forest green'), labels = c('9', '18'), name = "") +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fish, scales = 'free') +
  xlab("Predicted Relative Abundance") +
  ylab("Observed Relative Abundance") +
  ggtitle("Naive Model")


cowplot::plot_grid(p1, p2)

# ggsave("/Users/jsnow/Desktop/outofsample.pdf", cowplot::plot_grid(p1, p2), width = 12, height = 6)

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


# OOS relative abundance --------------------------------------------------


gamHatRel = array(NA, dim = c(nlks/2, 6, nsmps))
gamHatRel_naive = array(NA, dim = c(nlks/2, 6, nsmps))

for(i in 1:nsmps){
  
  gamHatRel[,,i] = gamHat[cumsum(apply(t(dat_test$each_lake), 1, sum)),,i]
  gamHatRel_naive[,,i] = gamHat_naive[cumsum(apply(t(dat_test$each_lake), 1, sum)),,i]
  
}


# DOW = dat_test$lake_id
gamRel = as_tibble(apply(gamHatRel, c(1,2), median), .name_repair = ~dat_test$fish_names) %>% 
  mutate(DOW = dat_test$lake_id) %>% 
  pivot_longer(-DOW, names_to = "fish", values_to = "mean") %>% 
  left_join(as_tibble(apply(gamHatRel, c(1,2), quantile, probs = 0.025), .name_repair = ~dat_test$fish_names) %>% 
              mutate(DOW = dat_test$lake_id) %>% 
              pivot_longer(-DOW, names_to = "fish", values_to = "lower"), by = c("DOW", "fish")) %>% 
  left_join(as_tibble(apply(gamHatRel, c(1,2), quantile, probs = 0.975), .name_repair = ~dat_test$fish_names) %>% 
              mutate(DOW = dat_test$lake_id) %>% 
              pivot_longer(-DOW, names_to = "fish", values_to = "upper"), by = c("DOW", "fish")) %>% 
  left_join(as_tibble(gammaTrue[cumsum(apply(t(dat_test$each_lake), 1, sum)),], .name_repair = ~dat_test$fish_names) %>% 
              mutate(DOW = dat_test$lake_id) %>% 
              pivot_longer(-DOW, names_to = "fish", values_to = "true"), by = c("DOW", "fish"))

gamRelNaive = as_tibble(apply(gamHatRel_naive, c(1,2), median), .name_repair = ~dat_test$fish_names) %>% 
  mutate(DOW = dat_test$lake_id) %>% 
  pivot_longer(-DOW, names_to = "fish", values_to = "mean") %>% 
  left_join(as_tibble(apply(gamHatRel_naive, c(1,2), quantile, probs = 0.025), .name_repair = ~dat_test$fish_names) %>% 
              mutate(DOW = dat_test$lake_id) %>% 
              pivot_longer(-DOW, names_to = "fish", values_to = "lower"), by = c("DOW", "fish")) %>% 
  left_join(as_tibble(apply(gamHatRel_naive, c(1,2), quantile, probs = 0.975), .name_repair = ~dat_test$fish_names) %>% 
              mutate(DOW = dat_test$lake_id) %>% 
              pivot_longer(-DOW, names_to = "fish", values_to = "upper"), by = c("DOW", "fish")) %>% 
  left_join(as_tibble(gammaTrue[cumsum(apply(t(dat_test$each_lake), 1, sum)),], .name_repair = ~dat_test$fish_names) %>% 
              mutate(DOW = dat_test$lake_id) %>% 
              pivot_longer(-DOW, names_to = "fish", values_to = "true"), by = c("DOW", "fish"))


gamRel = gamRel %>% 
  group_by(fish) %>% 
  mutate(mean_rank = dense_rank(mean),
         lower_rank = dense_rank(lower),
         upper_rank = dense_rank(upper),
         true_rank = dense_rank(true)) %>% 
  ungroup()

gamRelNaive = gamRelNaive %>% 
  group_by(fish) %>% 
  mutate(mean_rank = dense_rank(mean),
         true_rank = dense_rank(true)) %>% 
  ungroup()


rsqured = gamRel %>% 
  group_by(fish) %>% 
  summarize(r = summary(lm(mean_rank ~ true_rank))$r.squared)

rsquredNaive = gamRelNaive %>% 
  group_by(fish) %>% 
  summarize(r = summary(lm(mean_rank ~ true_rank))$r.squared)


xtable::xtable(as_tibble(rsqured) %>% 
  rename("varying catchability" = "r") %>% 
  left_join(as_tibble(rsquredNaive) %>% 
          rename("constant catchability" = "r")) %>% 
    rename("species" = "fish"), digits = 4)


ggplot(gamRel, aes(x = true_rank, y = mean_rank)) +
  geom_point() +
  facet_wrap(~fish) +
  xlab("True Rank") +
  ylab("Predicted Rank") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))


# ggsave('C:/Users/jsnow/Desktop/originalSimulationChains/LakeRanking.eps', dpi = 600, height = 8, width = 12, device=grDevices::cairo_ps)
# ggsave("C:/Users/jsnow/Desktop/originalSimulationChains/LakeRanking.png", width = 12, height = 8)


ggplot(gamRelNaive, aes(x = true_rank, y = mean_rank)) +
  geom_point() +
  facet_wrap(~fish) +
  xlab("True Rank") +
  ylab("Predicted Rank") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))


# ggsave('C:/Users/jsnow/Desktop/originalSimulationChains/LakeRankingNaive.eps', dpi = 600, height = 8, width = 12, device=grDevices::cairo_ps)
# ggsave("C:/Users/jsnow/Desktop/originalSimulationChains/LakeRankingNaive.png", width = 12, height = 8)



gamRatTrue = t(apply(gammaTrue[cumsum(apply(t(dat_test$each_lake), 1, sum)),], 1, function(x) x/sum(x))) 

gamHatRat = array(NA, dim = c(nlks/2, 6, nsmps))
gamHatRat_naive = array(NA, dim = c(nlks/2, 6, nsmps))

gamHatRatDiff = array(NA, dim = c(nlks/2, 6, nsmps))
gamHatRatDiff_naive = array(NA, dim = c(nlks/2, 6, nsmps))

for(i in 1:nsmps){
  
  gamHatRat[,,i] = t(apply(gamHatRel[,,i], 1, function(x) x/sum(x)))
  gamHatRat_naive[,,i] = t(apply(gamHatRel_naive[,,i], 1, function(x) x/sum(x)))
  
  
  # gamHatRatDiff[,,i] = t(apply(gammaTrue[cumsum(apply(t(dat_test$each_lake), 1, sum)),], 1, function(x) x/sum(x))) - gamHatRat[,,i]
  # gamHatRatDiff_naive[,,i] = t(apply(gammaTrue[cumsum(apply(t(dat_test$each_lake), 1, sum)),], 1, function(x) x/sum(x))) - gamHatRat_naive[,,i]
  
  gamHatRatDiff[,,i] = gamRatTrue - gamHatRat[,,i]
  gamHatRatDiff_naive[,,i] = gamRatTrue - gamHatRat_naive[,,i]
  
}


gamRat = t(apply(gammaTrue[cumsum(apply(t(dat_test$each_lake), 1, sum)),], 1, function(x) x/sum(x)))


crpsHat = array(NA, dim = c(nlks/2, 6))
crpsHat_naive = array(NA, dim = c(nlks/2, 6))

mseHat = array(NA, dim = c(nlks/2, 6))
mseHat_naive = array(NA, dim = c(nlks/2, 6))

for(i in 1:(nlks/2)){
  for( j in 1:6){
    crpsHat[i,j] = calc_crps(gamHatRat[i,j,], gamRatTrue[i,j])
    crpsHat_naive[i,j] = calc_crps(gamHatRat_naive[i,j,], gamRatTrue[i,j])
    
    mseHat[i,j] = mean((gamHatRat[i,j,] - gamRatTrue[i,j])^2)
    mseHat_naive[i,j] = mean((gamHatRat_naive[i,j,] - gamRatTrue[i,j])^2)
  }
}

tibble(names = dat_test$fish_names, our = apply(crpsHat, 2, mean), naive = apply(crpsHat_naive, 2, mean))
tibble(names = dat_test$fish_names, our = apply(mseHat, 2, mean), naive = apply(mseHat_naive, 2, mean))



lower = 0.025
upper = 0.975

gamRat = as_tibble(apply(gamHatRatDiff, c(1,2), median), .name_repair = ~dat_test$fish_names) %>% 
  mutate(DOW = dat_test$lake_id) %>% 
  pivot_longer(-DOW, names_to = "fish", values_to = "mean") %>% 
  left_join(as_tibble(apply(gamHatRatDiff, c(1,2), quantile, probs = lower), .name_repair = ~dat_test$fish_names) %>% 
              mutate(DOW = dat_test$lake_id) %>% 
              pivot_longer(-DOW, names_to = "fish", values_to = "lower"), by = c("DOW", "fish")) %>% 
  left_join(as_tibble(apply(gamHatRatDiff, c(1,2), quantile, probs = upper), .name_repair = ~dat_test$fish_names) %>% 
              mutate(DOW = dat_test$lake_id) %>% 
              pivot_longer(-DOW, names_to = "fish", values_to = "upper"), by = c("DOW", "fish")) %>% 
  mutate(CI = ifelse(lower < 0 & upper > 0, 1, 0),
         CI = factor(CI, levels = c(0, 1))) %>% 
  group_by(fish) %>% 
  mutate(ID = 1:n()) %>% 
  ungroup()

gamRatNaive = as_tibble(apply(gamHatRatDiff_naive, c(1,2), median), .name_repair = ~dat_test$fish_names) %>% 
  mutate(DOW = dat_test$lake_id) %>% 
  pivot_longer(-DOW, names_to = "fish", values_to = "mean") %>% 
  left_join(as_tibble(apply(gamHatRatDiff_naive, c(1,2), quantile, probs = lower), .name_repair = ~dat_test$fish_names) %>% 
              mutate(DOW = dat_test$lake_id) %>% 
              pivot_longer(-DOW, names_to = "fish", values_to = "lower"), by = c("DOW", "fish")) %>% 
  left_join(as_tibble(apply(gamHatRatDiff_naive, c(1,2), quantile, probs = upper), .name_repair = ~dat_test$fish_names) %>% 
              mutate(DOW = dat_test$lake_id) %>% 
              pivot_longer(-DOW, names_to = "fish", values_to = "upper"), by = c("DOW", "fish")) %>% 
  mutate(CI = ifelse(lower < 0 & upper > 0, 1, 0),
         CI = factor(CI, levels = c(0, 1))) %>% 
  group_by(fish) %>% 
  mutate(ID = 1:n()) %>% 
  ungroup()


ggplot(gamRat, aes(x = ID, y = mean)) +
  geom_point() +
  geom_errorbar(aes(x = ID, ymin = lower, ymax = upper)) +
  facet_wrap(~fish, scales = "free") +
  # scale_color_manual(values = c("red", "blue"), labels = c("No", "Yes"), name = "Credible Interval Covers Zero:") +
  ylab("True Minus Estimated Percent") +
  xlab("Lake ID") +
  theme(legend.position="bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave("C:/Users/jsnow/Desktop/originalSimulationChains/LakeComposition.png", width = 12, height = 8)
# ggsave('C:/Users/jsnow/Desktop/originalSimulationChains/LakeComposition.eps', dpi = 600, height = 8, width = 12, device=grDevices::cairo_ps)


ggplot(gamRatNaive, aes(x = ID, y = mean, color = CI)) +
  geom_point() +
  geom_errorbar(aes(x = ID, ymin = lower, ymax = upper)) +
  facet_wrap(~fish, scales = "free") +
  scale_color_manual(values = c("red", "black"), labels = c("No", "Yes"), name = "Credible Interval Covers Zero:") +
  ylab("True Minus Estimated Percent") +
  xlab("Lake ID") +
  theme(legend.position="bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave("C:/Users/jsnow/Desktop/originalSimulationChains/LakeCompositionNaive.png", width = 12, height = 8)
# ggsave('C:/Users/jsnow/Desktop/originalSimulationChains/LakeCompositionNaive.eps', dpi = 600, height = 8, width = 12, device=grDevices::cairo_ps)


# 
# gamRat = as_tibble(apply(gamHatRat, c(1,2), mean), .name_repair = ~dat_test$fish_names) %>%
#   mutate(DOW = dat_test$lake_id) %>%
#   pivot_longer(-DOW, names_to = "fish", values_to = "mean") %>%
#   left_join(as_tibble(apply(gamHatRat, c(1,2), quantile, probs = 0.025), .name_repair = ~dat_test$fish_names) %>%
#               mutate(DOW = dat_test$lake_id) %>%
#               pivot_longer(-DOW, names_to = "fish", values_to = "lower"), by = c("DOW", "fish")) %>%
#   left_join(as_tibble(apply(gamHatRat, c(1,2), quantile, probs = 0.975), .name_repair = ~dat_test$fish_names) %>%
#               mutate(DOW = dat_test$lake_id) %>%
#               pivot_longer(-DOW, names_to = "fish", values_to = "upper"), by = c("DOW", "fish")) %>%
#   left_join(as_tibble(t(apply(gammaTrue[cumsum(apply(t(dat_test$each_lake), 1, sum)),], 1, function(x) x/sum(x))), .name_repair = ~dat_test$fish_names) %>%
#               mutate(DOW = dat_test$lake_id) %>%
#               pivot_longer(-DOW, names_to = "fish", values_to = "true"), by = c("DOW", "fish"))
# 
# gamRatNaive = as_tibble(apply(gamHatRat_naive, c(1,2), mean), .name_repair = ~dat_test$fish_names) %>%
#   mutate(DOW = dat_test$lake_id) %>%
#   pivot_longer(-DOW, names_to = "fish", values_to = "mean") %>%
#   left_join(as_tibble(apply(gamHatRat_naive, c(1,2), quantile, probs = 0.025), .name_repair = ~dat_test$fish_names) %>%
#               mutate(DOW = dat_test$lake_id) %>%
#               pivot_longer(-DOW, names_to = "fish", values_to = "lower"), by = c("DOW", "fish")) %>%
#   left_join(as_tibble(apply(gamHatRat_naive, c(1,2), quantile, probs = 0.975), .name_repair = ~dat_test$fish_names) %>%
#               mutate(DOW = dat_test$lake_id) %>%
#               pivot_longer(-DOW, names_to = "fish", values_to = "upper"), by = c("DOW", "fish")) %>%
#   left_join(as_tibble(t(apply(gammaTrue[cumsum(apply(t(dat_test$each_lake), 1, sum)),], 1, function(x) x/sum(x))), .name_repair = ~dat_test$fish_names) %>%
#               mutate(DOW = dat_test$lake_id) %>%
#               pivot_longer(-DOW, names_to = "fish", values_to = "true"), by = c("DOW", "fish"))
# 
# 
# gamRat = gamRat %>%
#   mutate(CI = if_else((true-lower > 0) & (true - upper < 0), 1, 0),
#          CI = factor(CI, levels = c(0, 1))) %>%
#   group_by(fish) %>%
#   mutate(ID = 1:n()) %>%
#   ungroup()
# 
# gamRatNaive = gamRatNaive %>%
#   mutate(CI = if_else((true-lower > 0) & (true - upper < 0), 1, 0),
#          CI = factor(CI, levels = c(0, 1))) %>%
#   group_by(fish) %>%
#   mutate(ID = 1:n()) %>%
#   ungroup()
# 



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

