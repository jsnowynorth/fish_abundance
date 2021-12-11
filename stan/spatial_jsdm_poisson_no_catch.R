library(tidyverse)
library(rstan)
library(fields)
library(Matrix)


# load data ---------------------------------------------------------------

fish_dat = read_csv('data/fish_dat.csv')

fish_dat = fish_dat %>% 
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME)) %>% 
  arrange(DOW, year, COMMON_NAME)

# mean_covs = colnames(fish_dat)[c(7, 9, 23, 24, 13:15,17)]
# temporal_covs = colnames(fish_dat)[c(23, 24)]
# mean_covs_log = colnames(fish_dat)[c(7, 9)]
# mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
# catch_covs = colnames(fish_dat)[c(25, 27:30)] # temp doy interaction
# gear_types = colnames(fish_dat)[c(21, 22)]


mean_covs = colnames(fish_dat)[c(7, 9, 23, 24, 13:15)]
temporal_covs = colnames(fish_dat)[c(23, 24)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15)]
catch_covs = colnames(fish_dat)[c(25, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]



# set up stan parameters --------------------------------------------------

create_pars <- function(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs){
  
  
  fish_dat = fish_dat %>% 
    arrange(DOW, SURVEYDATE, COMMON_NAME, GN)
  
  K = nlevels(fish_dat$COMMON_NAME)
  levs = levels(fish_dat$COMMON_NAME)
  lake_index = (fish_dat %>% filter(COMMON_NAME == levs[1]))$DOW # lake_index == lake_id[2]
  lake_id = levels(fish_dat$DOW)
  n_lakes = length(levels(fish_dat$DOW))
  n_obs = (fish_dat %>% filter(COMMON_NAME == levs[1]) %>% select(TOTAL_CATCH))$TOTAL_CATCH
  
  pars = list()
  
  # data
  # pars$Y = list()
  # pars$X = list()
  # pars$Z = list()
  # pars$effort = list()
  
  
  # projection matrix
  P = fish_dat %>%
    distinct(DOW, .keep_all = T) %>%
    select(all_of(mean_covs), DOW) %>%
    select(-all_of(temporal_covs), DOW) %>%
    left_join(fish_dat %>%
                select(DOW, all_of(temporal_covs)) %>%
                group_by(DOW) %>%
                summarise_all(mean), by = 'DOW') %>%
    select(-DOW) %>%
    mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>%
    mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>%
    mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>%
    mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>%
    # mutate(secchi = secchi - mean(secchi)) %>%
    as.matrix()
  
  pars$P = diag(nrow(P)) - P %*% solve(t(P) %*% P) %*% t(P)
  
  
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
    mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN))
  
  pars$X = as.matrix(X)
  pars$Z = as.matrix(Z)
  pars$y_vec = c(pars$Y)
  
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

dat = create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs)

out = stan(file = 'stan/spatial_jdsm_poisson_no_catch.stan', data = dat, iter = 2000, warmup = 1000, chains = 1, cores = 1) # spatial cholesky
saveRDS(out, '/Users/joshuanorth/Desktop/stan_lake_random_effect_no_catch.rds')

# out = read_rds('/Users/joshuanorth/Desktop/stan_lake_random_effect_no_catch.rds')

chains = extract(out, permuted = T)
names(chains)
lapply(chains, dim)

b = t(apply(chains$beta, c(2,3), mean))
b_lower = t(apply(chains$beta, c(2,3), quantile, probs = 0.025))
b_upper = t(apply(chains$beta, c(2,3), quantile, probs = 0.975))
colnames(b) = colnames(b_lower) = colnames(b_upper) = colnames(head(dat$X))
rownames(b) = rownames(b_lower) = rownames(b_upper) = dat$fish_names
int = apply(chains$beta_0, 2, mean)
int_lower = apply(chains$beta_0, 2, quantile, probs = 0.025)
int_upper = apply(chains$beta_0, 2, quantile, probs = 0.975)
b = cbind(b, int)
b_lower = cbind(b_lower, int_lower)
b_upper = cbind(b_upper, int_upper)
b
b_lower
b_upper

sign(b_lower) == sign(b_upper)

# saveRDS(out, '/Users/joshuanorth/Desktop/stan_spatial.rds')

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$beta_0[,i], type = 'l', main = i)
}

par(mfrow = c(2,4))
for(i in 1:8){
  plot(chains$beta[,i,5], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$omega[,i,3], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$Sigma_species[,i,1], type = 'l', main = i)
}

# image.plot(apply(chains$Sigma_species, c(2,3), mean), col = two.colors(start = 'blue', end = 'red', middle = 'white'), zlim = c(-0.2, 0.2))

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$tau[,i], type = 'l', main = i)
}


plot(chains$lp__, type = 'l')


