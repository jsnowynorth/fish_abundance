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

fish_dat = fish_dat %>% 
  group_by(DOW) %>% 
  mutate(secchi = median(secchi)) %>% 
  ungroup()

fish_dat = fish_dat %>% 
  group_by(DOW) %>% 
  mutate(secchi_base = secchi[which.min(year)],
         DD5_base = DD5[which.min(year)],
         year_base = min(year)) %>% 
  mutate(secchi_trend = ifelse(is.na((secchi - secchi_base)/(year - year_base)), 0, (secchi - secchi_base)/(year - year_base)),
         DD5_trend = ifelse(is.na((DD5 - DD5_base)/(year - year_base)), 0, (DD5 - DD5_base)/(year - year_base))) %>% 
  ungroup()

# fish_dat %>% 
#   group_by(year) %>% 
#   summarise(n = n()) %>% 
#   arrange(n)
# 
# fish_dat = fish_dat %>% 
#   filter(year == 2008) %>% 
#   mutate(DOW = droplevels(DOW))

# fish_dat %>% 
#   filter(COMMON_NAME == 'bluegill' | COMMON_NAME == 'walleye') %>% 
#   mutate(DOW = droplevels(DOW),
#          COMMON_NAME = droplevels(COMMON_NAME))

# original
mean_covs = colnames(fish_dat)[c(7, 9, 23, 24, 13:15,17)]
temporal_covs = colnames(fish_dat)[c(23, 24)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
catch_covs = colnames(fish_dat)[c(25, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]

# with secchi and dd5 trend
# mean_covs = colnames(fish_dat)[c(7, 9, 31,32,34,35, 13:15,17)]
# temporal_covs = colnames(fish_dat)[c(23, 25)]
# mean_covs_log = colnames(fish_dat)[c(7, 9)]
# mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
# catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
# gear_types = colnames(fish_dat)[c(21, 22)]

mean_covs = colnames(fish_dat)[c(7, 9, 23, 24)]
temporal_covs = colnames(fish_dat)[c(23, 24)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = NULL
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
    mutate_at(vars(all_of(mean_covs_log)), ~ . - mean(.))
    # mutate(DS = DD5*secchi)
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

out = stan(file = 'stan/jsdm_poisson.stan', data = dat, iter = 1000, warmup = 500, chains = 1, cores = 1)
saveRDS(out, '/Users/joshuanorth/Desktop/stan_results/stan_non_spatial_new_temp_median_secchi.rds')
# saveRDS(out, '/Users/joshuanorth/Desktop/stan_non_spatial_no_land_no_secchi_no_temp.rds')

# out = read_rds('/Users/joshuanorth/Desktop/stan_non_spatial_trend.rds')
# out = read_rds('/Users/joshuanorth/Desktop/stan_non_spatial_no_ag.rds')
# out = read_rds('/Users/joshuanorth/Desktop/stan_results/stan_non_spatial_no_land_no_secchi_no_temp.rds')



out = read_rds('/Users/joshuanorth/Desktop/stan_results/stan_non_spatial_new_temp_no_land.rds')
out = read_rds('/Users/joshuanorth/Desktop/stan_results/stan_non_spatial_new_temp.rds')
out = read_rds('/Users/joshuanorth/Desktop/stan_results/stan_non_spatial_trend_new_temp.rds')


chains = extract(out, permuted = T)
names(chains)
lapply(chains, dim)

# saveRDS(out, '/Users/joshuanorth/Desktop/stan_spatial.rds')

b_stan = t(apply(chains$beta, c(2,3), mean))
colnames(b_stan) = colnames(head(dat$X))
rownames(b_stan) = dat$fish_names
int = apply(chains$beta_0, 2, mean)
b_stan = cbind(b_stan, int)

l = t(apply(chains$beta, c(2,3), quantile, probs = 0.025))
u = t(apply(chains$beta, c(2,3), quantile, probs = 0.975))
sign(u) == sign(l)

p_stan = t(apply(chains$phi, c(2,3), mean))
colnames(p_stan) = colnames(head(dat$Z))
rownames(p_stan) = dat$fish_names

l = t(apply(chains$phi, c(2,3), quantile, probs = 0.025))
u = t(apply(chains$phi, c(2,3), quantile, probs = 0.975))
sign(u) == sign(l)


par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$beta_0[,i], type = 'l', main = i)
}

par(mfrow = c(2,5))
for(i in 1:10){
  plot(chains$beta[,i,5], type = 'l', main = i)
}

par(mfrow = c(3,4))
for(i in 1:11){
  plot(chains$phi[,i,1], type = 'l', main = i)
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$omega[,i], type = 'l', main = i)
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



# no spatial original model -----------------------------------------------

# my sampler
run = read_rds('/Users/joshuanorth/Desktop/stan_results/full_model_non_spatial_1.rds')
b_run = apply(run$beta, c(1,2), mean)
colnames(b_run) = colnames(head(run$pars$X[[1]]))
rownames(b_run) = run$pars$fish_names
# int = apply(run$beta_0, 2, mean)

par(mfrow = c(2,4))
for(i in 1:8){
  plot(run$beta[5,i,], type = 'l', main = i)
}

# stan sampler
out = read_rds('/Users/joshuanorth/Desktop/stan_results/stan_non_spatial.rds')
chains = extract(out, permuted = T)
mean_covs = colnames(fish_dat)[c(7, 9, 23, 25, 13:15,17)]
temporal_covs = colnames(fish_dat)[c(23, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]
dat = create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs)

b_stan = t(apply(chains$beta, c(2,3), mean))
colnames(b_stan) = colnames(head(dat$X))
rownames(b_stan) = dat$fish_names
int = apply(chains$beta_0, 2, mean)
b_stan = cbind(b_stan, int)

# stan no ag
out = read_rds('/Users/joshuanorth/Desktop/stan_non_spatial_no_ag.rds')
chains = extract(out, permuted = T)
mean_covs = colnames(fish_dat)[c(7, 9, 23, 25, 14:15,17)]
temporal_covs = colnames(fish_dat)[c(23, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(14:15,17)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]
dat = create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs)

b_stan_no_ag = t(apply(chains$beta, c(2,3), mean))
colnames(b_stan_no_ag) = colnames(head(dat$X))
rownames(b_stan_no_ag) = dat$fish_names
int = apply(chains$beta_0, 2, mean)
b_stan_no_ag = cbind(b_stan_no_ag, int)










par(mfrow = c(2,4))
for(i in 1:8){
  plot(run$beta[,i,5], type = 'l', main = i)
}


par(mfrow = c(2,4))
for(i in 1:8){
  hist(dat$X[,i], main = colnames(dat$X)[i])
}
