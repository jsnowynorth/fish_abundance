library(tidyverse)
library(rstan)
library(fields)
library(Matrix)


# load data ---------------------------------------------------------------

fish_dat = read_csv('/home/jntmf/data/fish/data/fish_dat.csv')

fish_dat = fish_dat %>% 
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME)) %>% 
  arrange(DOW, year, COMMON_NAME)

# 308 lakes
# only lakes observed 4 or more times
fish_dat = fish_dat %>% 
  group_by(DOW) %>% 
  filter(n() >= 4*12) %>% 
  ungroup() %>% 
  mutate_at(vars(DOW), ~droplevels(.))

# length(table(fish_dat$DOW))


mean_covs = colnames(fish_dat)[c(7, 9, 23, 25, 13:15,17)]
temporal_covs = colnames(fish_dat)[c(23, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
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


out = stan(file = 'spatial_jsdm_poisson_cholesky.stan', data = dat, iter = 1000, warmup = 500, chains = 3, cores = 3) # spatial cholesky


# chains = extract(out)

# saveRDS(chains, '/home/jntmf/data/fish/results/spatial/stan_spatial.rds')
saveRDS(out, '/home/jntmf/data/fish/results/spatial/stan_spatial.rds')







