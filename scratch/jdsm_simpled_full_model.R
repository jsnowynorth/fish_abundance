

# load libraries ----------------------------------------------------------

library(tidyverse)
library(rstan)
library(fields)
library(Matrix)
library(lubridate)
library(stringr)



# load data ---------------------------------------------------------------


center_temp = read_csv('data/center_temp.csv')
center_temp = center_temp %>%
  mutate(DOW = as.factor(DOW)) %>%
  arrange(DOW, year)


fish_dat = read_csv("data/fish_dat.csv")
fish_dat = fish_dat %>%
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME)) %>%
  arrange(DOW, year, COMMON_NAME) %>% 
  filter(EFFORT != 0)


fish_dat %>% 
  group_by(DOW) %>% 
  summarise(nReps = n()/12) %>% 
  group_by(nReps) %>% 
  summarise(nLakes = n())


fish_dat = fish_dat %>% 
  right_join(fish_dat %>% 
              group_by(DOW) %>% 
              summarise(nReps = n()/12) %>% 
              filter(nReps >= 5) %>% 
              select(DOW)) %>% 
  arrange(DOW, year, COMMON_NAME) %>% 
  filter(EFFORT != 0) %>% 
  droplevels()
  


# select covariates from data
mean_covs = colnames(fish_dat)[c(7, 9, 23, 24, 14:15)]
temporal_covs = colnames(fish_dat)[c(23, 24)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(14:15)]
catch_covs = colnames(fish_dat)[c(25, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]


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


dat = create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, center_temp)


str(dat)



# look at data distribution -----------------------------------------------

ggplot(fish_dat, aes(x = DOY)) + 
  geom_histogram()


# lats = range(fish_dat$lat, na.rm = T)
# lons = range(fish_dat$lon, na.rm = T)
# 
# usa = st_as_sf(maps::map("state", fill= TRUE, plot = FALSE))
# 
# ggplot() +
#   geom_sf(data = usa) +
#   coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
#   geom_jitter(data = fish_dat, aes(x = lon, y = lat, color = DOY), width = 0.1, height = 0.1)
# 
# 
# ggsave("C:/Users/jsnow/Desktop/FourOrMore.pdf")



# set up stan parameters --------------------------------------------------


out = stan(file = 'scratch/jdsm_simples_v4.stan',
           data = dat, iter = 500, warmup = 250, chains = 3, cores = 3, refresh = 10) # our model

# saveRDS(out, "C:/Users/jsnow/Desktop/new_model_full_run.rds")
out = read_rds("C:/Users/jsnow/Desktop/new_model_full_run.rds")
# out = read_rds("C:/Users/jsnow/Desktop/stan_jsdm_full_v4.rds")


chains = extract(out, permuted = T)

names(chains)
lapply(chains, dim)


apply(chains$beta, c(2,3), mean)
apply(chains$phi, c(2,3), mean)
apply(chains$Sig, c(2,3), mean)
apply(chains$Sigma_species, c(2,3), function(x) x %*% t(x))


j=1
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$beta[,j,i], type = 'l', main = i)
}

j=2
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$phi[,j,i], type = 'l', main = i)
}


j = 80
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$omega[,j,i], type = 'l', main = i)
}

j = 10
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$z[,j,i], type = 'l', main = i)
}

j=1
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$Sig[,j,i], type = 'l', main = i)
}

j=6
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$Sigma_species[,j,i], type = 'l', main = i)
}

j=3
par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$tau[,i], type = 'l', main = i)
}


C = array(NA, dim = dim(chains$Sigma_species))
for(i in 1:dim(chains$Sigma_species)[1]){
  C[i,,] = chains$Sigma_species[i,,] %*% t(chains$Sigma_species[i,,])
}

j=2
par(mfrow = c(2,3))
for(i in 1:6){
  plot(C[,j,i], type = 'l', main = i)
}











