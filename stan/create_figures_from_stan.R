##----------------------------------------------------
## Name: Joshua North
##
## Date: 07/10/2021
##
## Project: Fish Abundance
##
## Objective: Create figures from the stan output
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
library(scales)
library(progress)
library(sparklyr)
library(stringr)
library(readxl)
library(GGally)
library(xtable)
library(Vizumap)
library(latex2exp)
library(ggcorrplot)
library(corrplot)
library(ggpubr)
library(Matrix)
library(spam)
library(abind)
library(zoo)
library(rstan)
library(ggsci)
# library(pubr)

# load in data ------------------------------------------------------------

static = read_csv('data/Static_MN_Data_JE.csv') %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
effort = read_csv('data/MN_catch_survey.csv') %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
temp = read_rds('data/daily_degree_days_MN_lakes_glm2.rds') %>% ungroup()

temp = temp %>%
  select(date, temp_0, MNDOW_ID, DD5, DOY) %>% # C5 temperature
  rename(SURVEYDATE = date) %>%
  mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>%
  select(-MNDOW_ID)

GDD = temp %>% 
  mutate(year = year(SURVEYDATE)) %>% 
  group_by(year, DOW) %>% 
  summarise(DD5 = max(DD5)) %>% 
  ungroup() %>% 
  group_by(DOW) %>%
  mutate(DD5 = rollmean(DD5, k = 7, align = 'right', fill = NA)) %>% 
  ungroup()


secchi_year = read_csv('data/annual_median_remote_sensing_secchi.csv') %>% 
  mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0"))) %>% 
  select(year, annual.median.rs, DOW) %>% 
  rename('secchi' = 'annual.median.rs')


land = read_xlsx('data/MN_lake_landuse.xlsx')

land = land %>%
  select(MNDOW_ID, year:grass, lake_elevation_m) %>% 
  mutate(DOW = str_split(MNDOW_ID, '_', simplify = T)[,2]) %>% 
  select(-MNDOW_ID) %>% 
  filter(complete.cases(.)) %>% 
  mutate(DOW = as.factor(DOW)) %>% 
  distinct(DOW, .keep_all = T)


# data cleaning -----------------------------------------------------------

# combine survey data with some static variables
fish_dat = effort %>% 
  left_join(static, by = 'DOW') %>%
  select(DOW, COMMON_NAME, TOTAL_CATCH, EFFORT, GEAR, CPUE, SURVEYDATE, MAX_DEPTH_FEET, mean.gdd, LAKE_AREA_GIS_ACRES, 
         LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>%
  mutate(SURVEYDATE = str_split(SURVEYDATE, ' ', simplify = T)[,1]) %>% 
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME),
         GEAR = as.factor(GEAR),
         SURVEYDATE = mdy(SURVEYDATE)) %>%
  filter(complete.cases(.)) %>% 
  filter(COMMON_NAME != 'white sucker',
         COMMON_NAME != 'smallmouth bass',
         COMMON_NAME != 'rock bass') %>%
  filter(SURVEYDATE >= '1993-01-01') %>% 
  mutate(GEAR = as.character(GEAR),
         GEAR = ifelse(GEAR == 'GDE' | GEAR == 'GSH', 'GN', GEAR),
         GEAR = as.factor(GEAR),
         ind = 1) %>% 
  group_by(DOW, COMMON_NAME, SURVEYDATE, GEAR) %>% 
  summarise_all(~sum(.)) %>% 
  ungroup() %>% 
  mutate_at(vars(MAX_DEPTH_FEET:LAKE_CENTER_UTM_NORTHING), ~ ./ind) %>% 
  select(-ind) %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW),
         GEAR = droplevels(GEAR)) %>% 
  arrange(DOW) %>% 
  mutate(year = year(SURVEYDATE))

# add in land use
fish_dat = fish_dat %>% 
  left_join(land %>% select(-year), by = c('DOW')) %>% 
  filter(complete.cases(.))

# all combos of fish, survey date, and gear
all <- fish_dat %>% 
  group_by(DOW) %>% 
  tidyr::expand(COMMON_NAME, SURVEYDATE, GEAR) %>% 
  ungroup() %>% 
  arrange(DOW)


fish_dat <- fish_dat %>% 
  right_join(all) %>%
  mutate(EFFORT = coalesce(EFFORT, 0L),
         TOTAL_CATCH = coalesce(TOTAL_CATCH, 0L),
         CPUE = coalesce(CPUE, 0L)) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:lake_elevation_m), list(~coalesce(., 0L))) %>% 
  group_by(GEAR, SURVEYDATE, DOW) %>% 
  mutate(EFFORT = max(EFFORT)) %>%
  ungroup() %>% 
  group_by(DOW, SURVEYDATE) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:lake_elevation_m), ~ mean(.[!is.na(.) & . != 0])) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:lake_elevation_m), ~ ifelse(is.na(.), 0, .)) %>% 
  ungroup() %>% 
  mutate(TN = ifelse(GEAR == 'TN', 1, 0),
         GN = ifelse(GEAR == 'GN', 1, 0)) %>% 
  select(-GEAR) %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW)) %>%
  mutate(CPUE = TOTAL_CATCH/EFFORT,
         CPUE = ifelse(is.na(CPUE), 0, CPUE)) %>% 
  arrange(DOW)

# add temperature 
# fish_dat = fish_dat %>% 
#   inner_join(GDD) %>% 
#   mutate(DD5 = (DD5 - mean(DD5))/sd(DD5)) %>% 
#   inner_join(temp %>% 
#                select(SURVEYDATE, temp_0, DOW) %>% 
#                mutate(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0))) %>% 
#   inner_join(secchi_year, by = c('DOW', 'year'))

# center by all lakes all time
# fish_dat = fish_dat %>% 
#   inner_join(GDD) %>% 
#   mutate(DD5 = (DD5 - mean(DD5))/sd(DD5)) %>% 
#   inner_join(temp %>% 
#                select(SURVEYDATE, temp_0, DOW)) %>% 
#   inner_join(secchi_year, by = c('DOW', 'year')) %>% 
#   filter(year >= 2000) %>% 
#   mutate(filter_date = ymd(format(SURVEYDATE, "2016-%m-%d"))) %>% 
#   filter(filter_date > ymd('2016-06-01'),
#          filter_date < ymd('2016-09-30')) %>% 
#   select(-filter_date) %>% 
#   mutate(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0)) %>% 
#   mutate(DOY = yday(SURVEYDATE),
#          DOY_sin_semi = sin(DOY/121 * 2*pi),
#          DOY_cos_semi = cos(DOY/121 * 2*pi),
#          DOY_sin_semi_temp = DOY_sin_semi * temp_0,
#          DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>% 
#   mutate(DOW = as.factor(DOW)) %>% 
#   arrange(DOW, year, COMMON_NAME)


# center temp by subset of lakes
center_temp = temp %>% 
  select(SURVEYDATE, temp_0, DOW) %>% 
  mutate(year = year(SURVEYDATE)) %>% 
  filter(year >= 2000) %>% 
  mutate(filter_date = ymd(format(SURVEYDATE, "2016-%m-%d"))) %>% 
  filter(filter_date > ymd('2016-06-01'),
         filter_date < ymd('2016-09-30')) %>% 
  select(-filter_date) %>% 
  group_by(DOW, year) %>% 
  mutate(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0)) %>% 
  ungroup()

# center by lake by year
fish_dat = fish_dat %>% 
  inner_join(GDD) %>% 
  inner_join(secchi_year, by = c('DOW', 'year')) %>% 
  filter(year >= 2000) %>% 
  mutate(filter_date = ymd(format(SURVEYDATE, "2016-%m-%d"))) %>% 
  filter(filter_date > ymd('2016-06-01'),
         filter_date < ymd('2016-09-30')) %>% 
  select(-filter_date) %>% 
  inner_join(center_temp %>% 
               select(SURVEYDATE, temp_0, DOW)) %>% 
  mutate(DD5 = (DD5 - mean(DD5))/sd(DD5)) %>% 
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>% 
  mutate(DOW = as.factor(DOW)) %>% 
  arrange(DOW, year, COMMON_NAME)


# load stan output ---------------------------------------------------------------

# fish_dat = read_csv('data/fish_dat.csv')
# 
# fish_dat = fish_dat %>% 
#   mutate(DOW = as.factor(DOW),
#          COMMON_NAME = as.factor(COMMON_NAME)) %>% 
#   arrange(DOW, year, COMMON_NAME)

fish_dat = fish_dat %>% 
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN)

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
    mutate(TN = 1) %>% 
    relocate(TN) %>%
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

# out = read_rds('/Users/joshuanorth/Desktop/stan_results/stan_lake_random_effect.rds')
# out = read_rds('/Users/joshuanorth/Desktop/stan_lake_random_effect_TN_int.rds')
# out = read_rds('/Users/joshuanorth/Desktop/stan_lake_random_effect_int.rds')
out = read_rds('/Users/joshuanorth/Desktop/stan_results/stan_lake_random_effect.rds')
out_no_catch = read_rds('/Users/joshuanorth/Desktop/stan_results/stan_lake_random_effect_no_catch.rds')

chains = extract(out)
chains_no_catch = extract(out_no_catch)
names(chains)
lapply(chains, dim)

beta_0 = chains$beta_0
beta = chains$beta
phi = chains$phi
omega = chains$omega
sigma_species = chains$Sigma_species
tau = chains$tau

# n_samps = dim(beta)[1]
# for(i in 1:n_samps){
#   
#   a = dat$Z %*% phi[i,,]
#   int_scale = var(c(a))/2 # shift theta star variance to intercept
#   
#   beta_0[i,] = beta_0[i,] + int_scale
#   phi[i,1,] = phi[i,1,] - int_scale
#   
# }



b_names = colnames(dat$X)
phi_names = colnames(dat$Z)

fnames = dat$fish_names
# fnames = dat$fish_names %>% str_replace(., ' ', '_')
# fnames = fnames %>% 
#   str_replace_all(., "_", " ") %>% 
#   str_to_title()

b_names = b_names %>% 
  str_replace_all("_", " ") %>% 
  str_to_title() %>% 
  str_replace_all('Dd5', 'GDD') %>% 
  str_replace_all('Gis', 'GIS')
b_names[1] = "Max Depth"
b_names[2] = "Lake Area"

b_names = c("Int", b_names)

phi_names = phi_names %>% 
  str_replace_all("_", " ") %>% 
  str_to_title() %>% 
  str_replace_all('Gn', 'GN') %>% 
  str_replace_all("Doy ", "") %>% 
  str_replace_all("Semi", "") %>% 
  str_replace_all("0", "") %>% 
  str_replace_all("Sin ", "Sin") %>% 
  str_replace_all("Cos ", "Cos") %>% 
  str_replace_all("Temp  GN", "Temp GN")
phi_names[2] = "Temp"


b0_hat = apply(beta_0, 2, mean)
b_hat = t(apply(beta, c(2,3), mean))
phi_hat = t(apply(phi, c(2,3), mean))


round(cbind(b0_hat, b_hat), 3) %>%
  as_tibble(.name_repair = ~b_names) %>%
  mutate(Species = dat$fish_names) %>%
  relocate(Species)

round(phi_hat, 3) %>%
  as_tibble(.name_repair = ~phi_names) %>%
  mutate(Species = dat$fish_names) %>%
  relocate(Species)



# set up plotting data ----------------------------------------------------

n_samps = dim(beta)[1]
N = dat$N
I = dat$n_lakes
K = dat$K
J = 2
r = ncol(dat$Z)
p = ncol(dat$X)
each_lake = dat$each_lake


theta_post = array(NA, dim = c(n_samps, N, K))
gamma_post = array(NA, dim = c(n_samps, N, K))
gamma_post_no_catch = array(NA, dim = c(n_samps, N, K))
lambda_post = array(NA, dim = c(n_samps, N, K))

for(i in 1:n_samps){
  a = dat$Z %*% phi[i,,]
  # int_scale = var(c(a))/2 # shift theta star variance to intercept
  # a = a - mean(a) - var(c(a))/2 # shift theta star variance to intercept
  # b = matrix(rep(beta_0[i,] + int_scale, N), N, byrow = T) + dat$X %*% beta[i,,] + each_lake %*% omega[i,,]  # shift theta star variance to intercept
  b = matrix(rep(beta_0[i,], N), N, byrow = T) + dat$X %*% beta[i,,] + each_lake %*% omega[i,,]
  
  theta_post[i,,] = exp(a)
  gamma_post[i,,] = exp(b)
  gamma_post_no_catch[i,,] = exp(matrix(rep(chains_no_catch$beta_0[i,], N), N, byrow = T) + dat$X %*% chains_no_catch$beta[i,-8,] + each_lake %*% chains_no_catch$omega[i,,])
  lambda_post[i,,] = exp(a+b)
}



theta = apply(theta_post, c(2,3), mean)
gamma = apply(gamma_post, c(2,3), mean)
gamma_no_catch = apply(gamma_post_no_catch, c(2,3), mean)
lambda = apply(lambda_post, c(2,3), mean)
omega_post = each_lake %*% apply(omega, c(2,3), mean)


fish_plot = rbind(fish_dat %>% filter(COMMON_NAME == 'black crappie') %>% mutate(theta = theta[,1], gamma = gamma[,1], lambda = lambda[,1], omega = omega_post[,1], gamma_no_catch = gamma_no_catch[,1]),
                  fish_dat %>% filter(COMMON_NAME == 'bluegill') %>% mutate(theta = theta[,2], gamma = gamma[,2], lambda = lambda[,2], omega = omega_post[,2], gamma_no_catch = gamma_no_catch[,2]),
                  fish_dat %>% filter(COMMON_NAME == 'largemouth bass') %>% mutate(theta = theta[,3], gamma = gamma[,3], lambda = lambda[,3], omega = omega_post[,3], gamma_no_catch = gamma_no_catch[,3]),
                  fish_dat %>% filter(COMMON_NAME == 'northern pike') %>% mutate(theta = theta[,4], gamma = gamma[,4], lambda = lambda[,4], omega = omega_post[,4], gamma_no_catch = gamma_no_catch[,4]),
                  fish_dat %>% filter(COMMON_NAME == 'walleye') %>% mutate(theta = theta[,5], gamma = gamma[,5], lambda = lambda[,5], omega = omega_post[,5], gamma_no_catch = gamma_no_catch[,5]),
                  fish_dat %>% filter(COMMON_NAME == 'yellow perch') %>% mutate(theta = theta[,6], gamma = gamma[,6], lambda = lambda[,6], omega = omega_post[,6], gamma_no_catch = gamma_no_catch[,6])) %>% 
  mutate(fish = str_replace_all(COMMON_NAME, " ", "_")) %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW') %>% 
  rename(lat = LAKE_CENTER_LAT_DD5,
         lon = LAKE_CENTER_LONG_DD5)


lats = range(fish_plot$lat, na.rm = T)
lons = range(fish_plot$lon, na.rm = T)

usa = st_as_sf(maps::map("state", fill= TRUE, plot = FALSE))




yr_choice = 2016
doy = ymd("2016-08-01")
gdd_center = fish_dat %>% 
  select(year, DOW) %>% 
  inner_join(GDD) %>% 
  summarise(mean = mean(DD5), sd = sd(DD5))

lks = fish_plot %>% 
  group_by(COMMON_NAME) %>% 
  distinct(DOW, .keep_all = T) %>% 
  ungroup() %>% 
  dplyr::select(-c(DD5, secchi, year))
  
fish_plot_year = GDD %>% 
  filter(year == yr_choice) %>% 
  left_join(secchi_year %>% 
              filter(year == yr_choice)) %>% 
  right_join(lks) %>% 
  mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd)

X_year = fish_plot_year %>% 
  distinct(DOW, .keep_all = T) %>% 
  select(MAX_DEPTH_FEET, LAKE_AREA_GIS_ACRES, DD5, secchi, ag, urban, wetlands) %>% 
  select(all_of(mean_covs)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  as.matrix()

nlakes = dim(X_year)[1]
gamma_post = array(NA, dim = c(n_samps, nlakes, K))
gamma_post_no_catch = array(NA, dim = c(n_samps, nlakes, K))

for(i in 1:n_samps){
  
  gamma_post[i,,] = exp(matrix(rep(beta_0[i,], nlakes), nlakes, byrow = T) + X_year %*% beta[i,,] + omega[i,,])
  gamma_post_no_catch[i,,] = exp(matrix(rep(chains_no_catch$beta_0[i,], nlakes), nlakes, byrow = T) + X_year %*% chains_no_catch$beta[i,-8,] + chains_no_catch$omega[i,,])
  
}

gamma_mean = apply(gamma_post, c(2,3), mean)
gamma_no_catch_mean = apply(gamma_post_no_catch, c(2,3), mean)



fish_plot_year = rbind(fish_plot_year %>% filter(COMMON_NAME == 'black crappie') %>% mutate(gamma = gamma_mean[,1], gamma_no_catch = gamma_no_catch_mean[,1]),
                  fish_plot_year %>% filter(COMMON_NAME == 'bluegill') %>% mutate(gamma = gamma_mean[,2], gamma_no_catch = gamma_no_catch_mean[,2]),
                  fish_plot_year %>% filter(COMMON_NAME == 'largemouth bass') %>% mutate(gamma = gamma_mean[,3], gamma_no_catch = gamma_no_catch_mean[,3]),
                  fish_plot_year %>% filter(COMMON_NAME == 'northern pike') %>% mutate(gamma = gamma_mean[,4], gamma_no_catch = gamma_no_catch_mean[,4]),
                  fish_plot_year %>% filter(COMMON_NAME == 'walleye') %>% mutate(gamma = gamma_mean[,5], gamma_no_catch = gamma_no_catch_mean[,5]),
                  fish_plot_year %>% filter(COMMON_NAME == 'yellow perch') %>% mutate(gamma = gamma_mean[,6], gamma_no_catch = gamma_no_catch_mean[,6]))


rm(theta_post, gamma_post, gamma_post_no_catch, lambda_post, theta, gamma, gamma_no_catch, lambda, omega_post, chains, chains_no_catch, temp, gamma_mean, gamma_no_catch_mean)

# xtable output -----------------------------------------------------------
# 
# fnames = fnames %>% 
#   str_replace_all(., "_", " ") %>% 
#   str_to_title()
# 
# b_names = b_names %>% 
#   str_replace_all("_", " ") %>% 
#   str_to_title() %>% 
#   str_replace_all('Dd5', 'DD5') %>% 
#   str_replace_all('Gis', 'GIS')
# 
# b_names = c("Int", b_names)
# 
# phi_names = phi_names %>% 
#   str_replace_all("_", " ") %>% 
#   str_to_title() %>% 
#   str_replace_all('Gn', 'GN')


# relative abundance

b_mean = t(apply(beta, c(2,3), mean))
b_lower = t(apply(beta, c(2,3), quantile, probs = 0.025))
b_upper = t(apply(beta, c(2,3), quantile, probs = 0.975))

b0_mean = apply(beta_0, c(2), mean)
b0_lower = apply(beta_0, c(2), quantile, probs = 0.025)
b0_upper = apply(beta_0, c(2), quantile, probs = 0.975)

b_mean = cbind(b0_mean, b_mean)
b_lower = cbind(b0_lower, b_lower)
b_upper = cbind(b0_upper, b_upper)

b_sig = sign(b_lower) == sign(b_upper)
b_sig = ifelse(b_sig, 1, 0)

colnames(b_mean) = b_names
rownames(b_mean) = fnames

b_mean_sig = round(b_mean, 3)
# cint = round((b_upper - b_lower)/2, 3)
# b_mean_sig = paste0(b_mean_sig, " (", cint, ")")

cint = paste0(" (", round(b_lower, 3), ',', round(b_upper, 3), ")")
b_mean_sig = paste0(b_mean_sig, cint)

b_col = ifelse(b_sig, paste0("\\cellcolor{red}", b_mean_sig), paste0("\\cellcolor{white}", b_mean_sig))
colnames(b_col) = b_names
rownames(b_col) = fnames
print(xtable(t(b_col)), sanitize.text.function = identity)

b_col = ifelse(b_sig, paste0("\\B", round(b_mean, 3)),  round(b_mean, 3))
colnames(b_col) = b_names
rownames(b_col) = fnames
print(xtable(t(b_col)), sanitize.text.function = identity)


# catchability

phi_mean = t(apply(phi, c(2,3), mean))
phi_lower = t(apply(phi, c(2,3), quantile, probs = 0.025))
phi_upper = t(apply(phi, c(2,3), quantile, probs = 0.975))

p_sig = sign(phi_lower) == sign(phi_upper)
p_sig = ifelse(p_sig, 1, 0)

p_mean_sig = round(phi_mean, 3)
# cint = round((phi_upper - phi_lower)/2, 3)
# p_mean_sig = paste0(p_mean_sig, " (", cint, ")")

cint = paste0(" (", round(phi_lower, 3), ',', round(phi_upper, 3), ")")
p_mean_sig = paste0(p_mean_sig, cint)

p_col = ifelse(p_sig, paste0("\\cellcolor{red}", p_mean_sig), paste0("\\cellcolor{white}", p_mean_sig))
colnames(p_col) = phi_names
rownames(p_col) = fnames
print(xtable(t(p_col)), sanitize.text.function = identity)

p_col = ifelse(p_sig, paste0("\\B", round(phi_mean, 3)),  round(phi_mean, 3))
colnames(p_col) = phi_names
rownames(p_col) = fnames
print(xtable(t(p_col)), sanitize.text.function = identity)



# relative abundance w/without scaling

chains_no_catch = extract(out_no_catch)

beta_0_no = chains_no_catch$beta_0
beta_no = chains_no_catch$beta

b_no_mean = t(apply(beta_no, c(2,3), mean))
b_no_lower = t(apply(beta_no, c(2,3), quantile, probs = 0.025))
b_no_upper = t(apply(beta_no, c(2,3), quantile, probs = 0.975))

b0_no_mean = apply(beta_0_no, c(2), mean)
b0_no_lower = apply(beta_0_no, c(2), quantile, probs = 0.025)
b0_no_upper = apply(beta_0_no, c(2), quantile, probs = 0.975)

b_no_mean = cbind(b0_no_mean, b_no_mean)
b_no_lower = cbind(b0_no_lower, b_no_lower)
b_no_upper = cbind(b0_no_upper, b_no_upper)

ind_array = array(NA, dim = dim(b_no_mean))
for(i in 1:nrow(b_no_mean)){
  for(j in 1:ncol(b_no_mean)){
    
    range1 = c(b_lower[i,j], b_upper[i,j])
    range2 = c(b_no_lower[i,j], b_no_upper[i,j])
    
    check1 = range1[1] < range2[1] & range1[2] > range2[1]
    check2 = range2[1] < range1[1] & range2[2] > range1[1]
    
    # check1 = range2[1] > range1[1] & range2[1] < range1[2]
    # check2 = range2[2] > range1[1] & range2[2] < range1[2]
    # check3 = sign(range1[1]) == sign(range2[1])
    # check4 = sign(range1[2]) == sign(range2[2])
    
    # ind_array[i,j] = ifelse((check1 | check2) & check3 & check4, 0, 1)
    ind_array[i,j] = ifelse((check1 | check2), 1, 0)
  }
}

b_no_sig = ifelse(ind_array == 1, 0, 1)

colnames(b_no_mean) = b_names
rownames(b_no_mean) = fnames

b_no_mean_sig = round(b_no_mean, 3)

b_no_col = ifelse(b_no_sig, paste0("\\B", round(b_no_mean, 3)),  round(b_no_mean, 3))
colnames(b_no_col) = b_names
rownames(b_no_col) = fnames
print(xtable(t(b_no_col)), sanitize.text.function = identity)


sign_ind_red = b_mean*b_no_sig > b_no_mean*b_no_sig
sign_ind_blue = b_mean*b_no_sig < b_no_mean*b_no_sig

b_no_color = ifelse(sign_ind_red, paste0("\\textcolor{red}{", b_no_mean_sig, "}"), b_no_mean_sig)
b_no_color = ifelse(sign_ind_blue, paste0("\\textcolor{blue}{", b_no_color, "}"), b_no_color)

b_no_col = ifelse(b_no_sig, paste0("\\B", b_no_color),  b_no_color)
colnames(b_no_col) = b_names
rownames(b_no_col) = fnames
print(xtable(t(b_no_col)), sanitize.text.function = identity)


rm(chains_no_catch)





colnames(b_no_lower) = colnames(b_no_upper) = b_names
rownames(b_no_lower) = rownames(b_no_upper) = fnames

colnames(b_lower) = colnames(b_upper) = b_names
rownames(b_lower) = rownames(b_upper) = fnames

t(b_no_lower)
t(b_no_upper)

t(b_lower)
t(b_upper)


# parameter estimate figures ----------------------------------------------

b_mean = t(apply(beta, c(2,3), mean))
b_lower = t(apply(beta, c(2,3), quantile, probs = 0.025))
b_upper = t(apply(beta, c(2,3), quantile, probs = 0.975))

b0_mean = apply(beta_0, c(2), mean)
b0_lower = apply(beta_0, c(2), quantile, probs = 0.025)
b0_upper = apply(beta_0, c(2), quantile, probs = 0.975)

b_mean = cbind(b0_mean, b_mean)
b_lower = cbind(b0_lower, b_lower)
b_upper = cbind(b0_upper, b_upper)

colnames(b_mean) = b_names
colnames(b_lower) = b_names
colnames(b_upper) = b_names
# rownames(b_mean) = fnames

b_mean = as_tibble(b_mean) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Mean")

b_lower = as_tibble(b_lower) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Lower")

b_upper = as_tibble(b_upper) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Upper")

b_plt = b_mean %>% 
  left_join(b_lower) %>% 
  left_join(b_upper) %>% 
  mutate(sig = as.factor(if_else(sign(Lower) == sign(Upper), 1, 0))) %>% 
  mutate(Variable = factor(Variable, 
                           levels = c('Int', 'Max Depth', 'Lake Area', 'Ag', 'Urban', 'Wetlands', 'GDD', 'Secchi'))) %>% 
  mutate(species = factor(species,
                          levels = c('black crappie', 'bluegill', 'largemouth bass', 'northern pike', 'walleye', 'yellow perch')))


ggplot(b_plt, aes(y = reorder(species, desc(species)))) +
  geom_point(aes(x = Mean, shape = sig, color = sig), size = 1.5) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper, color = sig), width = 0.3) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  facet_wrap(~Variable, scales = 'free_x', ncol = 4) +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual(values = c('darkgray', 'blue')) +
  xlab('Parameter Estimate') +
  ylab('') +
  guides(shape = 'none', color = 'none') 

# ggsave('results/spatial_results/parameter_estimate_plot.png', width = 12, height = 6)


chains_no_catch = extract(out_no_catch)
beta_0_no = chains_no_catch$beta_0
beta_no = chains_no_catch$beta

b_mean = t(apply(beta_no, c(2,3), mean))
b_lower = t(apply(beta_no, c(2,3), quantile, probs = 0.025))
b_upper = t(apply(beta_no, c(2,3), quantile, probs = 0.975))

b0_mean = apply(beta_0_no, c(2), mean)
b0_lower = apply(beta_0_no, c(2), quantile, probs = 0.025)
b0_upper = apply(beta_0_no, c(2), quantile, probs = 0.975)

b_mean = cbind(b0_mean, b_mean)
b_lower = cbind(b0_lower, b_lower)
b_upper = cbind(b0_upper, b_upper)

colnames(b_mean) = b_names
colnames(b_lower) = b_names
colnames(b_upper) = b_names
# rownames(b_mean) = fnames

b_mean_no = as_tibble(b_mean) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Mean")

b_lower_no = as_tibble(b_lower) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Lower")

b_upper_no = as_tibble(b_upper) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Upper")

colnames(b_no_sig) = b_names
b_no_sig = as_tibble(b_no_sig) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "sig")
  

b_plt_no = b_mean_no %>% 
  left_join(b_lower_no) %>% 
  left_join(b_upper_no) %>% 
  left_join(b_no_sig) %>% 
  # mutate(sig = as.factor(if_else(sign(Lower) == sign(Upper), 1, 0))) %>% 
  mutate(Variable = factor(Variable, 
                           levels = c('Int', 'Max Depth', 'Lake Area', 'Ag', 'Urban', 'Wetlands', 'GDD', 'Secchi'))) %>% 
  mutate(species = factor(species,
                          levels = c('black crappie', 'bluegill', 'largemouth bass', 'northern pike', 'walleye', 'yellow perch'))) %>% 
  mutate(sig = factor(sig))


for(i in 1:nrow(b_plt_no)){
  
  if(between(b_plt_no$Lower[i], b_plt$Lower[i], b_plt$Upper[i])){
    b_plt_no$sig[i] = 0
  }else if(between(b_plt_no$Upper[i], b_plt$Lower[i], b_plt$Upper[i])){
    b_plt_no$sig[i] = 0
  }else if(between(b_plt$Lower[i], b_plt_no$Lower[i], b_plt_no$Upper[i])){
    b_plt_no$sig[i] = 0
  }else if(between(b_plt$Upper[i], b_plt_no$Lower[i], b_plt_no$Upper[i])){
    b_plt_no$sig[i] = 0
  }else{
    b_plt_no$sig[i] = 1
  }
  
}


betas_plt = rbind(b_plt %>% mutate(model = 'Scaled Effort'), b_plt_no %>% mutate(model = 'Constant Effort'))

betas_plt %>% 
  filter(species == 'walleye')

ggplot(betas_plt, aes(y = reorder(species, desc(species)))) +
  geom_blank(aes(x = -Mean)) +
  geom_point(aes(x = Mean, color = model, shape = sig), size = 1.5) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper, color = model), width = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  facet_wrap(~Variable, scales = 'free_x', ncol = 4) +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual(values = c('red', 'blue')) +
  xlab('Parameter Estimate') +
  ylab('') +
  guides(shape = 'none') +
  theme(legend.position="bottom",
        legend.title=element_blank())


# ggsave('results/spatial_results/parameter_estimate_plot.png', width = 16, height = 8)





# covariance figure -------------------------------------------------------

cmat = apply(sigma_species, c(2,3), mean)

# fnames = c('Crappie', 'Bluegill', 'Bass', 'Pike', 'Walleye', 'Perch')

cnames = fnames
rnames = fnames

upper_mean = cmat
upper_mean[lower.tri(upper_mean)] = 0
upper_mean <- as_tibble(upper_mean)
colnames(upper_mean) = seq(length(cnames),1)
upper_mean = upper_mean %>% 
  mutate(Row = seq(1,length(cnames))) %>%
  pivot_longer(-Row, names_to = "Col", values_to = "Cov1")

lower_mean = round(cmat, 2)
lower_mean[upper.tri(lower_mean)] = NA
lower_mean <- as_tibble(lower_mean)
colnames(lower_mean) = seq(length(cnames), 1)
lower_mean = lower_mean %>% 
  mutate(Row = seq(1,length(cnames))) %>%
  pivot_longer(-Row, names_to = "Col", values_to = "Cov2")

fnames_df = tibble(name = fnames, ind = seq(1:6))

mean_complete <- upper_mean %>% 
  left_join(lower_mean, by = c('Row', 'Col')) %>% 
  mutate(Cov2 = sprintf( "%0.2f", Cov2)) %>% 
  rowwise() %>% 
  mutate(fish = fnames_df$name[which(Row[1] == fnames_df$ind)]) %>% 
  ungroup() %>% 
  mutate(Cov2 = ifelse(Cov2 == "NA", "", Cov2)) %>% 
  mutate(Cov2 = ifelse(Cov2 == "1.00", fish, Cov2)) %>% 
  select(-fish) %>% 
  mutate(Cov1 = ifelse(Cov1 == 1, NA, Cov1))


ggplot(data = mean_complete) +
  geom_tile(color = "black", aes(Row, Col, fill = Cov1, width=0.95, height=0.95), size = .25) +
  geom_text(aes(Row, Col, label = Cov2), color = "black", size = 8) +
  scale_fill_gradientn(colors = two.colors(n = 29, start = '#053061', end = '#67001f', middle = '#f7f7f7'), limits = c(-1, 1), 
                       guide = guide_colorbar(title = "",
                                              title.position = "bottom",
                                              barwidth = 25,
                                              barheight = 2.5),
                       na.value = 'grey80') +
  labs(x="", y="", title="") +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position="bottom",
        legend.box.margin=margin(-20,0,0,0),
        legend.text=element_text(size=20))

# ggsave('results/spatial_results/species_dependenc.png', width = 20, height = 12)



# with variance on diagonal



fnames_df = tibble(name = paste0(fnames, ": \n", round(apply(tau, 2, mean)^2, 2)), ind = seq(1:6))

mean_complete <- upper_mean %>% 
  left_join(lower_mean, by = c('Row', 'Col')) %>% 
  mutate(Cov2 = sprintf( "%0.2f", Cov2)) %>% 
  rowwise() %>% 
  mutate(fish = fnames_df$name[which(Row[1] == fnames_df$ind)]) %>% 
  ungroup() %>% 
  mutate(Cov2 = ifelse(Cov2 == "NA", "", Cov2)) %>% 
  mutate(Cov2 = ifelse(Cov2 == "1.00", fish, Cov2)) %>% 
  select(-fish) %>% 
  mutate(Cov1 = ifelse(Cov1 == 1, NA, Cov1))


ggplot(data = mean_complete) +
  geom_tile(color = "black", aes(Row, Col, fill = Cov1, width=0.95, height=0.95), size = .25) +
  geom_text(aes(Row, Col, label = Cov2), color = "black", size = 8) +
  scale_fill_gradientn(colors = two.colors(n = 29, start = '#053061', end = '#67001f', middle = '#f7f7f7'), limits = c(-1, 1),
                       guide = guide_colorbar(title = "",
                                              title.position = "bottom",
                                              barwidth = 25,
                                              barheight = 2.5),
                       na.value = 'grey80') +
  labs(x="", y="", title="") +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position="bottom",
        legend.box.margin=margin(-20,0,0,0),
        legend.text=element_text(size=20))


# ggsave('results/spatial_results/species_dependenc_var.png', width = 20, height = 12)


# lake random effect ---------------------------------------------------



ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot %>% 
               group_by(COMMON_NAME) %>% 
               distinct(DOW, .keep_all = T), 
             aes(x = lon, y = lat, color = ifelse(omega < -4, -4, ifelse(omega > 4, 4, omega))), size = 1.5) +
  scale_color_gradient2(low = 'blue', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~fish, 
             labeller = labeller(fish = c('black_crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'northern_pike' = 'northern pike',
                                          'yellow_perch' = 'yellow perch',
                                          'largemouth_bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle("") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave("results/spatial_results/lake_random_effect.png", height = 8, width = 12)



# sample date spatial plot ---------------------------------------------------

date_df = fish_plot %>% 
  mutate(doy = yday(SURVEYDATE)) %>% 
  group_by(DOW) %>% 
  distinct(SURVEYDATE, .keep_all = T) %>% 
  ungroup()

date_df_median = fish_plot %>% 
  mutate(doy = yday(SURVEYDATE)) %>% 
  group_by(DOW) %>% 
  summarise(doymed = median(doy)) %>% 
  left_join(date_df %>% distinct(DOW, .keep_all = T), by = "DOW")

date_df_sd = fish_plot %>% 
  mutate(doy = yday(SURVEYDATE)) %>% 
  group_by(DOW) %>% 
  summarise(doysd = sd(doy)) %>% 
  left_join(date_df %>% distinct(DOW, .keep_all = T), by = "DOW")


ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = date_df_median, 
             aes(x = lon, y = lat, color = doymed), size = 1.5) +
  scale_color_fermenter(breaks = round(unname(quantile(date_df_median$doymed, probs = seq(0.05, 0.95, 0.15))),2), 
                        palette = 'YlOrRd', direction = 1, name = 'Day of Year') +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/median_sample_date_spat.png', width = 12, height = 8)

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = date_df_sd, 
             aes(x = lon, y = lat, color = doysd), size = 1.5) +
  scale_color_fermenter(breaks = round(unname(quantile(date_df_sd$doysd, probs = seq(0.05, 0.95, 0.15))),2), 
                        palette = 'YlOrRd', direction = 1, name = 'Days') +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/sd_sample_date_spat.png', width = 12, height = 8)

# relative abundance ------------------------------------------------------

# quantile(fish_plot_year$gamma, probs = seq(0, 1, 0.01))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot_year, 
             aes(x = lon, y = lat, color = ifelse(gamma > 600, 600, gamma)), size = 1.5) +
  scale_color_gradient(low = 'green', high = 'purple', trans = pseudo_log_trans(sigma = 0.01), breaks = c(0, 1, 10, 100, 600), limits = c(0, 600)) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish, 
             labeller = labeller(fish = c('black_crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'northern_pike' = 'northern pike',
                                          'yellow_perch' = 'yellow perch',
                                          'largemouth_bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle("") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/relative_abundance.png', width = 12, height = 8)



# individual abundance plots
p1 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot %>% filter(fish == 'black_crappie'),
             aes(x = lon, y = lat, color = gamma), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', trans = pseudo_log_trans(sigma = 0.01), breaks = c(0, 1, 10, 50, 200), limits = c(0, 250)) +
  xlab("") +
  ylab("") +
  ggtitle("Black Crappie") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

p2 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot %>% filter(fish == 'bluegill'),
             aes(x = lon, y = lat, color = gamma), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', trans = pseudo_log_trans(sigma = 0.01), breaks = c(0, 1, 10, 100, 300), limits = c(0, 400)) +
  xlab("") +
  ylab("") +
  ggtitle("Bluegill") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

p3 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot %>% filter(fish == 'largemouth_bass'),
             aes(x = lon, y = lat, color = gamma), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', trans = pseudo_log_trans(sigma = 0.01), breaks = c(0, 1, 5, 10), limits = c(0, 15)) +
  xlab("") +
  ylab("") +
  ggtitle("Largemouth Bass") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))


p4 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot %>% filter(fish == 'northern_pike'),
             aes(x = lon, y = lat, color = gamma), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', trans = pseudo_log_trans(sigma = 0.01), breaks = c(0, 1, 2, 5), limits = c(0, 5)) +
  xlab("") +
  ylab("") +
  ggtitle("Northern Pike") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

p5 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot %>% filter(fish == 'walleye'),
             aes(x = lon, y = lat, color = gamma), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', trans = pseudo_log_trans(sigma = 0.01), breaks = c(0, 1, 5, 10), limits = c(0, 10)) +
  xlab("") +
  ylab("") +
  ggtitle("Walleye") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

p6 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot %>% filter(fish == 'yellow_perch'),
             aes(x = lon, y = lat, color = gamma), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', trans = pseudo_log_trans(sigma = 0.01), breaks = c(0, 1, 5, 10, 30), limits = c(0, 32)) +
  xlab("") +
  ylab("") +
  ggtitle("Yellow Perch") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))


rel_abun_grid = cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, align = 'hv')

# ggsave('results/spatial_results/relative_abundance_scale_free.png', rel_abun_grid, width = 15, height = 10)



# relative abundance ratio ------------------------------------------------

cut_pts = c(-100, -80, -70, -60, -50, -40, -30, -20, -10, -1, 1, 10, 20, 30, 40, 50, 60, 70, 80, 100)
cut_pts = c(-100, -25, -20, -15, -10, -5, -1, 1, 5, 10, 15, 20, 25, 100)
labs = c('<-25', '(-25, -20]', '(-20, -15]', '(-15, -10]', '(-10, -5]', '(-5, -1]', '(-1, 1]', '(1, 5]', '(5, 10]', '(10, 15]', '(15, 20]', '(20, 25]', '>25')


fish_ratio = fish_plot_year %>% 
  group_by(COMMON_NAME) %>% 
  mutate(quantile_scale = cut(gamma, breaks = quantile(gamma, probs = seq(0, 1, 0.01)), include.lowest = T, labels = 1:100),
         quantile_no_scale = cut(gamma_no_catch, breaks = quantile(gamma_no_catch, probs = seq(0, 1, 0.01)), include.lowest = T, labels = 1:100)) %>% 
  ungroup() %>%
  mutate(diff = as.numeric(quantile_scale) - as.numeric(quantile_no_scale)) %>% 
  pivot_longer(c(diff), names_to = 'scale', values_to = 'quantile') %>% 
  mutate(brks = cut(quantile, breaks = cut_pts))

labs = levels(cut(fish_ratio$quantile, breaks = cut_pts))

cols = fields::two.colors(n = length(labs), start = 'blue', end = 'red')


ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_jitter(data = fish_ratio,
              aes(x = lon, y = lat, color = brks), size = 1.5, width = 0.05, height = 0.05) +
  scale_color_manual(values = cols) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish,
             labeller = labeller(fish = c('black_crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'northern_pike' = 'northern pike',
                                          'yellow_perch' = 'yellow perch',
                                          'largemouth_bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle("") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 7, fill = cols)))

# ggplot() +
#   geom_sf(data = usa) +
#   coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
#   geom_jitter(data = fish_ratio,
#              aes(x = lon, y = lat, color = quantile), size = 1.5, width = 0.05, height = 0.05) +
#   scale_color_steps2(low = 'blue', high = 'red', breaks = cut_pts, labels = labs, limits = c(-26, 26)) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   facet_wrap(~ fish,
#              labeller = labeller(fish = c('black_crappie' = 'Black Crappie',
#                                           'bluegill' = 'Bluegill',
#                                           'northern_pike' = 'Northern Pike',
#                                           'yellow_perch' = 'Yellow Perch',
#                                           'largemouth_bass' = 'Largemouth Bass',
#                                           'walleye' = 'Walleye'))) +
#   ggtitle("Quantile Comparison for Scaled and Unscaled Effort") +
#   theme(legend.title=element_blank(),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 16, face = 'bold'),
#         title = element_text(size = 16, face = 'bold'),
#         strip.text = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         legend.key.height = unit(2, 'cm'))

# ggsave('results/spatial_results/quantile_comparison_relative_abundance.png', width = 12, height = 8)


fish_ra_comp_lake = fish_plot_year %>% 
  filter(GN == 0) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  mutate(perc_scale = gamma/sum(gamma),
         perc_no_scale = gamma_no_catch/sum(gamma_no_catch)) %>% 
  ungroup() %>% 
  mutate(diff = perc_scale - perc_no_scale) 

ggplot(fish_ra_comp_lake %>% mutate(COMMON_NAME = str_to_title(COMMON_NAME)), aes(y = diff, x = COMMON_NAME)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(c(-0.25, 0.4)) +
  xlab("") +
  ylab("Difference in Species Composition") +
  ggtitle('Species Composition Percentage Difference: Scaled - Unscaled') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18, face = 'bold'),
        title = element_text(size = 18, face = 'bold'),
        legend.position = 'bottom',
        legend.key.size = unit(1, "cm"),
        legend.title=element_blank(),
        strip.text = element_text(size=12))

# ggsave('results/spatial_results/lake_composition_percentages_box.png', width = 12, height = 10)

# quantile(fish_ra_comp_lake$diff, probs = seq(0,1,0.01))
ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_ra_comp_lake,
             aes(x = lon, y = lat, color = ifelse(abs(diff) > 0.3, sign(diff)*0.3, diff)), size = 1.5) +
  scale_color_gradient2(low = 'blue', high = 'red', limits = c(-0.31, 0.31)) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish,
             labeller = labeller(fish = c('black_crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'northern_pike' = 'northern pike',
                                          'yellow_perch' = 'yellow perch',
                                          'largemouth_bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle("") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(2, 'cm'))


# ggsave('results/spatial_results/lake_composition_percentages_map.png', width = 12, height = 8)


# cpue vs relative abundance --------------------------------------------------------

ggplot(fish_plot %>% 
         filter(GN == 1) %>% 
         group_by(fish) %>% 
         mutate_at(vars(CPUE, gamma), ~(. - mean(.))/sd(.)) %>% 
         ungroup(), aes(x = CPUE, y = gamma)) +
  geom_point(size = 3) +
  stat_cor(aes(label = ..r.label..), label.x.npc = 'center', label.y.npc = 'top', size = 10) +
  facet_wrap(~ fish, 
             scales = 'free',
             labeller = labeller(fish = c('black_crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'northern_pike' = 'northern pike',
                                          'yellow_perch' = 'yellow perch',
                                          'largemouth_bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ylab('Relative Abundance (gamma)') +
  xlab('Catch Per Unit Effort') +
  ggtitle('') + 
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/cpue_rel_GN_standardized.png', width = 12, height = 10)


ggplot(fish_plot %>% 
         filter(GN == 0) %>% 
         group_by(fish) %>% 
         mutate_at(vars(CPUE, gamma), ~(. - mean(.))/sd(.)) %>% 
         ungroup(), aes(x = CPUE, y = gamma)) +
  geom_point(size = 3) +
  stat_cor(aes(label = ..r.label..), label.x.npc = 'center', label.y.npc = 'top', size = 10) +
  facet_wrap(~ fish, 
             scales = 'free',
             labeller = labeller(fish = c('black_crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'northern_pike' = 'northern pike',
                                          'yellow_perch' = 'yellow perch',
                                          'largemouth_bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ylab('Relative Abundance (gamma)') +
  xlab('Catch Per Unit Effort') +
  ggtitle('') + 
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/cpue_rel_TN_standardized.png', width = 12, height = 10)


#### no axis

ggplot(fish_plot %>% filter(GN == 1), aes(x = CPUE, y = gamma)) +
  geom_point(size = 3) +
  stat_cor(aes(label = ..r.label..), label.x.npc = 'center', label.y.npc = 'top', size = 10) +
  facet_wrap(~ fish, 
             scales = 'free',
             labeller = labeller(fish = c('black_crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'northern_pike' = 'northern pike',
                                          'yellow_perch' = 'yellow perch',
                                          'largemouth_bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ylab('Relative Abundance (gamma)') +
  xlab('Catch Per Unit Effort') +
  ggtitle('') + 
  theme(legend.title=element_blank(),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'),
        axis.text=element_blank(),
        axis.ticks=element_blank())

# ggsave('results/spatial_results/cpue_rel_GN_no_axis.png', width = 12, height = 10)


ggplot(fish_plot %>% filter(TN == 1), aes(x = CPUE, y = gamma)) +
  geom_point(size = 3) +
  stat_cor(aes(label = ..r.label..), label.x.npc = 'center', label.y.npc = 'top', size = 10) +
  facet_wrap(~ fish, 
             scales = 'free',
             labeller = labeller(fish = c('black_crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'northern_pike' = 'northern pike',
                                          'yellow_perch' = 'yellow perch',
                                          'largemouth_bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ylab('Relative Abundance (gamma)') +
  xlab('Catch Per Unit Effort') +
  ggtitle('') + 
  theme(legend.title=element_blank(),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'),
        axis.text=element_blank(),
        axis.ticks=element_blank())

# ggsave('results/spatial_results/cpue_rel_TN_no_axis.png', width = 12, height = 10)




# spatial plot of TN vs GN ------------------------------------------------

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot %>% 
               filter(GN == 1), 
             aes(x = lon, y = lat, color = log(lambda)), size = 1.5) +
  scale_color_gradient2(low = 'blue', high = 'red', limits = c(-9, 7)) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish, 
             labeller = labeller(fish = c('black_crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'northern_pike' = 'northern pike',
                                          'yellow_perch' = 'yellow perch',
                                          'largemouth_bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle("Log theta*gamma for GN") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/gear_GN_equivalence.png', width = 12, height = 12)


ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot %>% 
               filter(GN == 0), 
             aes(x = lon, y = lat, color = log(lambda)), size = 1.5) +
  scale_color_gradient2(low = 'blue', high = 'red', limits = c(-9, 7)) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish, 
             labeller = labeller(fish = c('black_crappie' = 'Black Crappie',
                                          'bluegill' = 'Bluegill',
                                          'northern_pike' = 'Northern Pike',
                                          'yellow_perch' = 'Yellow Perch',
                                          'largemouth_bass' = 'Largemouth Bass',
                                          'walleye' = 'Walleye'))) +
  ggtitle("Log theta*gamma for TN") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/gear_TN_equivalence.png', width = 12, height = 12)



# spatial plot of TN vs GN for a given day --------------------------------

day_choice = ymd('2016-08-01')

gdd_center = fish_dat %>% 
  select(year, DOW) %>% 
  inner_join(GDD) %>% 
  summarise(mean = mean(DD5), sd = sd(DD5))

day_covs_TN = fish_plot %>%
  select(DOW, all_of(mean_covs)) %>% 
  group_by(DOW) %>% 
  summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
  select(-c(DD5, secchi)) %>% 
  left_join(secchi_year %>% 
              filter(year == year(day_choice)) %>% select(-year)) %>% 
  left_join(GDD %>% 
              filter(year == year(day_choice)) %>% select(-year)) %>% 
  relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
  relocate(secchi, .after = DD5)  %>% 
  mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
  left_join(center_temp %>%
              filter(SURVEYDATE == day_choice), by = 'DOW') %>% 
  select(-c(year)) %>% 
  mutate(TN = 1) %>% 
  relocate(temp_0, .after = TN) %>% 
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0,
         GN = 0) %>% 
  select(-DOY) %>% 
  mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  select(-c(SURVEYDATE)) %>% 
  arrange(DOW)

day_covs_GN = fish_plot %>%
  select(DOW, all_of(mean_covs)) %>% 
  group_by(DOW) %>% 
  summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
  select(-c(DD5, secchi)) %>% 
  left_join(secchi_year %>% 
              filter(year == year(day_choice)) %>% select(-year)) %>% 
  left_join(GDD %>% 
              filter(year == year(day_choice)) %>% select(-year)) %>% 
  relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
  relocate(secchi, .after = DD5) %>% 
  mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
  left_join(center_temp %>%
              filter(SURVEYDATE == day_choice), by = 'DOW') %>% 
  select(-c(year)) %>% 
  mutate(TN = 1) %>% 
  relocate(temp_0, .after = TN) %>% 
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0,
         GN = 1) %>% 
  select(-DOY) %>% 
  mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  select(-c(SURVEYDATE)) %>% 
  arrange(DOW)


day_TN = day_covs_TN %>% 
  select(-DOW) %>% 
  as.matrix()

day_GN = day_covs_GN %>% 
  select(-DOW) %>% 
  as.matrix()

n_lakes = dat$n_lakes

day_TN_post = day_GN_post = array(NA, dim = c(n_samps, n_lakes, K))


for(i in 1:n_samps){
  coefs = rbind(beta[i,,], phi[i,,])
  
  day_TN_post[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_TN %*% coefs + omega[i,,])
  day_GN_post[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_GN %*% coefs + omega[i,,])
  
}

TN_post = apply(day_TN_post, c(2,3), mean)
GN_post = apply(day_GN_post, c(2,3), mean)

fish_names = c("crappie", 'bluegill', 'bass', 'pike', 'walleye', 'perch')

colnames(TN_post) = paste0('TN_', fish_names)
colnames(GN_post) = paste0('GN_', fish_names)

lake_loc = fish_plot %>% 
  distinct(DOW, .keep_all = T) %>% 
  select(DOW, lat, lon)

day_comp = day_covs_TN %>% 
  select(DOW) %>% 
  cbind(TN_post) %>% 
  cbind(GN_post) %>% 
  tibble() %>% 
  left_join(lake_loc) %>% 
  pivot_longer(-c(DOW, lat, lon), names_to = c('Gear', 'fish'), names_sep = '_', values_to = 'value')


# predicted catch on 2016-08-01 for TN
ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = day_comp %>%
               filter(Gear == 'TN'),
             aes(x = lon, y = lat, color = ifelse(value > 300, 300, value)), size = 1.5) +
  scale_color_gradient(low = 'green', high = 'purple', trans = pseudo_log_trans(sigma = 0.01), breaks = c(0, 1, 10, 100, 300), limits = c(0, 300)) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish,
             labeller = labeller(fish = c('crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'pike' = 'northern pike',
                                          'perch' = 'yellow perch',
                                          'bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle('') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/gear_TN_equivalence.png', width = 12, height = 8)


# predicted catch on 2016-08-01 for GN
ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = day_comp %>%
               filter(Gear == 'GN'),
             aes(x = lon, y = lat, color = ifelse(value > 300, 300, value)), size = 1.5) +
  scale_color_gradient(low = 'green', high = 'purple', trans = pseudo_log_trans(sigma = 0.01), breaks = c(0, 1, 10, 100, 300), limits = c(0, 300)) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish,
             labeller = labeller(fish = c('crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'pike' = 'northern pike',
                                          'perch' = 'yellow perch',
                                          'bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle('') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/gear_GN_equivalence.png', width = 12, height = 8)

crappie_diff = day_comp %>% 
  filter(fish == 'crappie') %>% 
  group_by(DOW) %>% 
  mutate(diff = value - lag(value)) %>% 
  ungroup() %>% 
  filter(complete.cases(.)) %>% 
  select(-value)

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = crappie_diff,
             aes(x = lon, y = lat, color = diff), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', trans = pseudo_log_trans(sigma = 0.01), breaks = c(0, 1, 10, 50), limits = c(0, 60)) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle(paste0('Difference in Crappie Catch for 1 Unit Effort on ', day_choice, ': GN - TN')) +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/gear_equivalence_crappie_diff.png', width = 12, height = 12)


# spatial plot of TN vs GN for two days --------------------------------

day_choice1 = ymd('2016-06-15')
day_choice2 = ymd('2016-08-15')

gdd_center = fish_dat %>% 
  select(year, DOW) %>% 
  inner_join(GDD) %>% 
  summarise(mean = mean(DD5), sd = sd(DD5))

day_covs_TN1 = fish_plot %>%
  select(DOW, all_of(mean_covs)) %>% 
  group_by(DOW) %>% 
  summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
  select(-c(DD5, secchi)) %>% 
  left_join(secchi_year %>% 
              filter(year == year(day_choice1)) %>% select(-year)) %>% 
  left_join(GDD %>% 
              filter(year == year(day_choice1)) %>% select(-year)) %>% 
  relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
  relocate(secchi, .after = DD5)  %>% 
  mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
  left_join(center_temp %>%
              filter(SURVEYDATE == day_choice1), by = 'DOW') %>% 
  select(-c(year)) %>% 
  mutate(TN = 1) %>% 
  relocate(temp_0, .after = TN) %>% 
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0,
         GN = 0) %>% 
  select(-DOY) %>% 
  mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  select(-c(SURVEYDATE)) %>% 
  arrange(DOW)

day_covs_TN2 = fish_plot %>%
  select(DOW, all_of(mean_covs)) %>% 
  group_by(DOW) %>% 
  summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
  select(-c(DD5, secchi)) %>% 
  left_join(secchi_year %>% 
              filter(year == year(day_choice2)) %>% select(-year)) %>% 
  left_join(GDD %>% 
              filter(year == year(day_choice2)) %>% select(-year)) %>% 
  relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
  relocate(secchi, .after = DD5)  %>% 
  mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
  left_join(center_temp %>%
              filter(SURVEYDATE == day_choice2), by = 'DOW') %>% 
  select(-c(year)) %>% 
  mutate(TN = 1) %>% 
  relocate(temp_0, .after = TN) %>% 
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0,
         GN = 0) %>% 
  select(-DOY) %>% 
  mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  select(-c(SURVEYDATE)) %>% 
  arrange(DOW)

day_covs_GN1 = fish_plot %>%
  select(DOW, all_of(mean_covs)) %>% 
  group_by(DOW) %>% 
  summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
  select(-c(DD5, secchi)) %>% 
  left_join(secchi_year %>% 
              filter(year == year(day_choice1)) %>% select(-year)) %>% 
  left_join(GDD %>% 
              filter(year == year(day_choice1)) %>% select(-year)) %>% 
  relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
  relocate(secchi, .after = DD5) %>% 
  mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
  left_join(center_temp %>%
              filter(SURVEYDATE == day_choice1), by = 'DOW') %>% 
  select(-c(year)) %>% 
  mutate(TN = 1) %>% 
  relocate(temp_0, .after = TN) %>% 
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0,
         GN = 1) %>% 
  select(-DOY) %>% 
  mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  select(-c(SURVEYDATE)) %>% 
  arrange(DOW)

day_covs_GN2 = fish_plot %>%
  select(DOW, all_of(mean_covs)) %>% 
  group_by(DOW) %>% 
  summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
  select(-c(DD5, secchi)) %>% 
  left_join(secchi_year %>% 
              filter(year == year(day_choice2)) %>% select(-year)) %>% 
  left_join(GDD %>% 
              filter(year == year(day_choice2)) %>% select(-year)) %>% 
  relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
  relocate(secchi, .after = DD5) %>% 
  mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
  left_join(center_temp %>%
              filter(SURVEYDATE == day_choice2), by = 'DOW') %>% 
  select(-c(year)) %>% 
  mutate(TN = 1) %>% 
  relocate(temp_0, .after = TN) %>% 
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0,
         GN = 1) %>% 
  select(-DOY) %>% 
  mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  select(-c(SURVEYDATE)) %>% 
  arrange(DOW)


day_TN1 = day_covs_TN1 %>% 
  select(-DOW) %>% 
  as.matrix()

day_TN2 = day_covs_TN2 %>% 
  select(-DOW) %>% 
  as.matrix()

day_GN1 = day_covs_GN1 %>% 
  select(-DOW) %>% 
  as.matrix()

day_GN2 = day_covs_GN2 %>% 
  select(-DOW) %>% 
  as.matrix()

n_lakes = dat$n_lakes

day_TN_post1 = day_GN_post1 = day_TN_post2 = day_GN_post2 = array(NA, dim = c(n_samps, n_lakes, K))


for(i in 1:n_samps){
  coefs = rbind(beta[i,,], phi[i,,])
  
  day_TN_post1[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_TN1 %*% coefs + omega[i,,])
  day_GN_post1[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_GN1 %*% coefs + omega[i,,])
  day_TN_post2[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_TN2 %*% coefs + omega[i,,])
  day_GN_post2[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_GN2 %*% coefs + omega[i,,])
  
}

TN_post1 = apply(day_TN_post1, c(2,3), mean)
GN_post1 = apply(day_GN_post1, c(2,3), mean)
TN_post2 = apply(day_TN_post2, c(2,3), mean)
GN_post2 = apply(day_GN_post2, c(2,3), mean)

fish_names = c("crappie", 'bluegill', 'bass', 'pike', 'walleye', 'perch')

colnames(TN_post1) = paste0('TN_', fish_names)
colnames(GN_post1) = paste0('GN_', fish_names)
colnames(TN_post2) = paste0('TN_', fish_names)
colnames(GN_post2) = paste0('GN_', fish_names)

lake_loc = fish_plot %>% 
  distinct(DOW, .keep_all = T) %>% 
  select(DOW, lat, lon)

day_comp1 = day_covs_TN1 %>% 
  select(DOW) %>% 
  cbind(TN_post1) %>% 
  cbind(GN_post1) %>% 
  tibble() %>% 
  left_join(lake_loc) %>% 
  pivot_longer(-c(DOW, lat, lon), names_to = c('Gear', 'fish'), names_sep = '_', values_to = 'value1')

day_comp2 = day_covs_TN2 %>% 
  select(DOW) %>% 
  cbind(TN_post2) %>% 
  cbind(GN_post2) %>% 
  tibble() %>% 
  left_join(lake_loc) %>% 
  pivot_longer(-c(DOW, lat, lon), names_to = c('Gear', 'fish'), names_sep = '_', values_to = 'value2')

day_comp = day_comp1 %>% 
  mutate(doy1 = day_choice1) %>% 
  left_join(day_comp2 %>% 
              mutate(doy2 = day_choice2)) %>% 
  mutate(diff = value2 - value1) %>% 
  mutate(diff = if_else(abs(diff) > 70, 70*sign(diff), diff))

# quantile(day_comp$diff, probs = seq(0, 1, 0.01))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = day_comp %>%
               filter(Gear == 'TN'),
             aes(x = lon, y = lat, color = diff), size = 1.5) +
  scale_color_gradient2(low = 'blue', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish,
             labeller = labeller(fish = c('crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'pike' = 'northern pike',
                                          'perch' = 'yellow perch',
                                          'bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle(paste0('Difference in Predicted Catch for 1 Unit Effort on - TN')) +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

ggsave('results/spatial_results/gear_TN_difference.png', width = 12, height = 12)


ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = day_comp %>%
               filter(Gear == 'GN'),
             aes(x = lon, y = lat, color = diff), size = 1.5) +
  scale_color_gradient2(low = 'blue', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish,
             labeller = labeller(fish = c('crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'pike' = 'northern pike',
                                          'perch' = 'yellow perch',
                                          'bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle(paste0('Difference in Predicted Catch for 1 Unit Effort on - GN')) +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/gear_GN_difference.png', width = 12, height = 12)


# spatial plot of TN vs GN diff two days --------------------------------

day1 = ymd('2016-06-15')
day2 = ymd('2016-08-15')
yr = 2016
phis = phi
betas = beta
beta0s = beta_0

sample_day_diff = function(yr, phis, betas, beta0s, day1, day2){
  
  gdd_center = fish_dat %>% 
    select(year, DOW) %>% 
    inner_join(GDD) %>% 
    summarise(mean = mean(DD5), sd = sd(DD5))
  
  dates = center_temp %>% 
    filter(year == yr) %>% 
    distinct(SURVEYDATE) %>% 
    pull()
  
  n_lakes = dat$n_lakes
  n_days = length(dates)
  
  day_TN_mean = day_GN_mean = array(NA, dim = c(n_days, K))
  # day_TN_post = day_GN_post = array(NA, dim = c(n_days, n_samps, n_lakes, K))
  
  for(j in seq_along(dates)){
    
    day_TN_post = day_GN_post = array(NA, dim = c(n_samps, n_lakes, K))
    
    day_choice = dates[j]
    
    day_covs_TN = fish_plot %>%
      select(DOW, all_of(mean_covs)) %>% 
      group_by(DOW) %>% 
      summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
      select(-c(DD5, secchi)) %>% 
      left_join(secchi_year %>% 
                  filter(year == year(day_choice)) %>% select(-year)) %>% 
      left_join(GDD %>% 
                  filter(year == year(day_choice)) %>% select(-year)) %>% 
      relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
      relocate(secchi, .after = DD5)  %>% 
      mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
      left_join(center_temp %>%
                  filter(SURVEYDATE == day_choice), by = 'DOW') %>% 
      select(-c(year)) %>% 
      mutate(TN = 1) %>% 
      relocate(temp_0, .after = TN) %>% 
      mutate(DOY = yday(SURVEYDATE),
             DOY_sin_semi = sin(DOY/121 * 2*pi),
             DOY_cos_semi = cos(DOY/121 * 2*pi),
             DOY_sin_semi_temp = DOY_sin_semi * temp_0,
             DOY_cos_semi_temp = DOY_cos_semi * temp_0,
             GN = 0) %>% 
      select(-DOY) %>% 
      mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
      mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
      mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
      mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
      mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
      select(-c(SURVEYDATE)) %>% 
      arrange(DOW)
    
    day_covs_GN = fish_plot %>%
      select(DOW, all_of(mean_covs)) %>% 
      group_by(DOW) %>% 
      summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
      select(-c(DD5, secchi)) %>% 
      left_join(secchi_year %>% 
                  filter(year == year(day_choice)) %>% select(-year)) %>% 
      left_join(GDD %>% 
                  filter(year == year(day_choice)) %>% select(-year)) %>% 
      relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
      relocate(secchi, .after = DD5) %>% 
      mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
      left_join(center_temp %>%
                  filter(SURVEYDATE == day_choice), by = 'DOW') %>% 
      select(-c(year)) %>% 
      mutate(TN = 1) %>% 
      relocate(temp_0, .after = TN) %>% 
      mutate(DOY = yday(SURVEYDATE),
             DOY_sin_semi = sin(DOY/121 * 2*pi),
             DOY_cos_semi = cos(DOY/121 * 2*pi),
             DOY_sin_semi_temp = DOY_sin_semi * temp_0,
             DOY_cos_semi_temp = DOY_cos_semi * temp_0,
             GN = 1) %>% 
      select(-DOY) %>% 
      mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
      mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
      mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
      mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
      mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
      select(-c(SURVEYDATE)) %>% 
      arrange(DOW)
    
    
    day_TN = day_covs_TN %>% 
      select(-DOW) %>% 
      as.matrix()
    
    day_GN = day_covs_GN %>% 
      select(-DOW) %>% 
      as.matrix()
    
    for(i in 1:n_samps){
      coefs = rbind(beta[i,,], phi[i,,])
      
      day_TN_post[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_TN %*% coefs + omega[i,,])
      day_GN_post[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_GN %*% coefs + omega[i,,])
      
    }
    
    # day_TN_mean[j] = mean(day_TN_post)
    # day_GN_mean[j] = mean(day_GN_post)
    
    day_TN_mean[j,] = apply(day_TN_post, 3, mean)
    day_GN_mean[j,] = apply(day_GN_post, 3, mean)
    
  }
  
  # n_samps, n_lakes, K
  # nscale = n_samps * n_lakes * K
  # TN_mean = sum(day_TN_mean * nscale)/(n_days*nscale)
  # GN_mean = sum(day_GN_mean * nscale)/(n_days*nscale)
  
  nscale = n_samps * n_lakes
  TN_mean = apply(day_TN_mean, 2, function(x) sum(x * nscale)/(n_days*nscale))
  GN_mean = apply(day_GN_mean, 2, function(x) sum(x * nscale)/(n_days*nscale))
  
  
  day_covs_TN1 = fish_plot %>%
    select(DOW, all_of(mean_covs)) %>% 
    group_by(DOW) %>% 
    summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
    select(-c(DD5, secchi)) %>% 
    left_join(secchi_year %>% 
                filter(year == year(day1)) %>% select(-year)) %>% 
    left_join(GDD %>% 
                filter(year == year(day1)) %>% select(-year)) %>% 
    relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
    relocate(secchi, .after = DD5)  %>% 
    mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
    left_join(center_temp %>%
                filter(SURVEYDATE == day1), by = 'DOW') %>% 
    select(-c(year)) %>% 
    mutate(TN = 1) %>% 
    relocate(temp_0, .after = TN) %>% 
    mutate(DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/121 * 2*pi),
           DOY_cos_semi = cos(DOY/121 * 2*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0,
           GN = 0) %>% 
    select(-DOY) %>% 
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
    select(-c(SURVEYDATE)) %>% 
    arrange(DOW)
  
  day_covs_TN2 = fish_plot %>%
    select(DOW, all_of(mean_covs)) %>% 
    group_by(DOW) %>% 
    summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
    select(-c(DD5, secchi)) %>% 
    left_join(secchi_year %>% 
                filter(year == year(day2)) %>% select(-year)) %>% 
    left_join(GDD %>% 
                filter(year == year(day2)) %>% select(-year)) %>% 
    relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
    relocate(secchi, .after = DD5)  %>% 
    mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
    left_join(center_temp %>%
                filter(SURVEYDATE == day2), by = 'DOW') %>% 
    select(-c(year)) %>% 
    mutate(TN = 1) %>% 
    relocate(temp_0, .after = TN) %>% 
    mutate(DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/121 * 2*pi),
           DOY_cos_semi = cos(DOY/121 * 2*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0,
           GN = 0) %>% 
    select(-DOY) %>% 
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
    select(-c(SURVEYDATE)) %>% 
    arrange(DOW)
  
  day_covs_GN1 = fish_plot %>%
    select(DOW, all_of(mean_covs)) %>% 
    group_by(DOW) %>% 
    summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
    select(-c(DD5, secchi)) %>% 
    left_join(secchi_year %>% 
                filter(year == year(day_choice1)) %>% select(-year)) %>% 
    left_join(GDD %>% 
                filter(year == year(day_choice1)) %>% select(-year)) %>% 
    relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
    relocate(secchi, .after = DD5) %>% 
    mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
    left_join(center_temp %>%
                filter(SURVEYDATE == day_choice1), by = 'DOW') %>% 
    select(-c(year)) %>% 
    mutate(TN = 1) %>% 
    relocate(temp_0, .after = TN) %>% 
    mutate(DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/121 * 2*pi),
           DOY_cos_semi = cos(DOY/121 * 2*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0,
           GN = 1) %>% 
    select(-DOY) %>% 
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
    select(-c(SURVEYDATE)) %>% 
    arrange(DOW)
  
  day_covs_GN2 = fish_plot %>%
    select(DOW, all_of(mean_covs)) %>% 
    group_by(DOW) %>% 
    summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
    select(-c(DD5, secchi)) %>% 
    left_join(secchi_year %>% 
                filter(year == year(day_choice2)) %>% select(-year)) %>% 
    left_join(GDD %>% 
                filter(year == year(day_choice2)) %>% select(-year)) %>% 
    relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
    relocate(secchi, .after = DD5) %>% 
    mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd) %>% 
    left_join(center_temp %>%
                filter(SURVEYDATE == day_choice2), by = 'DOW') %>% 
    select(-c(year)) %>% 
    mutate(TN = 1) %>% 
    relocate(temp_0, .after = TN) %>% 
    mutate(DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/121 * 2*pi),
           DOY_cos_semi = cos(DOY/121 * 2*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0,
           GN = 1) %>% 
    select(-DOY) %>% 
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
    select(-c(SURVEYDATE)) %>% 
    arrange(DOW)
  
  
  day_TN1 = day_covs_TN1 %>% 
    select(-DOW) %>% 
    as.matrix()
  
  day_TN2 = day_covs_TN2 %>% 
    select(-DOW) %>% 
    as.matrix()
  
  day_GN1 = day_covs_GN1 %>% 
    select(-DOW) %>% 
    as.matrix()
  
  day_GN2 = day_covs_GN2 %>% 
    select(-DOW) %>% 
    as.matrix()
  
  n_lakes = dat$n_lakes
  
  day_TN_post1 = day_GN_post1 = day_TN_post2 = day_GN_post2 = array(NA, dim = c(n_samps, n_lakes, K))
  
  
  for(i in 1:n_samps){
    coefs = rbind(beta[i,,], phi[i,,])
    
    day_TN_post1[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_TN1 %*% coefs + omega[i,,])
    day_GN_post1[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_GN1 %*% coefs + omega[i,,])
    day_TN_post2[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_TN2 %*% coefs + omega[i,,])
    day_GN_post2[i,,] = exp(matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + day_GN2 %*% coefs + omega[i,,])
    
  }
  
  TN_post1 = apply(day_TN_post1, c(2,3), mean)
  GN_post1 = apply(day_GN_post1, c(2,3), mean)
  TN_post2 = apply(day_TN_post2, c(2,3), mean)
  GN_post2 = apply(day_GN_post2, c(2,3), mean)
  
  fish_names = c("crappie", 'bluegill', 'bass', 'pike', 'walleye', 'perch')
  
  colnames(TN_post1) = paste0('TN_', fish_names)
  colnames(GN_post1) = paste0('GN_', fish_names)
  colnames(TN_post2) = paste0('TN_', fish_names)
  colnames(GN_post2) = paste0('GN_', fish_names)
  
  lake_loc = fish_plot %>% 
    distinct(DOW, .keep_all = T) %>% 
    select(DOW, lat, lon)
  
  day_comp1 = day_covs_TN1 %>% 
    select(DOW) %>% 
    cbind(TN_post1) %>% 
    cbind(GN_post1) %>% 
    tibble() %>% 
    left_join(lake_loc) %>% 
    pivot_longer(-c(DOW, lat, lon), names_to = c('Gear', 'fish'), names_sep = '_', values_to = 'value1')
  
  day_comp2 = day_covs_TN2 %>% 
    select(DOW) %>% 
    cbind(TN_post2) %>% 
    cbind(GN_post2) %>% 
    tibble() %>% 
    left_join(lake_loc) %>% 
    pivot_longer(-c(DOW, lat, lon), names_to = c('Gear', 'fish'), names_sep = '_', values_to = 'value2')
  
  # day_comp = day_comp1 %>% 
  #   mutate(doy1 = day_choice1) %>% 
  #   left_join(day_comp2 %>% 
  #               mutate(doy2 = day_choice2)) %>% 
  #   mutate(diff = value2 - value1) %>% 
  #   mutate(diff = if_else(Gear == "TN", diff/TN_mean, diff/GN_mean)) %>% 
  #   select(-c(value1:doy2))
  
  day_comp = day_comp1 %>% 
    mutate(doy1 = day_choice1) %>% 
    left_join(day_comp2 %>% 
                mutate(doy2 = day_choice2)) %>% 
    mutate(diff = value2 - value1) %>% 
    group_by(DOW) %>% 
    mutate(diff = if_else(Gear == "TN", diff/TN_mean, diff/GN_mean)) %>% 
    ungroup() %>% 
    select(-c(value1:doy2))
  
  TN_mean
  GN_mean
  
  return(day_comp)
  
}

day_comp = sample_day_diff(yr, phis, betas, beta0s, day1, day2)

day_comp = day_comp %>% 
  mutate(fish = factor(fish, levels = c('crappie','bluegill', 'bass','pike','walleye', 'perch')))

quantile(day_comp$diff, probs = seq(0, 1, 0.01))

# Standardized Difference in Predicted Catch for 1 Unit Effort - TN
ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = day_comp %>%
               filter(Gear == 'TN'),
             aes(x = lon, y = lat, color = if_else(abs(diff) > 2, sign(diff)*2, diff)), size = 1.5) +
  scale_color_gradient2(low = 'blue', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish,
             labeller = labeller(fish = c('crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'pike' = 'northern pike',
                                          'perch' = 'yellow perch',
                                          'bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle('') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/gear_TN_difference.png', width = 12, height = 8)


# Standardized Difference in Predicted Catch for 1 Unit Effort - GN
ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = day_comp %>%
               filter(Gear == 'GN'),
             aes(x = lon, y = lat, color = if_else(abs(diff) > 2, sign(diff)*2, diff)), size = 1.5) +
  scale_color_gradient2(low = 'blue', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish,
             labeller = labeller(fish = c('crappie' = 'black crappie',
                                          'bluegill' = 'bluegill',
                                          'pike' = 'northern pike',
                                          'perch' = 'yellow perch',
                                          'bass' = 'largemouth bass',
                                          'walleye' = 'walleye'))) +
  ggtitle('') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/gear_GN_difference.png', width = 12, height = 8)



# catchability curve boxplots ------------------------------------------------------
# choose year!!!!
# post curves for all lakes for a year

phis = phi
fish_names = c("crappie", 'bluegill', 'bass', 'pike', 'walleye', 'perch')
day_choice = ymd('2016-08-01')
df = fish_plot

# center_temp
catchability_day = function(phis, fish_names, day_choice, df){
  
  tC = center_temp %>%
    filter(SURVEYDATE == day_choice)
  
  lakes = df %>% 
    distinct(DOW)
  
  n_samps = dim(phis)[1]
  n_lakes = nrow(lakes)
  catch_post_GN = catch_post_TN = array(NA, dim = c(6, n_lakes, n_samps))
  
  for(i in 1:n_lakes){
    
    lake = c(lakes[i,])$DOW
    tC_tmp = tC %>% 
      filter(DOW == lake)
    
    TN_df = tC_tmp %>% 
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
    
    GN_df = tC_tmp %>% 
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
    
    TN = as.matrix(TN_df)
    GN = as.matrix(GN_df)
    
    for(j in 1:n_samps){
      catch_post_GN[,i,j] = GN %*% phis[j,,]
      catch_post_TN[,i,j] = TN %*% phis[j,,]
    }
    
  }
  
  TN = t(apply(catch_post_TN, c(1,2), mean))
  GN = t(apply(catch_post_GN, c(1,2), mean))
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  return(tibble(DOW = lakes$DOW, as_tibble(TN), as_tibble(GN)))
}

post_curve = catchability_day(phis, fish_names, day_choice, df)

catch_plot = post_curve %>%
  pivot_longer(cols = -c(DOW), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "value") %>% 
  mutate(Fish = rep(fnames, n()/6))

ggplot(catch_plot, aes(x = Fish, y = exp(value))) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  facet_wrap(~Gear, scales = 'free_y') +
  xlab('') +
  ylab('Catchability') +
  ggtitle(paste0('Gear Type Scaling by Species - ', day_choice)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18, face = 'bold'),
        title = element_text(size = 18, face = 'bold'),
        legend.position = 'bottom',
        legend.key.size = unit(1, "cm"),
        legend.title=element_blank(),
        strip.text = element_text(size=12))

# ggsave(paste0('results/spatial_results/catchability_box_plot_exp.png'), width = 12, height = 8)
# ggsave(paste0('results/spatial_results/catchability_box_plot_log.png'), width = 12, height = 8)



# catchability curve ------------------------------------------------------
# choose year!!!!
# post curves for all lakes for a year

phis = phi
fish_names = c("crappie", 'bluegill', 'bass', 'pike', 'walleye', 'perch')
yr = 2016
df = fish_plot

# center_temp
catchability_mean = function(phis, fish_names, yr, df){
  
  tC = center_temp %>%
    filter(year(SURVEYDATE) == yr)

  lakes = df %>% 
    distinct(DOW)
  
  n_samps = dim(phis)[1]
  n_days = as.numeric(range(tC$SURVEYDATE)[2] - range(tC$SURVEYDATE)[1]+1)
  n_lakes = nrow(lakes)
  catch_post_GN = catch_post_TN = array(NA, dim = c(n_days, 6, 500, n_lakes))
  
  for(i in 1:n_lakes){
    
    lake = c(lakes[i,])$DOW
    tC_tmp = tC %>% 
      filter(DOW == lake)
    
    TN_df = tC_tmp %>% 
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
    
    GN_df = tC_tmp %>% 
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
    
    TN = as.matrix(TN_df)
    GN = as.matrix(GN_df)
    
    for(j in 501:n_samps){
      catch_post_GN[,,j-500,i] = GN %*% phis[j,,]
      catch_post_TN[,,j-500,i] = TN %*% phis[j,,]
    }
    
    print(i)
  }
  
  TN = apply(catch_post_TN, c(1,2), mean)
  GN = apply(catch_post_GN, c(1,2), mean)
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  lake = c(lakes[1,])$DOW
  tC_tmp = tC %>% 
    filter(DOW == lake)
  
  date_select = tC_tmp %>% 
    select(-c(temp_0, DOW, year))
  
  mean = tibble(as_tibble(TN), as_tibble(GN), date_select)
  
  return(mean)
}

post_curve = catchability_mean(phis, fish_names, yr, df)

catch_plot = post_curve %>%
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "value") %>% 
  mutate(Fish = rep(fnames, n()/6))


# Gear Type Scaling by Species
ggplot(catch_plot, aes(x = SURVEYDATE, y = exp(value))) +
  geom_line(aes(color = Fish), size = 1.2) +
  facet_wrap(~ Gear, scales = 'free_y', ncol = 1) +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd('2016-06-02'), ymd('2016-09-25'))) +
  scale_linetype_manual(values=c("solid", 'longdash', 'dotted')) +
  ggtitle('') +
  xlab('') +
  ylab('Catchability') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18, face = 'bold'),
        title = element_text(size = 18, face = 'bold'),
        legend.position = 'bottom',
        legend.key.size = unit(1, "cm"),
        legend.title=element_blank(),
        strip.text = element_text(size=12))

# ggsave(paste0('results/spatial_results/catchability_curves_exp_', yr, '.png'), width = 12, height = 12)
# ggsave(paste0('results/spatial_results/catchability_curves_log_', yr, '.png'), width = 12, height = 12)



# catchability curve every year -------------------------------------------
# choose lake
# post curves for a lake for all years

phis = phi
fish_names = c("crappie", 'bluegill', 'bass', 'pike', 'walleye', 'perch')
lk = '34007900'
lk = '53002800'
df = fish_plot

# lk = '53002800' # south
# lk = '16001900' # north
# lk = '18005000' # center
# 
# 
# center_temp %>% 
#   filter(DOW == lk) %>% 
#   mutate(year = as.factor(year))  %>% 
#   mutate(date = ymd(format(SURVEYDATE, paste0(2000, "-%m-%d"))))%>% 
#   ggplot(aes(x = date, y = temp_0, color = year)) +
#   geom_line()
# 
# 
# fish_plot %>%
#   group_by(DOW) %>%
#   summarize(n = n()) %>%
#   arrange(desc(n))

# center_temp
catchability_mean_lake = function(phis, fish_names, lk, df){
  
  tC = center_temp %>%
    filter(DOW == lk)
  
  # year_select = df %>% 
  #   filter(year(SURVEYDATE) == yr) %>% 
  #   select(SURVEYDATE, all_of(catch_covs), GN) %>% 
  #   mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN))
  
  lakes = df %>% 
    distinct(DOW)
  
  yr_choice = c(tC %>% distinct(year))$year
  
  n_yrs = c(tC %>% distinct(year) %>% dim())[1]
  n_samps = dim(phis)[1]
  n_days = c(tC %>% filter(year == 2000) %>% dim())[1]
  catch_post_GN = catch_post_TN = array(NA, dim = c(n_days, 6, n_yrs, n_samps))
  
  for(i in 1:n_yrs){
    
    tC_tmp = tC %>% 
      filter(year == yr_choice[i])
    
    TN_df = tC_tmp %>% 
      mutate(DOY = yday(SURVEYDATE)) %>% 
      relocate(temp_0, .after = DOY) %>% 
      mutate(DOY = yday(SURVEYDATE),
             DOY_sin_semi = sin(DOY/121 * 2*pi),
             DOY_cos_semi = cos(DOY/121 * 2*pi),
             DOY_sin_semi_temp = DOY_sin_semi * temp_0,
             DOY_cos_semi_temp = DOY_cos_semi * temp_0,
             GN = 0) %>% 
      mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
      select(-c(SURVEYDATE, DOY, DOW, year))
    
    GN_df = tC_tmp %>% 
      mutate(DOY = yday(SURVEYDATE)) %>% 
      relocate(temp_0, .after = DOY) %>% 
      mutate(DOY = yday(SURVEYDATE),
             DOY_sin_semi = sin(DOY/121 * 2*pi),
             DOY_cos_semi = cos(DOY/121 * 2*pi),
             DOY_sin_semi_temp = DOY_sin_semi * temp_0,
             DOY_cos_semi_temp = DOY_cos_semi * temp_0,
             GN = 1) %>% 
      mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
      select(-c(SURVEYDATE, DOY, DOW, year))
    
    TN = as.matrix(TN_df)
    GN = as.matrix(GN_df)
    
    for(j in 1:n_samps){
      catch_post_GN[,,i,j] = GN %*% phis[j,,]
      catch_post_TN[,,i,j] = TN %*% phis[j,,]
    }
    
  }
  
  TN = apply(catch_post_TN, c(1,2,3), mean)
  GN = apply(catch_post_GN, c(1,2,3), mean)
  
  date_select = c(tC %>% filter(year == yr_choice[1]) %>% select(-c(temp_0, DOW)))$SURVEYDATE
  TN_df = as_tibble(TN[,,1], .name_repair = ~paste0("TN_",fish_names)) %>% 
    mutate(year = yr_choice[1],
           date = date_select)
  GN_df = as_tibble(GN[,,1], .name_repair = ~paste0("GN_",fish_names)) %>% 
    mutate(year = yr_choice[1],
           date = date_select)
  
  for(i in 2:n_yrs){
    
    date_select = c(tC %>% filter(year == yr_choice[i]) %>% select(-c(temp_0, DOW)))$SURVEYDATE
    TN_tmp = as_tibble(TN[,,i], .name_repair = ~paste0("TN_",fish_names)) %>% 
      mutate(year = yr_choice[i],
             date = date_select)
    GN_tmp = as_tibble(GN[,,i], .name_repair = ~paste0("GN_",fish_names)) %>% 
      mutate(year = yr_choice[i],
             date = date_select)
    
    TN_df = rbind(TN_df, TN_tmp)
    GN_df = rbind(GN_df, GN_tmp)
    
  }
  
  return(TN_df %>% left_join(GN_df))
  
}

post_curve = catchability_mean_lake(phis, fish_names, lk, df)

catch_plot = post_curve %>%
  pivot_longer(cols = -c(year, date), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "value") %>% 
  mutate(Fish = rep(fnames, n()/6)) %>% 
  mutate(date = ymd(format(date, paste0(2000, "-%m-%d"))))



walleye_catch = catch_plot %>% 
  filter(Fish == "Walleye") %>% 
  filter(Gear == 'GN')

ggplot(walleye_catch, aes(x = date, y = value)) +
  geom_line(size = 0.9, aes(color = as.factor(year))) +
  facet_wrap(~ Gear, scales = 'free_y') +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d")

# ggsave(paste0('/Users/joshuanorth/Desktop/catchability_curves_walleye.png'), width = 12, height = 8)

ggplot(walleye_catch, aes(x = date, y = year, fill = value)) +
  geom_tile() +
  ggsci::scale_fill_material("grey", name = 'Effectiveness') +
  coord_fixed(ratio = 2) +
  xlab('') + ylab('') +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), breaks = c(2000, 2005, 2010, 2015, 2019)) +
  scale_x_date(date_breaks = '2 weeks', date_labels = "%b %d", expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.position = 'bottom',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank()) +
  guides(fill = guide_colorbar(title.position="top", title.hjust = 0.5, barwidth = 10))

# ggsave(paste0('/Users/joshuanorth/Desktop/catchability_grid_walleye.png'), width = 12, height = 8)



# catchability curve lake year ------------------------------------------------------


phis = phi
fish_names = c("crappie", 'bluegill', 'bass', 'pike', 'walleye', 'perch')
df = fish_plot

# center_temp
catchability_mean_year_lake = function(phis, fish_names, df){
  
  tC = center_temp
  
  yr_choice = c(tC %>% distinct(year))$year
  
  lakes = df %>% 
    distinct(DOW)
  
  n_yrs = c(tC %>% distinct(year) %>% dim())[1]
  n_samps = dim(phis)[1]
  n_days = c(tC %>% filter(year == 2000) %>% filter(DOW == lakes[i,]$DOW) %>% dim())[1]
  n_lakes = nrow(lakes)
  
  catch_post_GN_store = catch_post_TN_store = array(NA, dim = c(n_days, 6, n_lakes, n_yrs))
  
  for(k in 1:n_yrs){
    
    catch_post_GN = catch_post_TN = array(NA, dim = c(n_days, 6, n_lakes, n_samps))
    
    tC_yr_sub = tC %>% 
      filter(year == yr_choice[k])
    
    for(i in 1:n_lakes){
      
      catch_post_GN = catch_post_TN = array(NA, dim = c(n_days, 6, n_samps))
      
      lake = c(lakes[i,])$DOW
      tC_tmp = tC_yr_sub %>% 
        filter(DOW == lake)
  
      TN_df = tC_tmp %>% 
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
      
      GN_df = tC_tmp %>% 
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
      
      TN = as.matrix(TN_df)
      GN = as.matrix(GN_df)
      
      for(j in 1:n_samps){
        catch_post_GN[,,j] = GN %*% phis[j,,]
        catch_post_TN[,,j] = TN %*% phis[j,,]
      }
      
      # print(i)
      
      catch_post_TN_store[,,i,k] = apply(catch_post_TN, c(1,2), mean)
      catch_post_GN_store[,,i,k] = apply(catch_post_GN, c(1,2), mean)
      
    }
    
    print(k)
  }

  return(list(catch_post_TN_store = catch_post_TN_store, catch_post_GN_store = catch_post_GN_store))
  
}

post_curve = catchability_mean_year_lake(phis, fish_names, df)


ggplot() +
  geom_boxplot(aes(y = catch_post_TN_store[,1,,], x = "Black Crappie")) +
  geom_boxplot(aes(y = catch_post_TN_store[,2,,], x = "Bluegill")) +
  geom_boxplot(aes(y = catch_post_TN_store[,3,,], x = "Largemouth Bass")) +
  geom_boxplot(aes(y = catch_post_TN_store[,4,,], x = "Northern Pike")) +
  geom_boxplot(aes(y = catch_post_TN_store[,5,,], x = "Walleye")) +
  geom_boxplot(aes(y = catch_post_TN_store[,6,,], x = "Yellow Perch"))

ggplot() +
  geom_boxplot(aes(y = exp(catch_post_GN_store[,1,,]), x = "Black Crappie")) +
  geom_boxplot(aes(y = exp(catch_post_GN_store[,2,,]), x = "Bluegill")) +
  geom_boxplot(aes(y = exp(catch_post_GN_store[,3,,]), x = "Largemouth Bass")) +
  geom_boxplot(aes(y = exp(catch_post_GN_store[,4,,]), x = "Northern Pike")) +
  geom_boxplot(aes(y = exp(catch_post_GN_store[,5,,]), x = "Walleye")) +
  geom_boxplot(aes(y = exp(catch_post_GN_store[,6,,]), x = "Yellow Perch")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.position = 'bottom',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank()) +
  xlab('') + ylab('') +
  ggtitle('Effort Scaling Gill Nets') +
  theme_bw()

# write_rds(post_curve, '/Users/joshuanorth/Desktop/post_curve.rds')
# post_curve = read_rds('/Users/joshuanorth/Desktop/post_curve.rds')

catch_plot = post_curve %>%
  pivot_longer(cols = -c(year, date), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "value") %>% 
  mutate(Fish = rep(fnames, n()/6)) %>% 
  mutate(date = ymd(format(date, paste0(2000, "-%m-%d"))))


walleye_catch = catch_plot %>% 
  filter(Fish == "Walleye") %>% 
  filter(Gear == 'GN') %>% 
  group_by(year) %>% 
  distinct(date, .keep_all = T) %>%
  mutate(max = ifelse(value == max(value), 1, 0)) %>% 
  ungroup()

ggplot(walleye_catch, aes(x = date, y = value)) +
  geom_line(size = 0.9, aes(color = as.factor(year))) +
  geom_point(data = walleye_catch %>% filter(max == 1)) +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd('2000-07-15'), ymd('2000-09-01'))) +
  ylim(c(2.2, 2.5))

# ggsave(paste0('/Users/joshuanorth/Desktop/catchability_curves_walleye.png'), width = 12, height = 8)

ggplot(walleye_catch, aes(x = date, y = year, fill = value)) +
  geom_tile() +
  geom_point(data = walleye_catch %>% filter(max == 1)) +
  ggsci::scale_fill_material("grey", name = 'Effectiveness') +
  coord_fixed(ratio = 2) +
  xlab('') + ylab('') +
  theme_bw()  +
  scale_y_continuous(expand = c(0,0), breaks = c(2000, 2005, 2010, 2015, 2019)) +
  scale_x_date(date_breaks = '2 weeks', date_labels = "%b %d", expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.position = 'bottom',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank()) +
  guides(fill = guide_colorbar(title.position="top", title.hjust = 0.5, barwidth = 10))

# ggsave(paste0('/Users/joshuanorth/Desktop/catchability_grid_walleye.png'), width = 12, height = 8)


temp_cycle = center_temp %>% 
  group_by(SURVEYDATE) %>% 
  summarise(temp = mean(temp_0))

temp_cycle %>% 
  mutate(date = ymd(format(SURVEYDATE, paste0(2000, "-%m-%d"))),
         year = year(SURVEYDATE)) %>% 
  ggplot(aes(x = date, y = temp, color = as.factor(year))) +
  geom_line()


tmp = center_temp %>% 
  filter(DOW == '01000100') %>% 
  filter(year == 2000)

series = tmp$temp_0

fourier.method <- function(series, doy){
  H = length(series)
  
  #find the index associated with the annual and semi-annual cycle
  ind_an = 1 + round(H/120)
  ind_sa = 1 + round(H/(120/2))
  
  dft_tmp = fft(series)/(H/2)
  
  #amplitudes
  amps <- Mod(dft_tmp[cbind(ind_an,ind_sa)])
  
  #phases
  phases <- Arg(dft_tmp[cbind(ind_an,ind_sa)])
  
  times = 1:H
  return(tibble(temp = amps[1]*cos(2*pi*times/120 - phases[1]) + amps[2]*cos(4*pi*times/120 - phases[2]),
         date = doy))
  
}


temp_cycle = center_temp %>% 
  group_by(year, DOW) %>% 
  group_modify(~fourier.method(.x$temp_0, .x$SURVEYDATE)) %>% 
  ungroup() %>% 
  group_by(date) %>% 
  summarise(temp = mean(temp))

temp_cycle %>% 
  mutate(year = year(date),
         date = ymd(format(date, paste0(2000, "-%m-%d")))) %>% 
  ggplot(aes(x = date, y = temp, color = as.factor(year))) +
  geom_line()


# variance gamma vs random effect -----------------------------------------

n_lakes = dat$n_lakes

gdd_center = fish_dat %>% 
  select(year, DOW) %>% 
  inner_join(GDD) %>% 
  summarise(mean = mean(DD5), sd = sd(DD5))

day_choice = ymd('2016-07-15')

X_day_covs = fish_plot %>%
  select(DOW, all_of(mean_covs)) %>% 
  group_by(DOW) %>% 
  summarise_at(vars(all_of(mean_covs)), ~mean(.)) %>% 
  select(-c(DD5, secchi)) %>% 
  left_join(secchi_year %>% 
              filter(year == year(day_choice)) %>% select(-year)) %>% 
  left_join(GDD %>% 
              filter(year == year(day_choice)) %>% select(-year)) %>% 
  relocate(DD5, .after = LAKE_AREA_GIS_ACRES) %>% 
  relocate(secchi, .after = DD5) %>% 
  mutate(DD5 = (DD5 - gdd_center$mean)/gdd_center$sd)

X_day = X_day_covs %>% 
  select(-DOW) %>% 
  as.matrix()


day_post = day_post_no_omega = array(NA, dim = c(n_samps, n_lakes, K))
ratio_both = array(NA, dim = c(n_samps, n_lakes, K))

for(i in 1:n_samps){
  
  day_post[i,,] = matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + X_day %*% beta[i,,] + omega[i,,]
  day_post_no_omega[i,,] = matrix(rep(beta_0[i,], n_lakes), n_lakes, byrow = T) + X_day %*% beta[i,,]
  ratio_both[i,,] = day_post_no_omega[i,,]/day_post[i,,]
  
}

gamma_hat = apply(day_post, c(2,3), mean)
# gamma_hat = apply(day_post_no_omega, c(2,3), mean)
gamma_bar = apply(gamma_hat, 2, mean)
SSres = apply(apply(omega, c(2,3), mean)^2, 2, sum)
SStot = apply((gamma_hat - matrix(rep(gamma_bar, n_lakes), n_lakes, byrow = T))^2, 2, sum)

1-SSres/SStot

rb = apply(ratio_both, c(2,3), mean)
hist(rb, breaks = 200, xlim = c(0, 2))


cor(apply(day_post, c(2,3), mean), apply(omega, c(2,3), mean))

cor(day_post[,1,1], omega[,1,1])

corr_post = array(NA, dim = c(n_lakes, K))
for(i in 1:n_lakes){
  for(j in 1:K){
  corr_post[i,j] = cor(day_post[,i,j], omega[,i,j])
  }
}



fish_plot %>% 
  select(COMMON_NAME, gamma, omega) %>% 
  mutate(gamma = log(gamma)) %>% 
  group_by(COMMON_NAME) %>% 
  summarise(var = var(omega)/var(gamma))

apply(apply(omega, c(2,3), mean), 2, var)

apply(tau, 2, mean)^2


apply(apply(omega, c(2,3), mean), 2, var)

apply(day_post, c(2,3), mean)

apply(exp(gamma_hat), 2, var)


gamma_var_post = apply(day_post, c(2,3), var)
omega_var_post = apply(omega, c(2,3), var)
apply(ratio_both, c(2,3), mean)

round(quantile(apply(ratio_both, c(2,3), var), probs = seq(0,1,0.01)), 3)



sst = apply(day_post_no_omega, c(2,3), var)/apply(day_post, c(2,3), var)
sst = omega_var_post/gamma_var_post
sst = apply(ratio_both, c(2,3), mean)


apply(day_post, c(2,3), mean) - apply(day_post_no_omega, c(2,3), mean)

apply(omega, c(2,3), mean)/(mean(apply(day_post, c(2,3), mean)) - apply(day_post, c(2,3), mean))




fish_names = c("crappie", 'bluegill', 'bass', 'pike', 'walleye', 'perch')

colnames(sst) = paste0(fish_names)

lake_loc = fish_plot %>% 
  distinct(DOW, .keep_all = T) %>% 
  select(DOW, lat, lon)

sse_comp = X_day_covs %>% 
  select(DOW) %>% 
  cbind(sst) %>% 
  tibble() %>% 
  left_join(lake_loc) %>% 
  pivot_longer(-c(DOW, lat, lon), names_to = 'fish', values_to = 'value')



ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = sse_comp,
             aes(x = lon, y = lat, color = value), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ fish,
             labeller = labeller(fish = c('crappie' = 'Black Crappie',
                                          'bluegill' = 'Bluegill',
                                          'pike' = 'Northern Pike',
                                          'perch' = 'Yellow Perch',
                                          'bass' = 'Largemouth Bass',
                                          'walleye' = 'Walleye'))) +
  ggtitle(paste0('Predicted Catch for 1 Unit Effort on ', day_choice, ' - TN')) +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))
# summary of data ---------------------------------------------------------

fish_dat %>% 
  select(COMMON_NAME, TOTAL_CATCH, CPUE) %>% 
  rename('Total' = 'TOTAL_CATCH') %>% 
  group_by(COMMON_NAME) %>% 
  summarise_at(vars(Total, CPUE), list(min = ~min(., na.rm = T),
                                             max = ~max(., na.rm = T),
                                             median = ~median(., na.rm = T),
                                             mean = ~mean(., na.rm = T),
                                             sd = ~sd(., na.rm = T))) %>% 
  pivot_longer(-COMMON_NAME, names_to = c('Variable', 'Statistic'), names_sep = '_', values_to = 'Test') %>% 
  pivot_wider(names_from = 'Statistic', values_from = 'Test') %>% 
  xtable()


fish_dat %>% 
  summarise_at(vars(SURVEYDATE), list(min = ~min(., na.rm = T),
                                      max = ~max(., na.rm = T),
                                      median = ~median(., na.rm = T),
                                      mean = ~mean(., na.rm = T),
                                      sd = ~sd(., na.rm = T)))


fish_dat %>% 
  select(DOW, all_of(mean_covs)) %>% 
  rename('MaxDepth' = 'MAX_DEPTH_FEET',
         'LakeArea' = 'LAKE_AREA_GIS_ACRES',
         'GDD' = 'DD5') %>% 
  distinct(DOW, .keep_all = T) %>% 
  summarise_at(vars(MaxDepth:wetlands), list(min = ~min(., na.rm = T),
                                      max = ~max(., na.rm = T),
                                      median = ~median(., na.rm = T),
                                      mean = ~mean(., na.rm = T),
                                      sd = ~sd(., na.rm = T))) %>% 
  pivot_longer(MaxDepth_min:wetlands_sd, names_to = c('Variable', 'Statistic'), names_sep = '_', values_to = 'Test') %>% 
  pivot_wider(names_from = 'Statistic', values_from = 'Test') %>% 
  xtable()
  
fish_dat %>% 
  select(DOW, year) %>% 
  left_join(GDD) %>% 
  summarise_at(vars(DD5), list(min = ~min(., na.rm = T),
                               max = ~max(., na.rm = T),
                               median = ~median(., na.rm = T),
                               mean = ~mean(., na.rm = T),
                               sd = ~sd(., na.rm = T))) %>% 
  xtable()


fish_dat %>% 
  select(DOW, all_of(mean_covs)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  rename('MaxDepth' = 'MAX_DEPTH_FEET',
         'LakeArea' = 'LAKE_AREA_GIS_ACRES',
         'GDD' = 'DD5') %>% 
  distinct(DOW, .keep_all = T) %>% 
  summarise_at(vars(MaxDepth:wetlands), list(min = ~min(., na.rm = T),
                                             max = ~max(., na.rm = T),
                                             median = ~median(., na.rm = T),
                                             mean = ~mean(., na.rm = T),
                                             sd = ~sd(., na.rm = T))) %>% 
  pivot_longer(MaxDepth_min:wetlands_sd, names_to = c('Variable', 'Statistic'), names_sep = '_', values_to = 'Test') %>% 
  pivot_wider(names_from = 'Statistic', values_from = 'Test') %>% 
  xtable()


fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  summarize(number = n()) %>% 
  summarize(min = min(number), max = max(number), mean = mean(number), median = median(number), sd = sd(number)) %>% 
  xtable(caption = 'Number of repeat sampler per lake.')

# average survey date day-month
fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  mutate(SURVEYDATE=ymd(format(SURVEYDATE, "2016-%m-%d"))) %>% 
  ungroup() %>% 
  summarize(min = min(SURVEYDATE), max = max(SURVEYDATE), mean = mean(SURVEYDATE), median = median(SURVEYDATE), sd = sd(SURVEYDATE)) %>% 
  mutate_all(~format(., format="%m-%d")) %>% 
  xtable(caption = 'Summary of sampling calendar day when surveys were conducted.')

# average survey date doth  
fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  mutate(day=yday(SURVEYDATE)) %>% 
  ungroup() %>% 
  summarize(min = min(day), max = max(day), mean = mean(day), median = median(day), sd = sd(day)) %>% 
  xtable(caption = 'Summary of sampling day of the year when surveys were conducted.')
