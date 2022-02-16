
library(tidyverse)
library(rstan)
library(fields)
library(Matrix)
library(lubridate)
library(readxl)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(mvtnorm)
library(scales)
library(progress)
library(sparklyr)
library(stringr)
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
library(ggsci)

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
  left_join(static, by = 'DOW') %>% # join effort and static
  select(DOW, COMMON_NAME, TOTAL_CATCH, EFFORT, GEAR, CPUE, SURVEYDATE, MAX_DEPTH_FEET, mean.gdd, LAKE_AREA_GIS_ACRES, 
         LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING) %>% # select wanted columns
  mutate(SURVEYDATE = str_split(SURVEYDATE, ' ', simplify = T)[,1]) %>% 
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME),
         GEAR = as.factor(GEAR),
         SURVEYDATE = mdy(SURVEYDATE)) %>% # format variables
  filter(complete.cases(.)) %>% # only complete observations
  filter(COMMON_NAME != 'white sucker',
         COMMON_NAME != 'smallmouth bass',
         COMMON_NAME != 'rock bass') %>% # remove unused fish
  filter(SURVEYDATE >= '1993-01-01') %>% # subset to after 1993
  mutate(GEAR = as.character(GEAR),
         GEAR = ifelse(GEAR == 'GDE' | GEAR == 'GSH', 'GN', GEAR),
         GEAR = as.factor(GEAR),
         ind = 1) %>% # fix GN gears
  group_by(DOW, COMMON_NAME, SURVEYDATE, GEAR) %>% 
  summarise_all(~sum(.)) %>% # handle fringe case with multiple obs in year lake (~48 instances)
  ungroup() %>% 
  mutate_at(vars(MAX_DEPTH_FEET:LAKE_CENTER_UTM_NORTHING), ~ ./ind) %>% # correct fringe summaries
  select(-ind) %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW),
         GEAR = droplevels(GEAR)) %>% # drop unused factors
  arrange(DOW) %>% 
  mutate(year = year(SURVEYDATE)) # presentation

# add in land use
fish_dat = fish_dat %>% 
  left_join(land %>% select(-year), by = c('DOW')) %>% 
  filter(complete.cases(.))

# all combos of fish, survey date, and gear (when fish wasn't caught should be 0)
all <- fish_dat %>% 
  group_by(DOW) %>% 
  tidyr::expand(COMMON_NAME, SURVEYDATE, GEAR) %>% 
  ungroup() %>% 
  arrange(DOW)


fish_dat <- fish_dat %>% 
  right_join(all) %>% # expand data to include 0 caught cases
  mutate(EFFORT = coalesce(EFFORT, 0L),
         TOTAL_CATCH = coalesce(TOTAL_CATCH, 0L),
         CPUE = coalesce(CPUE, 0L)) %>% # change NA's from join to 0
  mutate_at(vars(MAX_DEPTH_FEET:lake_elevation_m), list(~coalesce(., 0L))) %>%  # change NA's from join to 0
  group_by(GEAR, SURVEYDATE, DOW) %>% 
  mutate(EFFORT = max(EFFORT)) %>% # make sure the effort is the same across sample numbers
  ungroup() %>% 
  group_by(DOW, SURVEYDATE) %>% 
  mutate_at(vars(MAX_DEPTH_FEET:lake_elevation_m), ~ mean(.[!is.na(.) & . != 0])) %>% # make sure variables are not zero
  mutate_at(vars(MAX_DEPTH_FEET:lake_elevation_m), ~ ifelse(is.na(.), 0, .)) %>%  # make sure variables are not zero
  ungroup() %>% 
  mutate(TN = ifelse(GEAR == 'TN', 1, 0),
         GN = ifelse(GEAR == 'GN', 1, 0)) %>% # number gear
  select(-GEAR) %>% 
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW)) %>%
  mutate(CPUE = TOTAL_CATCH/EFFORT,
         CPUE = ifelse(is.na(CPUE), 0, CPUE)) %>% 
  arrange(DOW)


# center temp by subset of lakes - used in covariate construction
center_temp = temp %>% 
  select(SURVEYDATE, temp_0, DOW) %>% 
  mutate(year = year(SURVEYDATE)) %>% 
  filter(year >= 2000) %>% 
  mutate(filter_date = ymd(format(SURVEYDATE, "2016-%m-%d"))) %>% 
  filter(filter_date > ymd('2016-06-01'),
         filter_date < ymd('2016-09-30')) %>% # center based on "fishing season"
  select(-filter_date) %>% 
  group_by(DOW, year) %>% 
  mutate(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0)) %>% # standardized
  ungroup()

# add in GDD, secchi, and centered temp
fish_dat = fish_dat %>% 
  inner_join(GDD) %>% # add GDD
  inner_join(secchi_year, by = c('DOW', 'year')) %>% # add secchi
  filter(year >= 2000) %>% # filter year
  mutate(filter_date = ymd(format(SURVEYDATE, "2016-%m-%d"))) %>% 
  filter(filter_date > ymd('2016-06-01'),
         filter_date < ymd('2016-09-30')) %>% # filter fishing season
  select(-filter_date) %>% 
  inner_join(center_temp %>% 
               select(SURVEYDATE, temp_0, DOW)) %>% # add temp
  mutate(DD5 = (DD5 - mean(DD5))/sd(DD5)) %>% # standardized GDD
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>% # create semi-annual fishing season
  mutate(DOW = as.factor(DOW)) %>% 
  arrange(DOW, year, COMMON_NAME)


# load data ---------------------------------------------------------------

# fish_dat = read_csv('data/fish_dat.csv')
# 
# fish_dat = fish_dat %>% 
#   mutate(DOW = as.factor(DOW),
#          COMMON_NAME = as.factor(COMMON_NAME)) %>% 
#   arrange(DOW, year, COMMON_NAME)


# without AG - select values to pass into pars
mean_covs = colnames(fish_dat)[c(7, 9, 23, 24, 14:15)]
temporal_covs = colnames(fish_dat)[c(23, 24)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(14:15)]
catch_covs = colnames(fish_dat)[c(25, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]


# set up stan parameters --------------------------------------------------


create_pars <- function(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, center_temp){
  
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

rm(all, effort, GDD, land, secchi_year, static, temp, create_pars)

out = stan(file = 'stan/spatial_jdsm_without_effort_scaling.stan', data = dat, iter = 2000, warmup = 1000, chains = 1, cores = 1)

# saveRDS(out, "C:/Users/jsnow/Desktop/stan_jdsm_without_effort_scaling.rds")

# out = read_rds('/Users/joshuanorth/Desktop/stan_jdsm_without_effort_scaling.rds')


chains = extract(out, permuted = T)
names(chains)
lapply(chains, dim)



b_names = colnames(dat$X)
b_hat = t(apply(chains$beta, c(2,3), mean))
fnames = dat$fish_names %>% str_replace(., ' ', '_')

round(b_hat, 3) %>%
  as_tibble(.name_repair = ~b_names) %>%
  mutate(Species = dat$fish_names) %>%
  relocate(Species)


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

par(mfrow = c(2,3))
for(i in 1:6){
  plot(chains$tau[,i], type = 'l', main = i)
}


plot(chains$lp__, type = 'l')


