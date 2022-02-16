
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


# load data ---------------------------------------------------------------

# fish_dat = read_csv('data/fish_dat.csv')
# 
# fish_dat = fish_dat %>% 
#   mutate(DOW = as.factor(DOW),
#          COMMON_NAME = as.factor(COMMON_NAME)) %>% 
#   arrange(DOW, year, COMMON_NAME)

# 308 lakes
# only lakes observed 4 or more times
# fish_dat = fish_dat %>%
#   group_by(DOW) %>%
#   filter(n() >= 6*12) %>%
#   ungroup() %>%
#   mutate_at(vars(DOW), ~droplevels(.))
# 
# length(table(fish_dat$DOW))


# mean_covs = colnames(fish_dat)[c(7, 9, 23, 24, 13:15,17)]
# temporal_covs = colnames(fish_dat)[c(23, 24)]
# mean_covs_log = colnames(fish_dat)[c(7, 9)]
# mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
# catch_covs = colnames(fish_dat)[c(25, 27:30)] # temp doy interaction
# gear_types = colnames(fish_dat)[c(21, 22)]


# with AG
# mean_covs = colnames(fish_dat)[c(7, 9, 23, 24, 13:15)]
# temporal_covs = colnames(fish_dat)[c(23, 24)]
# mean_covs_log = colnames(fish_dat)[c(7, 9)]
# mean_covs_logit = colnames(fish_dat)[c(13:15)]
# catch_covs = colnames(fish_dat)[c(25, 27:30)] # temp doy interaction
# gear_types = colnames(fish_dat)[c(21, 22)]


# without AG
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

out = stan(file = 'stan/spatial_jsdm_effort_scaling.stan', data = dat, iter = 2000, warmup = 1000, chains = 1, cores = 1, refresh = 10) # our model

# saveRDS(out, "C:/Users/jsnow/Desktop/stan_lake.rds")


out = read_rds("C:/Users/jsnow/Desktop/stan_lake.rds")

chains = extract(out, permuted = T)
names(chains)
lapply(chains, dim)


phis = chains$phi[100,,]
zm = dat$Zstar %*% chains$phi[sample(1:999, size = 1),,]
mu = mean(zm)
sig = var(c(zm))/2
mean(exp(zm - mu - sig))

apply(chains$omega, 3, mean)


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

par(mfrow = c(3,4))
for(i in 1:12){
  plot(chains$phi[,i,5], type = 'l', main = i)
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










# relative abundance ------------------------------------------------------

beta_0 = chains$beta_0
beta = chains$beta
phi = chains$phi
omega = chains$omega
sigma_species = chains$Sigma_species
tau = chains$tau

b_names = colnames(dat$X)
phi_names = colnames(dat$Z)

static = read_csv('data/Static_MN_Data_JE.csv') %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
effort = read_csv('data/MN_catch_survey.csv') %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))

fish_dat = fish_dat %>% 
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN)

b_hat = t(apply(beta, c(2,3), mean))
colnames(b_hat) = b_names

phi_hat = t(apply(phi, c(2,3), mean))
colnames(phi_hat) = phi_names

# construct omega
K = length(dat$fish_names)
omega_hat = apply(omega, c(2,3), mean)
ind_array = tibble(id = dat$lake_id, as_tibble(omega_hat, .name_repair = ~LETTERS[1:K]))
lake_array = tibble(id = dat$lake_index)
omega_hat = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))


# lam_hat = list()
# for(k in 1:K){
#   lam_hat[[k]] = exp(beta_0[k] + dat$X %*% b_hat[k,] + omega_hat[,k] + dat$Z %*% phi_hat[k,])
# }


b0 = apply(beta_0, 2, mean)
lam_hat = list()
for(k in 1:K){
  lam_hat[[k]] = exp(b0[k] + dat$X %*% b_hat[k,] + omega_hat[,k] + dat$Z %*% phi_hat[k,])
}


rel_abun = rbind(fish_dat %>% filter(COMMON_NAME == 'black crappie') %>% mutate(Abundance = lam_hat[[1]]),
                 fish_dat %>% filter(COMMON_NAME == 'bluegill') %>% mutate(Abundance = lam_hat[[2]]),
                 fish_dat %>% filter(COMMON_NAME == 'largemouth bass') %>% mutate(Abundance = lam_hat[[3]]),
                 fish_dat %>% filter(COMMON_NAME == 'northern pike') %>% mutate(Abundance = lam_hat[[4]]),
                 fish_dat %>% filter(COMMON_NAME == 'walleye') %>% mutate(Abundance = lam_hat[[5]]),
                 fish_dat %>% filter(COMMON_NAME == 'yellow perch') %>% mutate(Abundance = lam_hat[[6]]))


rel_abun = rel_abun %>% 
  mutate(Fish = str_replace_all(COMMON_NAME, " ", "_")) %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW')

lats = range(rel_abun$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(rel_abun$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = rel_abun, 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~Fish, 
             labeller = labeller(Fish = c('black_crappie' = 'Black Crappie',
                                          'bluegill' = 'Bluegill',
                                          'northern_pike' = 'Northern Pike',
                                          'yellow_perch' = 'Yellow Perch',
                                          'largemouth_bass' = 'Largemouth Bass',
                                          'walleye' = 'Walleye'))) +
  ggtitle("Relative Abundance") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))


rel_abun %>% 
  dplyr::select(-c(MAX_DEPTH_FEET:DOY_cos_semi_temp)) %>% 
  group_by(COMMON_NAME) %>% 
  filter(CPUE == max(CPUE))

rel_abun %>% 
  dplyr::select(-c(MAX_DEPTH_FEET:DOY_cos_semi_temp)) %>% 
  group_by(COMMON_NAME) %>% 
  filter(Abundance == max(Abundance))



rel_abun %>% 
  ggplot(., aes(x = CPUE, y = Abundance)) +
  geom_point() +
  facet_wrap(~COMMON_NAME, scales = 'free')


# spatial random effect ---------------------------------------------------



# construct omega
K = length(dat$fish_names)
omega_hat = apply(omega, c(2,3), mean)
ind_array = tibble(id = dat$lake_id, as_tibble(omega_hat, .name_repair = ~LETTERS[1:K]))
lake_array = tibble(id = dat$lake_index)
omega_hat = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))

spat_pat = rbind(fish_dat %>% filter(COMMON_NAME == 'black crappie') %>% mutate(Abundance = omega_hat[,1]),
                 fish_dat %>% filter(COMMON_NAME == 'bluegill') %>% mutate(Abundance = omega_hat[,2]),
                 fish_dat %>% filter(COMMON_NAME == 'largemouth bass') %>% mutate(Abundance = omega_hat[,3]),
                 fish_dat %>% filter(COMMON_NAME == 'northern pike') %>% mutate(Abundance = omega_hat[,4]),
                 fish_dat %>% filter(COMMON_NAME == 'walleye') %>% mutate(Abundance = omega_hat[,5]),
                 fish_dat %>% filter(COMMON_NAME == 'yellow perch') %>% mutate(Abundance = omega_hat[,6]))

spat_pat = spat_pat %>% 
  mutate(Fish = str_replace_all(COMMON_NAME, " ", "_")) %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW') %>% 
  group_by(COMMON_NAME) %>% 
  distinct(DOW, .keep_all = T) %>% 
  ungroup()

lats = range(spat_pat$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(spat_pat$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = spat_pat, 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance), size = 3) +
  scale_color_gradientn(colors = fields::two.colors(start = 'blue', end = 'red', middle = 'white')) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~Fish, 
             labeller = labeller(Fish = c('black_crappie' = 'Black Crappie',
                                          'bluegill' = 'Bluegill',
                                          'northern_pike' = 'Northern Pike',
                                          'yellow_perch' = 'Yellow Perch',
                                          'largemouth_bass' = 'Largemouth Bass',
                                          'walleye' = 'Walleye'))) +
  ggtitle("Spatial Random Effect") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))





# look at data summaries --------------------------------------------------


select(all_of(mean_covs)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.))

fish_dat %>% 
  select(all_of(mean_covs), COMMON_NAME, TOTAL_CATCH) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  filter(TOTAL_CATCH != 0) %>% 
  group_by(COMMON_NAME) %>% 
  summarise_at(vars(all_of(mean_covs)), mean)


fish_dat %>% 
  select(all_of(mean_covs), COMMON_NAME, TOTAL_CATCH) %>% 
  filter(TOTAL_CATCH != 0) %>% 
  group_by(COMMON_NAME) %>% 
  summarise_at(vars(all_of(mean_covs)), mean)



fish_dat %>% 
  select(DD5, COMMON_NAME, TOTAL_CATCH) %>% 
  filter(TOTAL_CATCH != 0) %>% 
  group_by(COMMON_NAME) %>% 
  summarise(mean(DD5))

fish_dat %>% 
  select(DD5, COMMON_NAME, TOTAL_CATCH) %>% 
  filter(TOTAL_CATCH != 0) %>% 
  group_by(COMMON_NAME) %>% 
  summarise(quantile(DD5, probs = c(0.025)))

fish_dat %>% 
  select(DD5, COMMON_NAME, TOTAL_CATCH) %>% 
  filter(TOTAL_CATCH != 0) %>% 
  group_by(COMMON_NAME) %>% 
  summarise(quantile(DD5, probs = c(0.975)))



fish_dat %>% 
  select(DD5, COMMON_NAME, TOTAL_CATCH) %>% 
  filter(TOTAL_CATCH != 0) %>% 
  group_by(COMMON_NAME) %>% 
  summarise(sd(DD5))






fish_dat %>% 
  select(all_of(mean_covs), COMMON_NAME, TOTAL_CATCH) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  filter(TOTAL_CATCH != 0) %>% 
  group_by(COMMON_NAME) %>% 
  summarise(mean = mean(MAX_DEPTH_FEET),
            sd = sd(MAX_DEPTH_FEET),
            lower = quantile(MAX_DEPTH_FEET, probs = 0.2), 
            upper = quantile(MAX_DEPTH_FEET, probs = 0.8),
            min = min(MAX_DEPTH_FEET),
            max = max(MAX_DEPTH_FEET))

fish_dat %>% 
  select(all_of(mean_covs), COMMON_NAME, TOTAL_CATCH) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  filter(TOTAL_CATCH != 0) %>% 
  group_by(COMMON_NAME) %>% 
  summarise(mean = mean(DD5),
            sd = sd(DD5),
            lower = quantile(DD5, probs = 0.2), 
            upper = quantile(DD5, probs = 0.8),
            min = min(DD5),
            max = max(DD5))

fish_dat %>% 
  select(all_of(mean_covs), COMMON_NAME, TOTAL_CATCH) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  summarise(mean = mean(DD5),
            sd = sd(DD5),
            lower = quantile(DD5, probs = 0.2), 
            upper = quantile(DD5, probs = 0.8),
            min = min(DD5),
            max = max(DD5))



out = read_rds('C:/Users/jsnow/Desktop/stan_lake_random_effect.rds')
chains = extract(out)
chains_no_catch = extract(out_no_catch)
names(chains)
lapply(chains, dim)

beta_0 = apply(chains$beta_0, 2, mean)
beta = apply(chains$beta, c(2,3), mean)

blue_mean = c(beta_0[2], beta[,2])


X = fish_dat %>% 
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  filter(COMMON_NAME == 'bluegill') %>% 
  select(DOW, year, all_of(mean_covs)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
  filter(DOW == '34007900') %>% 
  filter(year == 2016) %>% 
  select(-c(DOW, year)) %>% 
  as.matrix()

X = c(1, X[1,])

x1 = X # depth 0.656
x2 = X # depth -0.522
x3 = X # gdd 1
x4 = X # gdd -1
x1[2] = 0.656
x2[2] = -0.522
x3[4] = 1
x4[4] = -1


exp(sum(x1 * blue_mean)) # depth 0.656
exp(sum(x2 * blue_mean)) # depth -0.522
exp(sum(x3 * blue_mean)) # gdd 1
exp(sum(x4 * blue_mean)) # gdd -1

exp(sum(x1 * blue_mean))/exp(sum(x2 * blue_mean))
exp(sum(x3 * blue_mean))/exp(sum(x4 * blue_mean))

exp(sum(x4 * blue_mean))/exp(sum(x3 * blue_mean))






sort(table(fish_dat$DOW), decreasing = T)

fish_plot = fish_dat %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW') %>% 
  rename(lat = LAKE_CENTER_LAT_DD5,
         lon = LAKE_CENTER_LONG_DD5) 

fish_plot %>% 
  filter(DOW == '34007900')



lats = range(fish_plot$lat, na.rm = T)
lons = range(fish_plot$lon, na.rm = T)

usa = st_as_sf(maps::map("state", fill= TRUE, plot = FALSE))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_plot %>% 
               filter(DOW == '34007900'), 
             aes(x = lon, y = lat, color = TOTAL_CATCH), size = 1.5)



dat$fish_names
head(dat$X)














