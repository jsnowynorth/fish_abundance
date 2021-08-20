##----------------------------------------------------
## Name: Joshua North
##
## Date: 21/07/2021
##
## Project: Fish Abundance
##
## Objective: Multivariate fish abundance model. New survey information
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
# library(MCMCpack)
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
fish_dat = fish_dat %>% 
  inner_join(GDD) %>% 
  mutate(DD5 = (DD5 - mean(DD5))/sd(DD5)) %>% 
  inner_join(temp %>% 
               select(SURVEYDATE, temp_0, DOW) %>% 
               mutate(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0))) %>% 
  inner_join(secchi_year, by = c('DOW', 'year'))

fish_dat = fish_dat %>% 
  inner_join(GDD) %>% 
  mutate(DD5 = (DD5 - mean(DD5))/sd(DD5)) %>% 
  inner_join(temp %>% 
               select(SURVEYDATE, temp_0, DOW)) %>% 
  inner_join(secchi_year, by = c('DOW', 'year')) %>% 
  filter(year >= 2000) %>% 
  mutate(filter_date = ymd(format(SURVEYDATE, "2016-%m-%d"))) %>% 
  filter(filter_date > ymd('2016-06-01'),
         filter_date < ymd('2016-09-30')) %>% 
  select(-filter_date) %>% 
  mutate(temp_0 = (temp_0 - mean(temp_0))/sd(temp_0)) %>% 
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>% 
  mutate(DOW = as.factor(DOW)) %>% 
  arrange(DOW, year, COMMON_NAME)


# write_csv(fish_dat, 'data/fish_dat.csv')

# model -------------------------------------------------------------------

source('R/lewis_code/lewis_model.R')
# source('R/lewis_code/lewis_model_no_spatial.R')

# spatial covs=
mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:15,17, 25)]
temporal_covs = colnames(fish_dat)[c(23, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]


fish_dat %>% 
  select(all_of(mean_covs)) %>% 
  rename('depth' = 'MAX_DEPTH_FEET',
         'area' = 'LAKE_AREA_GIS_ACRES') %>% 
  summarise_all(list(min = ~min(.),
                     max = ~max(.), 
                     mean = ~mean(.),
                     median = ~median(.),
                     sd = ~sd(.))) %>% 
  pivot_longer(depth_min:secchi_sd, names_to = c('Variable', 'Summary'), names_sep = '_', values_to = 'Value') %>% 
  pivot_wider(names_from = 'Summary', values_from = 'Value')


fish_dat %>% 
  select(DOW:CPUE) %>% 
  rename('fish' = 'COMMON_NAME',
         'total' = 'TOTAL_CATCH') %>% 
  rename_all(~str_to_lower(.)) %>% 
  group_by(fish) %>% 
  summarise_at(vars(total, cpue), list(min = ~min(.),
                                       max = ~max(.), 
                                       mean = ~mean(.),
                                       median = ~median(.),
                                       sd = ~sd(.))) %>% 
  pivot_longer(total_min:cpue_sd, names_to = c('Variable', 'Summary'), names_sep = '_', values_to = 'Value') %>% 
  pivot_wider(names_from = 'Summary', values_from = 'Value')


fish_dat %>% 
  select(DOW:CPUE) %>% 
  rename('fish' = 'COMMON_NAME',
         'total' = 'TOTAL_CATCH') %>% 
  rename_all(~str_to_lower(.)) %>% 
  group_by(fish) %>% 
  summarise_at(vars(total, cpue), list(min = ~quantile(., probs = 0.02),
                                       max = ~quantile(., probs = 0.98), 
                                       mean = ~mean(.),
                                       median = ~median(.),
                                       sd = ~sd(.))) %>% 
  pivot_longer(total_min:cpue_sd, names_to = c('Variable', 'Summary'), names_sep = '_', values_to = 'Value') %>% 
  pivot_wider(names_from = 'Summary', values_from = 'Value')


fish_dat %>% 
  select(DOW:CPUE) %>% 
  rename('fish' = 'COMMON_NAME',
         'total' = 'TOTAL_CATCH') %>% 
  rename_all(~str_to_lower(.)) %>% 
  ggplot(., aes(x = total, y = effort)) +
  geom_point() +
  facet_wrap(~fish, scales = 'free')

fish_dat %>% 
  select(DOW:CPUE) %>% 
  rename('fish' = 'COMMON_NAME',
         'total' = 'TOTAL_CATCH') %>% 
  rename_all(~str_to_lower(.)) %>% 
  ggplot(., aes(x = total, y = effort)) +
  geom_point() +
  facet_wrap(~fish)


fish_dat %>% 
  select(DOW:CPUE) %>% 
  rename('fish' = 'COMMON_NAME',
         'total' = 'TOTAL_CATCH') %>% 
  group_by(fish) %>% 
  filter(CPUE < quantile(CPUE, probs = 0.99)) %>% 
  ungroup() %>% 
  rename_all(~str_to_lower(.)) %>% 
  ggplot(., aes(x = total, y = effort)) +
  geom_point() +
  facet_wrap(~fish, scales = 'free')


fish_dat %>% 
  select(DOW:CPUE) %>% 
  rename('fish' = 'COMMON_NAME',
         'total' = 'TOTAL_CATCH') %>% 
  group_by(fish) %>% 
  filter(total < quantile(total, probs = 0.97)) %>% 
  ungroup() %>% 
  rename_all(~str_to_lower(.)) %>% 
  ggplot(., aes(x = total, y = effort)) +
  geom_point() +
  facet_wrap(~fish, scales = 'free')