
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
# library(Vizumap)
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
  mutate(DD5 = rollmean(DD5, k = 5, align = 'right', fill = NA)) %>% 
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
  select(-filter_date)

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
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/121 * 2*pi),
         DOY_cos_semi = cos(DOY/121 * 2*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>% # create semi-annual fishing season
  mutate(DOW = as.factor(DOW)) %>% 
  arrange(DOW, year, COMMON_NAME)


fish_dat = fish_dat %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW') %>% 
  rename(lat = LAKE_CENTER_LAT_DD5,
         lon = LAKE_CENTER_LONG_DD5)


# save data ---------------------------------------------------------------

# write_csv(center_temp, 'data/raw_temp.csv')
# write_csv(fish_dat, 'data/fish_dat_raw.csv')

