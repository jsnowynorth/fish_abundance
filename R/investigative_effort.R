##----------------------------------------------------
## Name: Joshua North
##
## Date: 10/19/2020
##
## Project: Fish Abundance
##
## Objective: Investigative plots
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



# load in data ------------------------------------------------------------


static = read_csv('data/Static_MN_Data_JE.csv')
effort = read_csv('data/MN_fish_catch_effort_all_lakes_years_surveys_JE.csv')



# join data ---------------------------------------------------------------


static = static %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
effort = effort %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))

dat = effort %>% 
  left_join(static, by = 'DOW') %>% 
  mutate(SURVEYDATE = mdy(SURVEYDATE)) %>% 
  select(-c(SURVEY_ID, LKNAME, CTY:SURVEYTYPE, AVG_MEASURED_WT_LBS:ALT_LAKE_NAME)) %>% 
  group_by(DOW, COMMON_NAME, GEAR) %>% 
  summarise_at(vars(CPUE:LAKE_CENTER_UTM_NORTHING), list(mean)) %>% 
  ungroup()

# dat = effort %>% 
#   left_join(static, by = 'DOW') %>% 
#   mutate(SURVEYDATE = mdy(SURVEYDATE)) %>% 
#   select(-c(SURVEY_ID, LKNAME, CTY:SURVEYTYPE, AVG_MEASURED_WT_LBS:ALT_LAKE_NAME)) %>% 
#   group_by(DOW, COMMON_NAME, GEAR) %>% 
#   summarise_all(~mean(., na.rm = T))



# looking at area to littoral surface -------------------------------------

ggplot(data = dat, aes(x = log(LAKE_AREA_GIS_ACRES), fill = GEAR)) +
  geom_density(alpha = 0.5) +
  ggtitle('Log Lake Area GIS acres by Gear type') +
  xlab('') +
  facet_wrap(~GEAR)

ggplot(data = dat, aes(x = log(LITTORAL_AREA_ACRES), fill = GEAR)) +
  geom_density(alpha = 0.5) +
  ggtitle('Log Littoral area acres by Gear type') +
  xlab('') +
  facet_wrap(~GEAR)


sub_dat <- dat %>% 
  select(LAKE_AREA_GIS_ACRES, LITTORAL_AREA_ACRES, GEAR) %>% 
  rename(Lake = LAKE_AREA_GIS_ACRES,
         Littoral = LITTORAL_AREA_ACRES) %>% 
  pivot_longer(-GEAR, names_to = 'Measure', values_to = 'Area')

ggplot(data = sub_dat, aes(x = log(Area), fill = Measure)) +
  geom_density(alpha = 0.5) +
  ggtitle('Log Lake and Log Littoral Area by Gear Type') +
  xlab('') +
  ylab('Density') +
  facet_wrap(~GEAR)

ggplot(data = dat, aes(x = log(LAKE_AREA_GIS_ACRES), y = log(EFFORT), color = GEAR)) +
  geom_point(size = 1) +
  facet_wrap(~GEAR)

  


# density plots -----------------------------------------------------------

ggplot(data = dat, aes(x = log(CPUE), fill = COMMON_NAME)) +
  geom_density(alpha = 0.5) +
  ggtitle('Log CPUE by Species') +
  xlab('')

ggplot(data = dat, aes(x = log(CPUE), fill = GEAR)) +
  geom_density(alpha = 0.5) +
  ggtitle('Log CPUE by gear') +
  xlab('')

ggplot(data = dat, aes(x = log(TOTAL_CATCH), fill = COMMON_NAME)) +
  geom_density(alpha = 0.5) +
  ggtitle('Log total catch by Species') +
  xlab('')

ggplot(data = dat, aes(x = PDIST_cum, fill = COMMON_NAME)) +
  geom_density(alpha = 0.5) +
  ggtitle('PDIST by Species') +
  xlab('')

ggplot(data = dat, aes(x = log(EFFORT), fill = COMMON_NAME)) +
  geom_density(alpha = 0.5) +
  ggtitle('Log effort by Species') +
  xlab('')

ggplot(data = dat, aes(x = log(MAX_DEPTH_FEET), fill = COMMON_NAME)) +
  geom_density(alpha = 0.5) +
  ggtitle('Log max depth by Species') +
  xlab('')

ggplot(data = dat, aes(x = log(SHORE_LENGTH_MILES), fill = COMMON_NAME)) +
  geom_density(alpha = 0.5) +
  ggtitle('Log shore length by Species') +
  xlab('')

ggplot(data = dat, aes(x = secchi.m, fill = COMMON_NAME)) +
  geom_density(alpha = 0.5) +
  ggtitle('Secchi depth by Species') +
  xlab('')


# 2d density plots --------------------------------------------------------

ggplot(data = dat, aes(x = log(MEAN_DEPTH_FEET), y = log(SHORE_LENGTH_MILES), z = log(CPUE), color = COMMON_NAME)) +
  stat_density2d(size = 1)

ggplot(data = dat, aes(x = log(MEAN_DEPTH_FEET), y = log(SHORE_LENGTH_MILES), z = log(CPUE), color = GEAR)) +
  stat_density2d(size = 1)

ggplot(data = dat, aes(x = PDIST_cum, y = log(LITTORAL_AREA_ACRES), z = log(CPUE), color = COMMON_NAME)) +
  stat_density2d(size = 1)


# plots by longitude and latitude -----------------------------------------

lats = range(dat$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(dat$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = log(CPUE))) +
  scale_color_viridis(name = "Log CPUE", option = 'inferno', direction = 1, na.value = NA) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ GEAR)


ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = dat %>% filter(COMMON_NAME == 'northern pike'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = log(TOTAL_CATCH))) +
  scale_color_viridis(name = "Log Total Catch", option = 'inferno', direction = 1, na.value = NA) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ GEAR)


ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = log(CPUE))) +
  scale_color_viridis(name = "Log CPUE", option = 'inferno', direction = 1, na.value = NA) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ COMMON_NAME)

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = secchi.m)) +
  scale_color_viridis(name = "Secchi", option = 'inferno', direction = 1, na.value = NA) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ COMMON_NAME)

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = log(MAX_DEPTH_FEET))) +
  scale_color_viridis(name = "log max depth", option = 'inferno', direction = 1, na.value = NA) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ COMMON_NAME)

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = log(LITTORAL_AREA_ACRES))) +
  scale_color_viridis(name = "log littoral area", option = 'inferno', direction = 1, na.value = NA) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ COMMON_NAME)


# smaller subset ----------------------------------------------------------

# select only northern pike
# fish_dat = dat %>%
#   filter(COMMON_NAME == 'northern pike') %>%
#   select(-c(SURVEY_ID, LKNAME, COMMON_NAME, EFFORT, TOTAL_CATCH, CTY:SURVEYTYPE, FISHERIES_WATERBODY_ID.x:ALT_LAKE_NAME)) %>%
#   mutate(SURVEYDATE = mdy(SURVEYDATE),
#          DOW = as.factor(DOW),
#          GEAR = as.factor(GEAR)) %>%
#   group_by(.dots=c('DOW', 'SURVEYDATE', 'GEAR')) %>%
#   summarise_all( ~ mean(., na.rm = T))

# select only northern pike
fish_dat = dat %>%
  filter(COMMON_NAME == 'northern pike') %>%
  select(-c(SURVEYDATE, SURVEY_ID, LKNAME, COMMON_NAME, EFFORT, TOTAL_CATCH, CTY:SURVEYTYPE, FISHERIES_WATERBODY_ID.x:ALT_LAKE_NAME)) %>%
  mutate(DOW = as.factor(DOW),
         GEAR = as.factor(GEAR)) %>%
  group_by(.dots=c('DOW', 'GEAR')) %>%
  summarise_all( ~ mean(., na.rm = T))
  



# plots by longitude and latitude -----------------------------------------

lats = range(fish_dat$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(fish_dat$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = log(CPUE))) +
  scale_color_viridis(name = "Log CPUE", option = 'inferno', direction = 1) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ GEAR)


ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = secchi.m)) +
  scale_color_viridis(name = "Secchi", option = 'inferno', direction = 1, na.value = NA) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~ GEAR)




# density plots -----------------------------------------------------------

ggplot(data = fish_dat, aes(x = log(CPUE), fill = GEAR)) +
  geom_density(alpha = 0.5) +
  ggtitle('Fish Density by log CPUE') +
  xlab('Log CPUE')

ggplot(data = fish_dat, aes(x = log(CPUE), fill = GEAR)) +
  geom_density(alpha = 0.5) +
  ggtitle('Fish Density by log CPUE') +
  xlab('Log CPUE')

ggplot(data = fish_dat, aes(x = secchi.m, fill = GEAR)) +
  geom_density(alpha = 0.5) +
  ggtitle('Fish Density by Secchi Depth') +
  xlab('Secchi Depth')




