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
library(readxl)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(fields)
library(viridis)



# load in data ------------------------------------------------------------


static = read_csv('data/Static_MN_Data_JE.csv')
static_del = read_csv('data/Static_MN_Data_JE_deleted_columns.csv')
dynamic = read_csv('data/Dynamic_MN_Data_JE.csv')
fish = read_csv('data/Coldwater_fish_PA_JE.csv')
effort = read_csv('data/MN_fish_catch_effort_all_lakes_years_surveys_JE.csv')



# join data ---------------------------------------------------------------


static = static %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
effort = effort %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
fish = fish %>% mutate(DOW = (str_pad(DOWNumber, 8, side = "left", pad = "0"))) %>% select(-DOWNumber)

dat = fish %>% left_join(static, by = c('DOW', 'LAKE_NAME')) #%>% left_join(effort, by = 'DOW')

# plots by longitude and latitude -----------------------------------------

lats = range(dat$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(dat$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = LAKE_AREA_DOW_ACRES)) +
  scale_color_viridis(name = "Area", limits = c(0, 2000), option = 'inferno', direction = -1) +
  xlab("Longitude") +
  ylab("Latitude")


fish_lon = dat %>% 
  select(c(cisco:laketrout, LAKE_CENTER_LONG_DD5, LAKE_CENTER_LAT_DD5, LAKE_AREA_DOW_ACRES)) %>% 
  pivot_longer(-c(LAKE_CENTER_LONG_DD5:LAKE_AREA_DOW_ACRES), values_to = 'Fish_ind', names_to = 'Fish') %>% 
  rename(lon = LAKE_CENTER_LONG_DD5, lat = LAKE_CENTER_LAT_DD5, area = LAKE_AREA_DOW_ACRES) %>% 
  mutate(lon = round(lon, 1), lat = round(lat, 1)) %>% 
  group_by(lon, Fish) %>% 
  summarise(Fish_ct = sum(Fish_ind)) %>% 
  ungroup()

fish_lat = dat %>% 
  select(c(cisco:laketrout, LAKE_CENTER_LONG_DD5, LAKE_CENTER_LAT_DD5, LAKE_AREA_DOW_ACRES)) %>% 
  pivot_longer(-c(LAKE_CENTER_LONG_DD5:LAKE_AREA_DOW_ACRES), values_to = 'Fish_ind', names_to = 'Fish') %>% 
  rename(lon = LAKE_CENTER_LONG_DD5, lat = LAKE_CENTER_LAT_DD5, area = LAKE_AREA_DOW_ACRES) %>% 
  mutate(lon = round(lon, 1), lat = round(lat, 1)) %>% 
  group_by(lat, Fish) %>% 
  summarise(Fish_ct = sum(Fish_ind)) %>% 
  ungroup()

table(round(fish_lat_lon$lon, 1))

p1 <- ggplot(fish_lon, aes(x = lon, y = Fish_ct, color = Fish)) +
  geom_point()

p2 <- ggplot(fish_lat, aes(x = lat, y = Fish_ct, color = Fish)) +
  geom_point()

cowplot::plot_grid(p1, p2)



# plots by secchi ---------------------------------------------------------

fish_trait = dat %>% 
  select(c(cisco:laketrout, secchi.m, LITTORAL_AREA_ACRES, littoral.zone, SHORE_LENGTH_MILES, LAKE_CENTER_LONG_DD5, LAKE_CENTER_LAT_DD5, LAKE_AREA_DOW_ACRES)) %>% 
  pivot_longer(c(cisco:laketrout), values_to = 'Fish_ind', names_to = 'Fish') %>% 
  rename(secchi = secchi.m, 
         littoral_area = LITTORAL_AREA_ACRES, 
         littoral_zone = littoral.zone,
         shore_len = SHORE_LENGTH_MILES,
         lon = LAKE_CENTER_LONG_DD5, 
         lat = LAKE_CENTER_LAT_DD5, 
         area = LAKE_AREA_DOW_ACRES) %>% 
  mutate(Fish_pres = ifelse(Fish_ind == 1, Fish, 'None')) %>% 
  filter(Fish_pres != 'None')


p1 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_trait, aes(x = lon, y = lat, color = secchi)) +
  scale_color_continuous(low = 'blue', high = 'darkred', guide = 'colorbar', na.value = NA, name = 'Secchi') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~Fish_pres, ncol = 1) +
  theme(legend.position = 'bottom')


p2 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_trait, aes(x = lon, y = lat, color = log(littoral_area))) +
  scale_color_continuous(low = 'blue', high = 'darkred', guide = 'colorbar', na.value = NA, name = 'Littoral Area') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~Fish_pres, ncol = 1) +
  theme(legend.position = 'bottom')

p3 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_trait, aes(x = lon, y = lat, color = littoral_zone)) +
  scale_color_continuous(low = 'blue', high = 'darkred', guide = 'colorbar', na.value = NA, name = 'Littoral Zone') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~Fish_pres, ncol = 1) +
  theme(legend.position = 'bottom')

p4 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fish_trait, aes(x = lon, y = lat, color = log(shore_len))) +
  scale_color_continuous(low = 'blue', high = 'darkred', guide = 'colorbar', na.value = NA, name = 'Shore Length') +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~Fish_pres, ncol = 1) +
  theme(legend.position = 'bottom')

cowplot::plot_grid(p1, p2, p3, p4, nrow = 1)


# fish_trait = dat %>% 
#   select(c(cisco:laketrout, secchi.m, LITTORAL_AREA_ACRES, littoral.zone, SHORE_LENGTH_MILES, LAKE_CENTER_LONG_DD5, LAKE_CENTER_LAT_DD5, LAKE_AREA_DOW_ACRES)) %>% 
#   pivot_longer(c(cisco:laketrout), values_to = 'Fish_ind', names_to = 'Fish') %>% 
#   rename(secchi = secchi.m, 
#          littoral_area = LITTORAL_AREA_ACRES, 
#          littoral_zone = littoral.zone,
#          shore_len = SHORE_LENGTH_MILES,
#          lon = LAKE_CENTER_LONG_DD5, 
#          lat = LAKE_CENTER_LAT_DD5, 
#          area = LAKE_AREA_DOW_ACRES) %>% 
#   mutate(Fish_pres = ifelse(Fish_ind == 1, Fish, 'None')) %>% 
#   filter(Fish_pres != 'None') %>% 
#   pivot_longer(c(secchi:shore_len), names_to = 'Metric', values_to = 'Val') %>% 
#   select(-c(Fish, Fish_ind))
# 
# 
# ggplot() +
#   geom_sf(data = usa) +
#   coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
#   geom_point(data = fish_trait, aes(x = lon, y = lat, color = Val)) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   facet_wrap(~ Fish_pres + Metric) +
#   scale_color_continuous(low = 'blue', high = 'darkred', guide = 'colorbar', na.value = 'black') 




p1 <- ggplot(data = fish_trait, aes(x = secchi, fill = Fish_pres)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(name = 'Fish', values = c('red', 'green', 'blue')) +
  ggtitle('Fish Density by Secchi Depth') +
  xlab('Secchi Depth')

p2 <- ggplot(data = fish_trait, aes(x = log(littoral_area), fill = Fish_pres)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(name = 'Fish', values = c('red', 'green', 'blue')) +
  ggtitle('Fish Density by Log Littoral Area') +
  xlab('Log Littoral Area')

p3 <- ggplot(data = fish_trait, aes(x = littoral_zone, fill = Fish_pres)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(name = 'Fish', values = c('red', 'green', 'blue')) +
  ggtitle('Fish Density by Littoral Zone') +
  xlab('Littoral Zone')

p4 <- ggplot(data = fish_trait, aes(x = log(shore_len), fill = Fish_pres)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(name = 'Fish', values = c('red', 'green', 'blue')) +
  ggtitle('Fish Density by Log Shore Length') +
  xlab('Log Shore Length')

cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)














