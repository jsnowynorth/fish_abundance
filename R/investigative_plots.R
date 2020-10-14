##----------------------------------------------------
## Name: Joshua North
##
## Date: 10/13/2020
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


# load data ---------------------------------------------------------------


cold = read_csv("data/coldwater_fish_PA.csv")
effort = read_csv("data/MN_fish_catch_effort_all_lakes_years_surveys.csv")
lake_list = read_csv("data/mn_lake_list.csv")
secchi = read_csv("data/Secchi_annual_CDOM_lakedepth.csv")
preds = read_csv("data/Static_lake_predictors.csv")
temp = read_excel("data/temperature_summarized_surface_lakeyear.xlsx")




# plots -------------------------------------------------------------------

lats = range(lake_list$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(lake_list$LAKE_CENTER_LONG_DD5, na.rm = T)

world <- ne_countries(scale = "small", returnclass = "sf", country = 'United States of America')
usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = lake_list, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = LAKE_AREA_MN_ACRES)) +
  scale_color_viridis(name = "Area", limits = c(0, 2000), option = 'inferno', direction = -1) +
  xlab("Longitude") +
  ylab("Latitude")



ggplot(lake_list, aes())

range(cold$DOWNumber)
range(lake_list$DOW, na.rm = T)
range(secchi$DOWLKNUM, na.rm = T)

secchi_sel = secchi %>% 
  select(DOWLKNUM, avgSecchi.m) %>% 
  rename(dow = DOWLKNUM)

cold_sel = cold %>% 
  select(DOWNumber, cisco, whitefish, laketrout) %>% 
  rename(dow = DOWNumber)

lake_sel = lake_list %>% 
  select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5, LAKE_AREA_MN_ACRES) %>% 
  rename(dow = DOW, lat = LAKE_CENTER_LAT_DD5, lon = LAKE_CENTER_LONG_DD5, acre = LAKE_AREA_MN_ACRES) %>% 
  mutate(dow = as.numeric(dow))

secchi_join = secchi_sel %>% 
  inner_join(cold_sel, by = 'dow')

lake_join = lake_sel %>% 
  inner_join(cold_sel, by = 'dow') %>% 
  pivot_longer(-c(dow, lat, lon, acre), names_to = 'fish', values_to = 'ind')

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = lake_join %>% filter(ind == 1), aes(x = lon, y = lat, color = fish)) +
  scale_color_discrete(name = "Fish") +
  xlab("Longitude") +
  ylab("Latitude")









