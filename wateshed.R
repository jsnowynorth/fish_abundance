library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(foreign)

# tmp = read_csv('/Users/joshuanorth/Desktop/fish_abundance/data/Static_MN_Data_JE.csv')
dat = read_csv('/Users/joshuanorth/Desktop/fish_abundance/data/mn_lake_list.csv') %>% 
  mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0"))) %>% 
  mutate(Bnum = str_sub(DOW,-2,-1))

fdat = read_csv('/Users/joshuanorth/Desktop/fish_abundance/data/fish_dat.csv')

fdat = fdat %>% 
  distinct(DOW, .keep_all = T) %>% 
  left_join(dat)

lk_bound = read_sf('/Users/joshuanorth/Desktop/hydrologic_units_WBDHU8_mn_3945305_01/hydrologic_units_wbdhu8_a_mn.shp')
st_geometry_type(lk_bound)
st_crs(lk_bound)
st_bbox(lk_bound)
s_inds = which(lk_bound$name == 'Snake')
lk_bound$name[s_inds[1]] = 'Snake_A'
lk_bound$name[s_inds[2]] = 'Snake_B'


tmp = fdat %>% 
  rename(x = LAKE_CENTER_LAT_DD5,
         y = LAKE_CENTER_LONG_DD5)

pnts_sf = st_as_sf(tmp, coords = c('y', 'x'), crs = 'NAD83')
# %>% st_transform(crs = 'NAD83')

fdat = fdat %>% 
  mutate(watershed = lk_bound$name[unlist(st_intersects(pnts_sf, lk_bound))])

# write_csv(fdat, 'data/watershed_id.csv')


lats = range(fdat$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(fdat$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot = FALSE))

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = fdat, 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = watershed)) +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(colour = "none")

length(table(fdat$watershed))


# st_centroid(lk_bound)$geometry[,1]
# lk_bound %>%
#   st_centroid() %>%
#   st_geometry() %>% 
#   stat_sf_coordinates()

water_centroids = tibble(watershed = lk_bound$name, long = st_coordinates(st_centroid(lk_bound$geometry))[,1], lat = st_coordinates(st_centroid(lk_bound$geometry))[,2])
# write_csv(water_centroids, 'data/watershed_centroids.csv')