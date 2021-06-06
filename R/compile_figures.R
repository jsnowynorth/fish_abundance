##----------------------------------------------------
## Name: Joshua North
##
## Date: 03/01/2020
##
## Project: Fish Abundance
##
## Objective: Multivariate fish abundance model. Land use covariates, GDD, secchi.
##            Separate function for each component.
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


# load in data ------------------------------------------------------------

static = read_csv('data/Static_MN_Data_JE.csv') %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
effort = read_csv('data/MN_fish_catch_effort_all_lakes_years_surveys_JE.csv') %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
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
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME),
         GEAR = as.factor(GEAR),
         SURVEYDATE = mdy(SURVEYDATE)) %>%
  filter(complete.cases(.)) %>% 
  filter(COMMON_NAME != 'white sucker',
         COMMON_NAME != 'smallmouth bass') %>%
  filter(SURVEYDATE >= '1993-01-01') %>% 
  filter(GEAR == 'GN' | GEAR == 'TN') %>%
  mutate(COMMON_NAME = droplevels(COMMON_NAME),
         DOW = droplevels(DOW),
         GEAR = droplevels(GEAR)) %>% 
  arrange(DOW) %>% 
  mutate(year = year(SURVEYDATE))

# add in land use
fish_dat = fish_dat %>% 
  left_join(land %>% select(-year), by = c('DOW')) %>% 
  filter(complete.cases(.))

# all <- fish_dat %>% 
#   group_by(SURVEYDATE, DOW) %>% 
#   tidyr::expand(COMMON_NAME, SURVEYDATE, GEAR) %>% 
#   ungroup() %>% 
#   arrange(DOW)

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
  mutate(DOY = yday(SURVEYDATE),
         DOY_sin_semi = sin(DOY/365 * 4*pi),
         DOY_cos_semi = cos(DOY/365 * 4*pi),
         DOY_sin_semi_temp = DOY_sin_semi * temp_0,
         DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>% 
  mutate(DOW = as.factor(DOW)) %>% 
  arrange(DOW)


colnames(fish_dat)
mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:19, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9, 13:19)]
catch_covs = colnames(fish_dat)[c(24, 27, 28)] # no temp doy interaction
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]
# need columns: DOW, COMMON_NAME, EFFORT, SURVEYDATE

# model -------------------------------------------------------------------

source('R/lewis_model.R')

mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:15,17, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]

# load in sampler data ----------------------------------------------------

run_1 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/full_model_test_1.rds')
run_2 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/full_model_test_2.rds')
run_1 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/full_model.rds')
run_2 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/full_model_2.rds')
run_3 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/full_model_3.rds')
run_4 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/full_model_4.rds')
# run_5 = read_rds(file = '/Users/joshuanorth/Desktop/lewis_results/full_model_5.rds')

beta = run$beta

phi = run$phi

omega = run$omega

sigma_species = run$sigma_species


beta = abind(run_1$beta, 
             run_2$beta, along = 3)

phi = abind(run_1$phi, 
            run_2$phi, along = 3)

omega = abind(run_1$omega, 
              run_2$omega, along = 2)

sigma_species = abind(run_1$sigma_species, 
                      run_2$sigma_species, along = 3)


beta = abind(run_1$beta, 
             run_2$beta, 
             run_3$beta, 
             run_4$beta, along = 3)

phi = abind(run_1$phi, 
            run_2$phi, 
            run_3$phi, 
            run_4$phi, along = 3)

omega = abind(run_1$omega, 
              run_2$omega, 
              run_3$omega, 
              run_4$omega, along = 2)

sigma_species = abind(run_1$sigma_species, 
                      run_2$sigma_species, 
                      run_3$sigma_species, 
                      run_4$sigma_species, along = 3)

rm(run_1, run_2, run_3, run_4)

pars <- create_pars(fish_dat, mean_covs, mean_covs_log, mean_covs_logit, catch_covs)
b_names = colnames(pars$X[[1]])
phi_names = colnames(pars$Z[[1]])

par(mfrow = c(3,3))
for(i in 1:9){
  plot(beta[6,i,], type = 'l', main = b_names[i])
  abline(h = 0)
}


par(mfrow = c(3,4))
for(i in 1:11){
  plot(phi[1,i,], type = 'l', main = phi_names[i])
}


par(mfrow = c(3,6))
for(i in 1001:1018){
  j = i %% 6
  if(j == 0){j=6}
  plot(omega[i,], type = 'l', main = pars$fish_names[j])
}

par(mfrow = c(2,3))
for(i in 1:6){
  plot(sigma_species[3,i,], type = 'l', main = pars$fish_names[i])
}


# xtable output -----------------------------------------------------------

fnames = pars$fish_names %>% 
  str_to_title()

b_names = b_names %>% 
  str_replace_all("_", " ") %>% 
  str_to_title() %>% 
  str_replace_all('Dd5', 'DD5') %>% 
  str_replace_all('Gis', 'GIS')

phi_names = phi_names %>% 
  str_replace_all("_", " ") %>% 
  str_to_title() %>% 
  str_replace_all('Gn', 'GN')


# relative abundance

b_mean = apply(beta[,,-c(1:2000)], c(1,2), mean)
b_lower = apply(beta[,,-c(1:2000)], c(1,2), quantile, probs = 0.025)
b_upper = apply(beta[,,-c(1:2000)], c(1,2), quantile, probs = 0.975)

b_sig = array(0, dim = dim(b_mean))
b_sig = ifelse((b_upper > 0) & (b_lower < 0), 1, 0)

colnames(b_mean) = b_names
rownames(b_mean) = fnames

b_mean_sig = round(b_mean, 3)
cint = round((b_upper - b_lower)/2, 3)

b_mean_sig = paste0(b_mean_sig, " (", cint, ")")
b_col = ifelse(b_sig, paste0("\\cellcolor{red}", b_mean_sig), paste0("\\cellcolor{white}", b_mean_sig))
colnames(b_col) = b_names
rownames(b_col) = fnames
print(xtable(t(b_col)), sanitize.text.function = identity)


# catchability

phi_mean = apply(phi[,,-c(1:2000)], c(1,2), mean)
phi_lower = apply(phi[,,-c(1:2000)], c(1,2), quantile, probs = 0.025)
phi_upper = apply(phi[,,-c(1:2000)], c(1,2), quantile, probs = 0.975)

p_sig = array(0, dim = dim(phi_mean))
p_sig = ifelse((phi_upper > 0) & (phi_lower < 0), 1, 0)

p_mean_sig = round(phi_mean, 3)
cint = round((phi_upper - phi_lower)/2, 3)

p_mean_sig = paste0(p_mean_sig, " (", cint, ")")
p_col = ifelse(p_sig, paste0("\\cellcolor{red}", p_mean_sig), paste0("\\cellcolor{white}", p_mean_sig))
colnames(p_col) = phi_names
rownames(p_col) = fnames
print(xtable(t(p_col)), sanitize.text.function = identity)


# covariance figure -------------------------------------------------------

cmat = cov2cor(apply(sigma_species[,,-c(1:2000)], c(1,2), mean))

fnames = c('Crappie', 'Bluegill', 'Bass', 'Pike', 'Walleye', 'Perch')

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

fnames = tibble(name = c('Crappie', 'Bluegill', 'Bass', 'Pike', 'Walleye', 'Perch'), ind = seq(1:6))

mean_complete <- upper_mean %>% 
  left_join(lower_mean, by = c('Row', 'Col')) %>% 
  mutate(Cov2 = sprintf( "%0.2f", Cov2)) %>% 
  rowwise() %>% 
  mutate(fish = fnames$name[which(Row[1] == fnames$ind)]) %>% 
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

# ggsave('results/spatial_results/species_dependenc.png', width = 10, height = 10)


# relative abundance ------------------------------------------------------

b_hat = apply(beta[,,-c(1:2000)], c(1,2), mean)
colnames(b_hat) = b_names

phi_hat = apply(phi[,,-c(1:2000)], c(1,2), mean)
colnames(phi_hat) = phi_names

# construct omega
K = length(pars$fish_names)
omega_hat = matrix(apply(omega[,-c(1:2000)], 1, mean), ncol = K, byrow = F)
ind_array = tibble(id = pars$lake_id, as_tibble(omega_hat, .name_repair = ~LETTERS[1:K]))
lake_array = tibble(id = pars$lake_index)
omega_hat = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))


lam_hat = list()
for(k in 1:pars$K){
  lam_hat[[k]] = exp(pars$X[[k]] %*% b_hat[k,] + omega_hat[k] + pars$Z[[k]] %*% phi_hat[k,])
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

# ggsave('results/spatial_results/relative_abund.png', width = 12, height = 12)

p1 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = rel_abun %>% filter(Fish == 'black_crappie'), 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple') +
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
  geom_point(data = rel_abun %>% filter(Fish == 'bluegill'), 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple') +
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
  geom_point(data = rel_abun %>% filter(Fish == 'largemouth_bass'), 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple') +
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
  geom_point(data = rel_abun %>% filter(Fish == 'northern_pike'), 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple') +
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
  geom_point(data = rel_abun %>% filter(Fish == 'walleye'), 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple') +
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
  geom_point(data = rel_abun %>% filter(Fish == 'yellow_perch'), 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abundance), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple') +
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


rel_abun_grid = cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2)

# ggsave('results/spatial_results/relative_abund.png', rel_abun_grid, width = 15, height = 10)



# cpue vs cpue hat --------------------------------------------------------

rel_abun = rbind(fish_dat %>% filter(COMMON_NAME == 'black crappie') %>% mutate(Abundance = c(lam_hat[[1]])),
                 fish_dat %>% filter(COMMON_NAME == 'bluegill') %>% mutate(Abundance = c(lam_hat[[2]])),
                 fish_dat %>% filter(COMMON_NAME == 'largemouth bass') %>% mutate(Abundance = c(lam_hat[[3]])),
                 fish_dat %>% filter(COMMON_NAME == 'northern pike') %>% mutate(Abundance = c(lam_hat[[4]])),
                 fish_dat %>% filter(COMMON_NAME == 'walleye') %>% mutate(Abundance = c(lam_hat[[5]])),
                 fish_dat %>% filter(COMMON_NAME == 'yellow perch') %>% mutate(Abundance = c(lam_hat[[6]])))

rel_abun = rel_abun %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW')

# CPUE by Gear type
cpue_rel_GN = rel_abun %>% 
  filter(GN == 1) %>% 
  filter(CPUE < quantile(CPUE, probs = 0.99)) %>% 
  mutate(yr = year(SURVEYDATE)) %>% 
  select(COMMON_NAME, Abundance, CPUE, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
  group_by(COMMON_NAME) %>% 
  mutate_at(vars(CPUE, Abundance), ~ (. - mean(.))/sd(.)) %>% 
  ungroup()

cpue_rel_TN = rel_abun %>% 
  filter(TN == 1) %>% 
  filter(CPUE < quantile(CPUE, probs = 0.99)) %>% 
  mutate(yr = year(SURVEYDATE)) %>% 
  select(COMMON_NAME, Abundance, CPUE, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
  group_by(COMMON_NAME) %>% 
  mutate_at(vars(CPUE, Abundance), ~ (. - mean(.))/sd(.)) %>% 
  ungroup()

ggplot(cpue_rel_GN, aes(x = CPUE, y = Abundance)) +
  geom_point(size = 3) +
  stat_cor(aes(label = ..r.label..), label.x.npc = 'center', label.y.npc = 'top', size = 10) +
  facet_wrap(~COMMON_NAME, 
             scales = 'free', 
             labeller = labeller(COMMON_NAME = c('black crappie' = 'Black Crappie',
                                                 'bluegill' = 'Bluegill',
                                                 'largemouth bass' = 'Largemouth Bass',
                                                 'northern pike' = 'Northern Pike',
                                                 'walleye' = 'Walleye',
                                                 'yellow perch' = 'Yellow Perch'))) +
  ylab('Relative Abuncance') +
  xlab('Catch Per Unit Effort') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/cpue_rel_GN.png', width = 12, height = 10)


ggplot(cpue_rel_TN, aes(x = CPUE, y = Abundance)) +
  geom_point(size = 3) +
  stat_cor(aes(label = ..r.label..), label.x.npc = 'center', label.y.npc = 'top', size = 10) +
  facet_wrap(~COMMON_NAME, 
             scales = 'free', 
             labeller = labeller(COMMON_NAME = c('black crappie' = 'Black Crappie',
                                                 'bluegill' = 'Bluegill',
                                                 'largemouth bass' = 'Largemouth Bass',
                                                 'northern pike' = 'Northern Pike',
                                                 'walleye' = 'Walleye',
                                                 'yellow perch' = 'Yellow Perch'))) +
  ylab('Relative Abuncance') +
  xlab('Catch Per Unit Effort') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/cpue_rel_TN.png', width = 12, height = 10)


# spatial residual plot ---------------------------------------------------

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
  ggtitle("Relative Log Abundance") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# relative abundance by CPUE ----------------------------------------------


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

# CPUE by Gear type
cpue_rel_GN = rel_abun %>% 
  filter(GN == 1) %>% 
  filter(CPUE < 80) %>%
  mutate(yr = year(SURVEYDATE)) %>% 
  select(COMMON_NAME, Abundance, CPUE, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5)

cpue_rel_TN = rel_abun %>% 
  filter(TN == 1) %>% 
  filter(CPUE < 80) %>%
  mutate(yr = year(SURVEYDATE)) %>% 
  select(COMMON_NAME, Abundance, CPUE, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5)

ggplot(cpue_rel_GN, aes(x = CPUE, y = Abundance)) +
  geom_point(size = 3) +
  stat_cor(aes(label = ..r.label..), label.x.npc = 'center', label.y.npc = 'top', size = 10) +
  facet_wrap(~COMMON_NAME, 
             scales = 'free', 
             labeller = labeller(COMMON_NAME = c('black crappie' = 'Black Crappie',
                                                 'bluegill' = 'Bluegill',
                                                 'largemouth bass' = 'Largemouth Bass',
                                                 'northern pike' = 'Northern Pike',
                                                 'walleye' = 'Walleye',
                                                 'yellow perch' = 'Yellow Perch'))) +
  ylab('Relative Abuncance') +
  xlab('Catch Per Unit Effort') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/cpue_rel_GN.png', width = 12, height = 10)


ggplot(cpue_rel_TN, aes(x = CPUE, y = Abundance)) +
  geom_point(size = 3) +
  stat_cor(aes(label = ..r.label..), label.x.npc = 'center', label.y.npc = 'top', size = 10) +
  facet_wrap(~COMMON_NAME, 
             scales = 'free', 
             labeller = labeller(COMMON_NAME = c('black crappie' = 'Black Crappie',
                                                 'bluegill' = 'Bluegill',
                                                 'largemouth bass' = 'Largemouth Bass',
                                                 'northern pike' = 'Northern Pike',
                                                 'walleye' = 'Walleye',
                                                 'yellow perch' = 'Yellow Perch'))) +
  ylab('Relative Abuncance') +
  xlab('Catch Per Unit Effort') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/cpue_rel_TN.png', width = 12, height = 10)


# cpue spatial comparison -------------------------------------------------

cpue_rel_spat = rel_abun %>% 
  filter(CPUE < 80) %>%
  mutate(yr = year(SURVEYDATE)) %>% 
  group_by(DOW, COMMON_NAME, yr) %>% 
  mutate(CPUE = sum(CPUE, na.rm = T)) %>% 
  distinct(DOW, .keep_all = T) %>% 
  ungroup() %>% 
  select(COMMON_NAME, Abundance, CPUE, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
  pivot_longer(c(Abundance, CPUE), names_to = 'Est', values_to = 'Abun')



p1 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'black crappie'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Black Crappie - CPUE") +
  guides(color = F)

p2 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'black crappie'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Black Crappie - Abundance") +
  guides(color = F)

p3 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'bluegill'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Bluegill - CPUE") +
  guides(color = F)

p4 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'bluegill'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Bluegill - Abundance") +
  guides(color = F)

p5 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'largemouth bass'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Largemouth Bass - CPUE") +
  guides(color = F)

p6 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'largemouth bass'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Largemouth Bass - Abundance") +
  guides(color = F)

p7 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'northern pike'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Northern Pike - CPUE") +
  guides(color = F)

p8 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'northern pike'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Northern Pike - Abundance") +
  guides(color = F)

p9 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'walleye'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Walleye - CPUE") +
  guides(color = F)

p10 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'walleye'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Walleye - Abundance") +
  guides(color = F)

p11 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'CPUE') %>% 
               filter(COMMON_NAME == 'yellow perch'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Yellow Perch - CPUE") +
  guides(color = F)

p12 = ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = cpue_rel_spat %>% 
               filter(Est == 'Abundance') %>% 
               filter(COMMON_NAME == 'yellow perch'), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 3) +
  scale_color_gradient(low = 'green', high = 'purple', na.value = 'white') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle("Yellow Perch - Abundance") +
  guides(color = F)

p1_grid = cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, byrow = F)
p2_grid = cowplot::plot_grid(p7, p8, p9, p10, p11, p12, nrow = 2, byrow = F)
# ggsave('results/spatial_results/cpue_rel_spat1.png',p1_grid, width = 15, height = 10)
# ggsave('results/spatial_results/cpue_rel_spat2.png',p2_grid, width = 15, height = 10)


# catchability - temperature ----------------------------------------

d_plot_lake = function(phis, fish_names, yr, fish_dat, lake_dow){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  tC = temp %>%
    filter(year(SURVEYDATE) == yr) %>%
    filter(DOW == lake_dow) %>% 
    filter(DOY != 366) %>% 
    select(-c(DOW, DD5, DOY))
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == yr) %>% 
    select(SURVEYDATE, all_of(catch_covs), GN) %>% 
    mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN))
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:TN) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == yr) %>% 
    pivot_longer(GN:TN, names_to = "Gear", values_to = "Ind") %>% 
    mutate(Gear = factor(Gear)) %>% 
    filter(Ind == 1) %>% 
    group_by(Gear, SURVEYDATE) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    ungroup() %>% 
    select(-EFFORT) %>% 
    spread(key = Gear, value = Ind, drop = F, fill = 0)
  
  
  TN_df = year_select %>% 
    filter(GN != 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 0, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>%
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  GN_df = year_select %>% 
    filter(GN == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 1, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>%
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  TN = as.matrix(TN_df)
  GN = as.matrix(GN_df)
  
  n_runs = dim(phis)[3]
  
  post_stores_TN = post_stores_GN = array(NA, dim = c(nrow(tC), length(fnames), n_runs))
  
  for(i in 1:n_runs){
    post_stores_TN[,,i] = TN %*% t(phis[,,i])
    post_stores_GN[,,i] = GN %*% t(phis[,,i])
  }
  
  TN = apply(post_stores_TN, c(1,2), mean)
  TN_lower = apply(post_stores_TN, c(1,2), quantile, probs = 0.025)
  TN_upper = apply(post_stores_TN, c(1,2), quantile, probs = 0.975)
  
  GN = apply(post_stores_GN, c(1,2), mean)
  GN_lower = apply(post_stores_GN, c(1,2), quantile, probs = 0.025)
  GN_upper = apply(post_stores_GN, c(1,2), quantile, probs = 0.975)
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  colnames(TN_lower) = paste0("TN_",fish_names)
  colnames(GN_lower) = paste0("GN_",fish_names)
  
  colnames(TN_upper) = paste0("TN_",fish_names)
  colnames(GN_upper) = paste0("GN_",fish_names)
  
  
  date_select = year_select %>% 
    distinct(SURVEYDATE) %>% 
    mutate(sample_day = SURVEYDATE) %>% 
    right_join(tC, by = c('SURVEYDATE')) %>% 
    left_join(gear_ind, by = c('SURVEYDATE')) %>% 
    select(-c(temp_0)) %>% 
    mutate_at(vars(GN:TN), ~ if_else(. == 1, SURVEYDATE, NaN)) %>% 
    select(-sample_day) %>% 
    arrange(SURVEYDATE) %>% 
    rename_at(vars(-SURVEYDATE), ~ paste0(., '_survey'))
  
  mean = tibble(as_tibble(TN), as_tibble(GN), date_select)
  lower = tibble(as_tibble(TN_lower), as_tibble(GN_lower), date_select)
  upper = tibble(as_tibble(TN_upper), as_tibble(GN_upper), date_select)
  
  return(list(mean = mean,
              lower = lower,
              upper = upper))
}

yr = 2016
fnames = c('crappie', 'bluegill', 'bass', 'pike', 'walleye', 'perch')
lake_dow_south = 53002800 # south
lake_dow_north = 16001900 # north
lake_dow_center = 18005000 # center
plt_dat_south = d_plot_lake(phi, fnames, yr, fish_dat, lake_dow_south)
plt_dat_north = d_plot_lake(phi, fnames, yr, fish_dat, lake_dow_north)
plt_dat_center = d_plot_lake(phi, fnames, yr, fish_dat, lake_dow_center)


plt_dat_south$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_South") %>% 
  left_join(plt_dat_north$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_North"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat_center$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_Cetner"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  pivot_longer(cols = -c(SURVEYDATE:Fish), names_to = c("Type", "Lake"), names_sep = "_", values_to = "Effectiveness") %>% 
  select(-Type) %>% 
  left_join(plt_dat_south$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-sample_day) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(aes(linetype = Lake), size = 0.9) +
  facet_wrap(~ Gear) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-04-01')), ymd(paste0(yr, '-10-01')))) +
  scale_linetype_manual(values=c("solid", 'longdash', "dotted")) +
  ggtitle(paste0('Gear Type Effectiveness by Region, ', yr)) +
  xlab('') +
  ylab('Relative Effectiveness') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ylim(c(-25, 25))

# ggsave(paste0('results/spatial_results/catchability_lake_comp_', yr, '.png'), width = 10, height = 6)


plt_dat = plt_dat_south
plt_dat = plt_dat_north
plt_dat = plt_dat_center

plt_dat$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness") %>% 
  left_join(plt_dat$upper %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_upper"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat$lower %>% select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_lower"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>%
  left_join(plt_dat$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-c(sample_day)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(size = 0.8) +
  facet_wrap(~ Gear) +
  geom_ribbon(aes(ymin = Effectiveness_lower, ymax = Effectiveness_upper), alpha = 0.2) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-04-01')), ymd(paste0(yr, '-10-01')))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ggtitle(paste0('Gear Type Effectiveness for ', yr)) +
  xlab('') +
  ylab('Relative Effectiveness') +
  ylim(c(-30, 30))


# ggsave(paste0('results/spatial_results/catchability_south_lake_', yr, '.png'), width = 10, height = 6)
# ggsave(paste0('results/spatial_results/catchability_north_lake_', yr, '.png'), width = 10, height = 6)
# ggsave(paste0('results/spatial_results/catchability_cental_lake_', yr, '.png'), width = 10, height = 6)


# catchability by temperature by gear -------------------------------------

### by species and gear
dft = plt_dat_south$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_South") %>% 
  left_join(plt_dat_north$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_North"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat_center$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Effectiveness_Cetner"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  pivot_longer(cols = -c(SURVEYDATE:Fish), names_to = c("Type", "Lake"), names_sep = "_", values_to = "Effectiveness") %>% 
  select(-Type) %>% 
  left_join(plt_dat_south$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-sample_day)

pa = ggplot(dft %>% filter(Fish == 'crappie' | Fish == 'bluegill' | Fish ==  'bass'), aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(aes(linetype = Lake), size = 1) +
  facet_wrap(~ Gear + Fish, ncol = 3,
             labeller = labeller(Fish = c('crappie' = 'Black Crappie',
                                          'bluegill' = 'Bluegill',
                                          'bass' = 'Largemouth Bass',
                                          'pike' = 'Northern Pike',
                                          'walleye' = 'Walleye',
                                          'perch' = 'Yellow Perch'))) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-04-01')), ymd(paste0(yr, '-10-01')))) +
  scale_linetype_manual(values=c("solid", 'longdash', "dotted")) +
  ggtitle(paste0('Gear Type Effectiveness by Region, ', yr)) +
  xlab('') +
  ylab('Relative Effectiveness') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size=14)) +
  ylim(c(-25, 25)) +
  guides(color = F)

# ggsave(paste0('results/spatial_results/catchability_lake_comp_species_a_', yr, '.png'), pa, width = 12, height = 8)


pb = ggplot(dft %>% filter(Fish != 'crappie' & Fish != 'bluegill' & Fish !=  'bass'), aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(aes(linetype = Lake), size = 1) +
  facet_wrap(~ Gear + Fish, ncol = 3,
             labeller = labeller(Fish = c('crappie' = 'Black Crappie',
                                          'bluegill' = 'Bluegill',
                                          'bass' = 'Largemouth Bass',
                                          'pike' = 'Northern Pike',
                                          'walleye' = 'Walleye',
                                          'perch' = 'Yellow Perch'))) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-04-01')), ymd(paste0(yr, '-10-01')))) +
  scale_linetype_manual(values=c("solid", 'longdash', "dotted")) +
  ggtitle(paste0('Gear Type Effectiveness by Region, ', yr)) +
  xlab('') +
  ylab('Relative Effectiveness') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size=14)) +
  ylim(c(-25, 25)) +
  guides(color = F)

# ggsave(paste0('results/spatial_results/catchability_lake_comp_species_b_', yr, '.png'), pb, width = 12, height = 8)


pfull = ggplot(dft, aes(x = SURVEYDATE, y = Effectiveness, color = Fish, fill = Fish)) +
  geom_line(aes(linetype = Lake), size = 1) +
  facet_wrap(~ Gear + Fish, ncol = 3,
             labeller = labeller(Fish = c('crappie' = 'Black Crappie',
                                          'bluegill' = 'Bluegill',
                                          'bass' = 'Largemouth Bass',
                                          'pike' = 'Northern Pike',
                                          'walleye' = 'Walleye',
                                          'perch' = 'Yellow Perch'))) +
  geom_rug(aes(x = date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr, '-04-01')), ymd(paste0(yr, '-10-01')))) +
  scale_linetype_manual(values=c("solid", 'longdash', "dotted")) +
  ggtitle(paste0('Gear Type Effectiveness by Region, ', yr)) +
  xlab('') +
  ylab('Relative Effectiveness') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size=14)) +
  ylim(c(-25, 25)) +
  guides(color = F)

ggsave(paste0('results/spatial_results/catchability_lake_comp_species_', yr, '.png'), pfull, width = 14, height = 16)


# difference in catchability cold/hot -------------------------------------

d_plot_lake = function(phis, fish_names, yr, fish_dat, lake_dow){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  tC = temp %>%
    filter(year(SURVEYDATE) == yr) %>%
    filter(DOW == lake_dow) %>% 
    filter(DOY != 366) %>% 
    select(-c(DOW, DD5, DOY))
  
  year_select = fish_dat %>% 
    filter(year(SURVEYDATE) == yr) %>% 
    select(SURVEYDATE, all_of(catch_covs), GN) %>% 
    mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN))
  
  gear_ind = fish_dat %>% 
    select(SURVEYDATE, EFFORT, GN:TN) %>% 
    filter(EFFORT != 0) %>% 
    filter(year(SURVEYDATE) == yr) %>% 
    pivot_longer(GN:TN, names_to = "Gear", values_to = "Ind") %>% 
    mutate(Gear = factor(Gear)) %>% 
    filter(Ind == 1) %>% 
    group_by(Gear, SURVEYDATE) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    ungroup() %>% 
    select(-EFFORT) %>% 
    spread(key = Gear, value = Ind, drop = F, fill = 0)
  
  
  TN_df = year_select %>% 
    filter(GN != 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 0, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>%
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  GN_df = year_select %>% 
    filter(GN == 1) %>% 
    distinct(SURVEYDATE, .keep_all = T) %>% 
    right_join(tC, by = c('SURVEYDATE', 'temp_0')) %>% 
    arrange(SURVEYDATE) %>% 
    mutate(GN = 1, 
           DOY = yday(SURVEYDATE),
           DOY_sin_semi = sin(DOY/365 * 4*pi),
           DOY_cos_semi = cos(DOY/365 * 4*pi),
           DOY_sin_semi_temp = DOY_sin_semi * temp_0,
           DOY_cos_semi_temp = DOY_cos_semi * temp_0) %>%
    mutate_at(vars(temp_0:DOY_cos_semi_temp), .funs = list(GN = ~.*GN)) %>% 
    select(-c(SURVEYDATE, DOY))
  
  TN = as.matrix(TN_df)
  GN = as.matrix(GN_df)
  
  n_runs = dim(phis)[3]
  
  post_stores_TN = post_stores_GN = array(NA, dim = c(nrow(tC), length(fnames), n_runs))
  
  for(i in 1:n_runs){
    post_stores_TN[,,i] = TN %*% t(phis[,,i])
    post_stores_GN[,,i] = GN %*% t(phis[,,i])
  }
  
  TN = apply(post_stores_TN, c(1,2), mean)
  TN_lower = apply(post_stores_TN, c(1,2), quantile, probs = 0.025)
  TN_upper = apply(post_stores_TN, c(1,2), quantile, probs = 0.975)
  
  GN = apply(post_stores_GN, c(1,2), mean)
  GN_lower = apply(post_stores_GN, c(1,2), quantile, probs = 0.025)
  GN_upper = apply(post_stores_GN, c(1,2), quantile, probs = 0.975)
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  colnames(TN_lower) = paste0("TN_",fish_names)
  colnames(GN_lower) = paste0("GN_",fish_names)
  
  colnames(TN_upper) = paste0("TN_",fish_names)
  colnames(GN_upper) = paste0("GN_",fish_names)
  
  
  date_select = year_select %>% 
    distinct(SURVEYDATE) %>% 
    mutate(sample_day = SURVEYDATE) %>% 
    right_join(tC, by = c('SURVEYDATE')) %>% 
    left_join(gear_ind, by = c('SURVEYDATE')) %>% 
    select(-c(temp_0)) %>% 
    mutate_at(vars(GN:TN), ~ if_else(. == 1, SURVEYDATE, NaN)) %>% 
    select(-sample_day) %>% 
    arrange(SURVEYDATE) %>% 
    rename_at(vars(-SURVEYDATE), ~ paste0(., '_survey'))
  
  mean = tibble(as_tibble(TN), as_tibble(GN), date_select)
  lower = tibble(as_tibble(TN_lower), as_tibble(GN_lower), date_select)
  upper = tibble(as_tibble(TN_upper), as_tibble(GN_upper), date_select)
  
  return(list(mean = mean,
              lower = lower,
              upper = upper))
}

yr_1 = 1993
yr_2 = 2016
fnames = c('crappie', 'bluegill', 'bass', 'pike', 'walleye', 'perch')
lake_dow_south = 53002800 # south
lake_dow_north = 16001900 # north
lake_dow_center = 18005000 # center
plt_dat_south_cold = d_plot_lake(phi, fnames, yr_1, fish_dat, lake_dow_south)
plt_dat_north_cold = d_plot_lake(phi, fnames, yr_1, fish_dat, lake_dow_north)
plt_dat_center_cold = d_plot_lake(phi, fnames, yr_1, fish_dat, lake_dow_center)
plt_dat_south_hot = d_plot_lake(phi, fnames, yr_2, fish_dat, lake_dow_south)
plt_dat_north_hot = d_plot_lake(phi, fnames, yr_2, fish_dat, lake_dow_north)
plt_dat_center_hot = d_plot_lake(phi, fnames, yr_2, fish_dat, lake_dow_center)


cold = plt_dat_south_cold$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "South") %>% 
  left_join(plt_dat_north_cold$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "North"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat_center_cold$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Center"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  pivot_longer(cols = -c(SURVEYDATE:Fish), names_to = c("Lake"), values_to = "Effectiveness") %>% 
  left_join(plt_dat_south_cold$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "Date"), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-sample_day) %>% 
  mutate(Year = yr_1)

hot = plt_dat_south_hot$mean %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "South") %>% 
  left_join(plt_dat_north_hot$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "North"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  left_join(plt_dat_center_hot$mean %>% 
              select(-c(GN_survey:TN_survey)) %>% 
              pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "Center"), by = c('SURVEYDATE', 'Gear', 'Fish')) %>% 
  pivot_longer(cols = -c(SURVEYDATE:Fish), names_to = c("Lake"), values_to = "Effectiveness") %>% 
  left_join(plt_dat_south_hot$mean %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = 'Date'), by = c('Gear', 'SURVEYDATE')) %>% 
  select(-sample_day) %>% 
  mutate(Year = yr_2)

cold = cold %>% 
  mutate(SURVEYDATE=ymd(format(SURVEYDATE, paste0(yr_2, "-%m-%d"))),
         Date=ymd(format(Date, paste0(yr_2, "-%m-%d"))))


cold %>% 
  rbind(hot) %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Year, fill = Year)) +
  geom_line(aes(linetype = Lake), size = 0.9) +
  facet_wrap(~ Fish + Gear) +
  geom_rug(aes(x = Date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr_2, '-04-01')), ymd(paste0(yr_2, '-10-01')))) +
  scale_linetype_manual(values=c("solid", 'longdash', "dotted")) +
  ggtitle(paste0('Gear Type Effectiveness by Region')) +
  xlab('') +
  ylab('Relative Effectiveness') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ylim(c(-25, 25))

# ggsave(paste0('results/spatial_results/catchability_lake_comp_species_', yr, '.png'), pfull, width = 14, height = 16)


cold %>% 
  rbind(hot) %>% 
  mutate(Year = factor(Year)) %>% 
  filter(Fish == 'perch') %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Year, fill = Year)) +
  geom_line(aes(linetype = Lake), size = 1) +
  facet_wrap(~ Gear) +
  geom_rug(aes(x = Date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr_2, '-04-01')), ymd(paste0(yr_2, '-10-01')))) +
  scale_linetype_manual(values = c("solid", 'longdash', "dotted")) +
  ggtitle(paste0('Gear Type Effectiveness by Region')) +
  xlab('') +
  ylab('Relative Effectiveness') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ylim(c(-25, 25))

# summaries ---------------------------------------------------------------



fish_dat %>% 
  select(ag:grass) %>% 
  pivot_longer(ag:grass, names_to = 'variable', values_to = 'value') %>% 
  group_by(variable) %>% 
  summarise_all(list(mean = mean,
                     lower = ~quantile(., probs = 0.025),
                     upper = ~quantile(., probs = 0.975)))

fish_dat %>% 
  select(ag:grass) %>% 
  pivot_longer(ag:grass, names_to = 'variable', values_to = 'value') %>% 
  mutate(value = value/100) %>% 
  mutate(value = log(value/(1-value))) %>% 
  ggplot(., aes(x = value)) +
  geom_histogram() +
  facet_wrap(~variable, scales = 'free_y')


fish_dat %>% 
  select(all_of(colnames(fish_dat)[c(7, 9)])) %>% 
  pivot_longer(all_of(colnames(fish_dat)[c(7, 9)]), names_to = 'variable', values_to = 'value') %>% 
  mutate(value = log(value)) %>% 
  ggplot(., aes(x = value)) +
  geom_histogram() +
  facet_wrap(~variable, scales = 'free_y')

