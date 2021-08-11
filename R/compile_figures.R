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
library(zoo)


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
  arrange(DOW, year, COMMON_NAME)

# all_mean = colnames(fish_dat)[c(7, 9, 23, 13:19, 25)]
# all_mean_log = colnames(fish_dat)[c(7, 9)]
# all_mean_logit = colnames(fish_dat)[c(13:19)]
# 
# p = fish_dat %>% 
#   select(DD5, temp_0) %>% 
#   ggpairs()

# model -------------------------------------------------------------------

fish_dat = read_csv('data/fish_dat.csv')

fish_dat = fish_dat %>% 
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME)) %>% 
  arrange(DOW, year, COMMON_NAME)

fish_dat = fish_dat %>%
  filter(year >= 2000) %>% 
  mutate(DOW = droplevels(DOW))

source('R/lewis_code/lewis_model.R')
# source('R/lewis_code/lewis_model_no_spatial.R')

# spatial covs
mean_covs = colnames(fish_dat)[c(7, 9, 23, 13:15,17, 25)]
temporal_covs = colnames(fish_dat)[c(23, 25)]
mean_covs_log = colnames(fish_dat)[c(7, 9)]
mean_covs_logit = colnames(fish_dat)[c(13:15,17)]
catch_covs = colnames(fish_dat)[c(24, 27:30)] # temp doy interaction
gear_types = colnames(fish_dat)[c(21, 22)]


# load in sampler data ----------------------------------------------------

# par(mfrow = c(3,3))
# for(i in 1:9){
#   plot(beta[1,i,], type = 'l')
#   abline(h = 0)
# }
# 
# par(mfrow = c(3,4))
# for(i in 1:11){
#   plot(phi[6,i,], type = 'l')
# }

run = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_mean_5_tmp.rds')
# run = read_rds(file = '/Users/joshuanorth/Desktop/full_model_non_spatial_3.rds')
beta = run$beta
phi = run$phi
omega = run$omega
sigma_species = run$sigma_species


# run_1 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_mean_1_tmp.rds')
# run_2 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_mean_2_tmp.rds')
# run_3 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_mean_3_tmp.rds')
# run_4 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_mean_4_tmp.rds')
# run_5 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_mean_5_tmp.rds')


run_1 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_mean_3_tmp.rds')
run_2 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_mean_4_tmp.rds')
run_3 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_mean_5_tmp.rds')

# run_1 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_1_tmp.rds')
# run_2 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_2_tmp.rds')
# run_3 = read_rds(file = '/Users/joshuanorth/Desktop/full_model_spatial_3_tmp.rds')
beta = abind(run_1$beta, 
             run_2$beta,
             run_3$beta, along = 3)

phi = abind(run_1$phi, 
            run_2$phi,
            run_3$phi, along = 3)

omega = abind(run_1$omega, 
              run_2$omega,
              run_3$omega, along = 3)

sigma_species = abind(run_1$sigma_species, 
                      run_2$sigma_species,
                      run_3$sigma_species,along = 3)

rm(run_1, run_2, run_3)



beta = abind(run_1$beta, 
             run_2$beta,
             run_3$beta,
             run_4$beta,
             run_5$beta, along = 3)

phi = abind(run_1$phi, 
            run_2$phi,
            run_3$phi,
            run_4$phi,
            run_5$phi, along = 3)

omega = abind(run_1$omega, 
              run_2$omega, 
              run_3$omega,
              run_4$omega,
              run_5$omega, along = 2)

sigma_species = abind(run_1$sigma_species, 
                      run_2$sigma_species,
                      run_3$sigma_species,
                      run_4$sigma_species,
                      run_5$sigma_species, along = 3)

rm(run_1, run_2, run_3, run_4, run_5)

# pars <- create_pars(fish_dat, mean_covs, mean_covs_log, mean_covs_logit, catch_covs) # non-spatial
pars <- create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs) # spatial
b_names = colnames(pars$X[[1]])
phi_names = colnames(pars$Z[[1]])

b_hat = apply(beta, c(1,2), mean)

round(b_hat, 3) %>%
  as_tibble(.name_repair = ~b_names) %>%
  mutate(Species = pars$fish_names) %>%
  relocate(Species)

# round(b_hat, 3) %>% 
#   as_tibble(.name_repair = ~b_names) %>% 
#   mutate(Species = pars$fish_names) %>% 
#   relocate(Species) %>% 
#   write_csv("/Users/joshuanorth/Desktop/b_hat_spatial.csv")
# 
# round(b_hat, 3) %>% 
#   as_tibble(.name_repair = ~b_names) %>% 
#   mutate(Species = pars$fish_names) %>% 
#   relocate(Species) %>% 
#   write_csv("/Users/joshuanorth/Desktop/b_hat_no_spatial.csv")
# 
# 
# round(b_hat, 3) %>% 
#   as_tibble(.name_repair = ~b_names) %>% 
#   mutate(Species = pars$fish_names) %>% 
#   relocate(Species) %>% 
#   write_csv("/Users/joshuanorth/Desktop/b_hat_no_land.csv")

fnames = pars$fish_names %>% str_replace(., ' ', '_')

for(j in 1:6){
  png(paste0('/Users/joshuanorth/Desktop/mean_chains/', fnames[j], '.png'), width = 800, height = 600)
  # png(paste0('/Users/joshuanorth/Desktop/full_chains/', fnames[j], '.png'), width = 800, height = 600)
  par(mfrow = c(3,3))
  for(i in 1:9){
    plot(beta[j,i,], type = 'l', main = b_names[i])
    abline(h = 0)
  }
  dev.off()
}


par(mfrow = c(3,3))
for(i in 1:9){
  plot(beta[6,i,], type = 'l', main = b_names[i])
  abline(h = 0)
}


par(mfrow = c(3,4))
for(i in 1:11){
  plot(phi[6,i,], type = 'l', main = phi_names[i])
}


par(mfrow = c(3,6))
for(i in 1001:1018){
  j = i %% 6
  if(j == 0){j=6}
  plot(omega[i,], type = 'l', main = pars$fish_names[j])
}

par(mfrow = c(3,6))
for(i in 101:103){
  plot(omega[i,1,], type = 'l', main = pars$fish_names[1])
  plot(omega[i,2,], type = 'l', main = pars$fish_names[2])
  plot(omega[i,3,], type = 'l', main = pars$fish_names[3])
  plot(omega[i,4,], type = 'l', main = pars$fish_names[4])
  plot(omega[i,5,], type = 'l', main = pars$fish_names[5])
  plot(omega[i,6,], type = 'l', main = pars$fish_names[6])
}


par(mfrow = c(2,3))
for(i in 1:6){
  plot(sigma_species[5,i,], type = 'l', main = pars$fish_names[i])
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

b_mean = apply(beta, c(1,2), mean)
b_lower = apply(beta, c(1,2), quantile, probs = 0.025)
b_upper = apply(beta, c(1,2), quantile, probs = 0.975)

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

phi_mean = apply(phi, c(1,2), mean)
phi_lower = apply(phi, c(1,2), quantile, probs = 0.025)
phi_upper = apply(phi, c(1,2), quantile, probs = 0.975)

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

cmat = cov2cor(apply(sigma_species, c(1,2), mean))

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


# sample date spatial plot ---------------------------------------------------

sdate = fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  mutate(day=yday(SURVEYDATE)) %>% 
  ungroup() %>% 
  group_by(DOW) %>% 
  summarize(min = min(day), max = max(day), mean = mean(day), median = median(day), sd = sd(day)) %>% 
  mutate(sd = ifelse(is.na(sd), 0, sd))

sdate = fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  mutate(day=yday(SURVEYDATE)) %>% 
  ungroup() %>% 
  group_by(DOW) %>% 
  summarize(min = min(SURVEYDATE), max = max(SURVEYDATE), mean = mean(SURVEYDATE), median = median(SURVEYDATE), sd = sd(day)) %>% 
  mutate(sd = ifelse(is.na(sd), 0, sd))

sdate_pat = sdate %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW')

lats = range(sdate_pat$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(sdate_pat$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))



as.Date(quantile(unclass(sdate_pat$median), seq(0.05, 0.95, 0.1)), origin = "1970-01-01")

ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = sdate_pat, 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = median), size = 3) +
  scale_color_fermenter(breaks = round(unname(quantile((sdate_pat)$median, probs = seq(0.05, 0.95, 0.15))),2), 
                        palette = 'YlOrRd', direction = 1, name = 'Day of Year') +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Median Sample Date") +
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
  geom_point(data = sdate_pat, 
             aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = sd), size = 3) +
  scale_color_fermenter(breaks = round(unname(quantile((sdate_pat)$sd, probs = seq(0.05, 0.95, 0.15))),2), palette = 'YlOrRd', direction = 1) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Standard Deviation in Sample Date") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/sd_sample_date_spat.png', width = 12, height = 8)


fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  mutate(day=yday(SURVEYDATE),
         SURVEYDATE=ymd(format(SURVEYDATE, "2016-%m-%d"))) %>% 
  ungroup() %>% 
  group_by(DOW) %>% 
  summarize(min = min(SURVEYDATE), max = max(SURVEYDATE), mean = mean(SURVEYDATE), median = median(SURVEYDATE), sd = sd(day)) %>% 
  mutate(sd = ifelse(is.na(sd), 0, sd)) %>% 
  ggplot(., aes(x = median)) +
  geom_histogram(bins = 50, color = 'black') +
  scale_x_date(date_breaks = 'week', date_labels = "%b %d", limits = c(ymd('2016-05-30'), ymd('2016-09-12'))) +
  ylab('Count') +
  xlab('') +
  ggtitle('Median Sample Date by Lake') +
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))


fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  mutate(day=yday(SURVEYDATE),
         SURVEYDATE=ymd(format(SURVEYDATE, "2016-%m-%d"))) %>% 
  ungroup() %>% 
  ggplot(., aes(x = SURVEYDATE)) +
  geom_histogram(bins = 50, color = 'black') +
  scale_x_date(date_breaks = 'week', date_labels = "%b %d", limits = c(ymd('2016-05-30'), ymd('2016-09-12'))) +
  ylab('Count') +
  xlab('') +
  ggtitle('Sample Date by Lake') +
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave('results/spatial_results/sample_date_hist.png', width = 12, height = 8)


# relative abundance ------------------------------------------------------

b_hat = apply(beta, c(1,2), mean)
colnames(b_hat) = b_names

phi_hat = apply(phi, c(1,2), mean)
colnames(phi_hat) = phi_names

# construct omega
K = length(pars$fish_names)
# omega_hat = pars$P %*% matrix(apply(omega, 1, mean), ncol = K, byrow = F)
omega_hat = pars$P %*% apply(omega, c(1,2), mean)
ind_array = tibble(id = pars$lake_id, as_tibble(omega_hat, .name_repair = ~LETTERS[1:K]))
lake_array = tibble(id = pars$lake_index)
omega_hat = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))


lam_hat = list()
for(k in 1:pars$K){
  lam_hat[[k]] = exp(pars$X[[k]] %*% b_hat[k,] + omega_hat[,k] + pars$Z[[k]] %*% phi_hat[k,])
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


# construct omega
K = length(pars$fish_names)
# omega_hat = pars$P %*% matrix(apply(omega, 1, mean), ncol = K, byrow = F)
omega_hat = pars$P %*% apply(omega, c(1,2), mean)
ind_array = tibble(id = pars$lake_id, as_tibble(omega_hat, .name_repair = ~LETTERS[1:K]))
lake_array = tibble(id = pars$lake_index)
omega_hat = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))
# omega_hat2 = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))
# omega_hat = omega_hat - omega_hat2

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

# ggsave('results/spatial_results/spatial_random_effect.png', width = 12, height = 10)

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


# cpue spatial comparison centered ----------------------------------------

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

cpue_rel_spat = rel_abun %>% 
  filter(CPUE < 80) %>%
  mutate(yr = year(SURVEYDATE)) %>% 
  group_by(DOW, COMMON_NAME, yr) %>% 
  mutate(CPUE = sum(CPUE, na.rm = T)) %>% 
  distinct(DOW, .keep_all = T) %>% 
  ungroup() %>% 
  select(COMMON_NAME, Abundance, CPUE, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
  pivot_longer(c(Abundance, CPUE), names_to = 'Est', values_to = 'Abun')

bc_no_cpue = bc %>%
  filter(CPUE == 0)

create_spatial_abun <- function(fish, abun, title, save_fig = FALSE){
  
  p = ggplot() +
    geom_sf(data = usa) +
    coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
    geom_point(data = cpue_rel_spat %>% 
                 filter(Est == abun) %>% 
                 filter(COMMON_NAME == fish), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 1.5) +
    # scale_color_gradientn(colors = c('red', 'white', 'blue')) +
    scale_color_fermenter(breaks = round(unname(quantile((cpue_rel_spat %>% 
                                                            filter(Est == abun) %>% 
                                                            filter(COMMON_NAME == fish))$Abun, probs = seq(0.05, 0.95, 0.1))),2), palette = 'RdBu') + #
    xlab("Longitude") +
    ylab("Latitude") +
    theme(legend.title=element_blank(),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = 'bold'),
          title = element_text(size = 16, face = 'bold')) +
    ggtitle(title)
  
  if(save_fig == TRUE){
    path = paste0('results/spatial_results/', str_replace_all(fish, ' ', '_'), '_', abun, '.png')
    ggsave(path, p, width = 8, height = 8)
  }else{
    return(p)
  }
  
}

create_spatial_abun('black crappie', 'CPUE', 'Black Crappie - CPUE', save_fig = F)
create_spatial_abun('black crappie', 'Abundance', 'Black Crappie - Abundance', save_fig = F)

create_spatial_abun('black crappie', 'CPUE', 'Black Crappie - CPUE', save_fig = T)
create_spatial_abun('black crappie', 'Abundance', 'Black Crappie - Abundance', save_fig = T)
create_spatial_abun('bluegill', 'CPUE', 'Bluegill - CPUE', save_fig = T)
create_spatial_abun('bluegill', 'Abundance', 'Bluegill - Abundance', save_fig = T)
create_spatial_abun('largemouth bass', 'CPUE', 'Largemouth Bass - CPUE', save_fig = T)
create_spatial_abun('largemouth bass', 'Abundance', 'Largemouth Bass - Abundance', save_fig = T)
create_spatial_abun('northern pike', 'CPUE', 'Northern Pike - CPUE', save_fig = T)
create_spatial_abun('northern pike', 'Abundance', 'Northern Pike - Abundance', save_fig = T)
create_spatial_abun('walleye', 'CPUE', 'Walleye - CPUE', save_fig = T)
create_spatial_abun('walleye', 'Abundance', 'Walleye - Abundance', save_fig = T)
create_spatial_abun('yellow perch', 'CPUE', 'Yellow Perch - CPUE', save_fig = T)
create_spatial_abun('yellow perch', 'Abundance', 'Yellow Perch - Abundance', save_fig = T)


p1 = create_spatial_abun('black crappie', 'CPUE', 'Black Crappie - CPUE', save_fig = F)
p2 = create_spatial_abun('black crappie', 'Abundance', 'Black Crappie - Abundance', save_fig = F)
p3 = create_spatial_abun('bluegill', 'CPUE', 'Bluegill - CPUE', save_fig = F)
p4 = create_spatial_abun('bluegill', 'Abundance', 'Bluegill - Abundance', save_fig = F)
p5 = create_spatial_abun('largemouth bass', 'CPUE', 'Largemouth Bass - CPUE', save_fig = F)
p6 = create_spatial_abun('largemouth bass', 'Abundance', 'Largemouth Bass - Abundance', save_fig = F)
p7 = create_spatial_abun('northern pike', 'CPUE', 'Northern Pike - CPUE', save_fig = F)
p8 = create_spatial_abun('northern pike', 'Abundance', 'Northern Pike - Abundance', save_fig = F)
p9 = create_spatial_abun('walleye', 'CPUE', 'Walleye - CPUE', save_fig = F)
p10 = create_spatial_abun('walleye', 'Abundance', 'Walleye - Abundance', save_fig = F)
p11 = create_spatial_abun('yellow perch', 'CPUE', 'Yellow Perch - CPUE', save_fig = F)
p12 = create_spatial_abun('yellow perch', 'Abundance', 'Yellow Perch - Abundance', save_fig = F)

rel_spat_1 = cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, align = 'hv', byrow = F)
rel_spat_2 = cowplot::plot_grid(p7, p8, p9, p10, p11, p12, nrow = 2, align = 'hv', byrow = F)

# ggsave('results/spatial_results/cpue_rel_spat1.png', rel_spat_1, width = 15, height = 8)
# ggsave('results/spatial_results/cpue_rel_spat2.png', rel_spat_2, width = 15, height = 8)


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

yr_1 = 2000
yr_2 = 2019
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
         Date=ymd(format(Date, paste0(yr_2, "-%m-%d")))) %>% 
  filter(complete.cases(SURVEYDATE))



#################3
ch_plot = cold %>% 
  rbind(hot) %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Year, fill = Year)) +
  geom_line(aes(linetype = Lake), size = 0.9) +
  facet_wrap(~ Fish + Gear) +
  geom_rug(aes(x = Date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr_2, '-06-01')), ymd(paste0(yr_2, '-10-01')))) +
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

# ggsave(paste0('results/spatial_results/catchability_lake_comp_species_', 'cold_hot', '.png'), ch_plot, width = 14, height = 10)


cold %>% 
  rbind(hot) %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = exp(Effectiveness), color = Year, fill = Year)) +
  geom_line(aes(linetype = Lake), size = 0.9) +
  facet_wrap(~ Fish + Gear) +
  geom_rug(aes(x = Date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr_2, '-06-01')), ymd(paste0(yr_2, '-10-01')))) +
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
  ylim(c(0, 180000))

tmp_dat = cold %>% 
  rbind(hot) %>% 
  mutate(Eff = exp(Effectiveness))


catch_df = cold %>% 
  mutate(SURVEYDATE=ymd(format(SURVEYDATE, "2019-%m-%d")),
         Date=ymd(format(Date, "2019-%m-%d"))) %>% 
  rbind(hot) %>% 
  mutate(Year = factor(Year),
         Lake = factor(Lake, levels = c('North', 'Center', 'South')))
f_plots = unique(catch_df$Fish)
f_plot_title = str_to_title(pars$fish_names)

for(i in seq_along(f_plots)){
  
  curr_df = catch_df %>% 
    filter(Fish == f_plots[i])
  
  ggplot(curr_df, aes(x = SURVEYDATE, y = Effectiveness, color = Lake)) +
    geom_line(aes(linetype = Year), size = 1.3) +
    facet_wrap(~ Gear, ncol = 2) +
    geom_rug(aes(x = Date), sides="b") +
    scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr_2, '-06-01')), ymd(paste0(yr_2, '-10-01')))) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_color_manual(values = c('blue', 'red', 'black')) +
    ggtitle(paste0('Gear Type Effectiveness by Region - ', f_plot_title[i])) +
    xlab('') +
    ylab('Relative Effectiveness') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = 'bold'),
          title = element_text(size = 16, face = 'bold'),
          legend.position = "bottom",
          legend.key.width = unit(1.5, 'cm'),
          strip.text = element_text(size = 14),
          legend.spacing = unit(5, 'cm'))
  
  fpath = paste0('results/spatial_results/catchability_lake_comp_', f_plots[i], '.png')
  ggsave(fpath, width = 14, height = 8)
  
}

lake_dow_south = 53002800 # south
lake_dow_north = 16001900 # north
lake_dow_center = 18005000 # center

lake_loc = tibble(DOW = c(lake_dow_south, lake_dow_north, lake_dow_center), lake = c('South', 'North', 'Center')) %>% 
  mutate(DOW = as.factor(DOW))

lake_loc = effort %>% left_join(static, by = 'DOW') %>% 
            select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
            mutate(DOW = factor(DOW)) %>% 
            distinct(DOW, .keep_all = T) %>% 
  right_join(lake_loc)


lats = c(43.50852, 48.72190)
lons = c(-96.75060, -90.06664)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot(data = lake_loc) +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5), size = 3) +
  geom_text(aes(label=lake, x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5), hjust=0.5, vjust=-1.1, size = 8) + 
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Lake Location") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave(paste0('results/spatial_results/lake_location.png'), width = 8, height = 8)

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


# average number of lakes surveyed per year
fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  summarize(number = n()) %>% 
  summarize(min = min(number), max = max(number), mean = mean(number), median = median(number), sd = sd(number)) %>% 
  xtable(caption = 'Number of repeat sampler per lake per year.')

# hist of number of lakes per year
fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  summarize(number = n()) %>% 
  arrange(desc(number)) %>% 
  ggplot(., aes(x = number)) +
  geom_histogram(bins = 23)

fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  summarize(number = n()) %>% 
  filter(number == 1)

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


# transformed covariate summary
fish_dat %>% 
  select(all_of(mean_covs)) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
  mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
  mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
  rename('depth' = 'MAX_DEPTH_FEET',
         'area' = 'LAKE_AREA_GIS_ACRES') %>% 
  summarise_all(list(min = ~min(.),
                     max = ~max(.), 
                     mean = ~mean(.),
                     median = ~median(.),
                     sd = ~sd(.))) %>% 
  pivot_longer(depth_min:secchi_sd, names_to = c('Variable', 'Summary'), names_sep = '_', values_to = 'Value') %>% 
  pivot_wider(names_from = 'Summary', values_from = 'Value') %>% 
  xtable(caption = 'Summary of covariates used. All covariates are transformed as they would appear within the model.')

# raw covariate summary
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
  pivot_wider(names_from = 'Summary', values_from = 'Value') %>% 
  xtable(caption = 'Summary of covariates used. All covariates are on their original, untransformed scales.')



# fish catch summaries
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
  pivot_wider(names_from = 'Summary', values_from = 'Value') %>% 
  xtable(caption = 'Summaries of fish caught by lake for the whole data set, comparing total caught to CPUE.')


# fish catch summaries
fish_dat %>% 
  select(DOW:CPUE) %>% 
  rename('fish' = 'COMMON_NAME',
         'total' = 'TOTAL_CATCH') %>% 
  rename_all(~str_to_lower(.)) %>% 
  group_by(fish) %>% 
  summarise_at(vars(total, cpue), list(min = ~quantile(., probs = 0.01),
                                       max = ~quantile(., probs = 0.99), 
                                       mean = ~mean(.),
                                       median = ~median(.),
                                       sd = ~sd(.))) %>% 
  pivot_longer(total_min:cpue_sd, names_to = c('Variable', 'Summary'), names_sep = '_', values_to = 'Value') %>% 
  pivot_wider(names_from = 'Summary', values_from = 'Value')



tmp = fish_dat %>% 
  group_by(COMMON_NAME) %>% 
  filter(TOTAL_CATCH == max(TOTAL_CATCH))

77018100

tmp2 = fish_dat %>% 
  filter(DOW == '77018100') %>% 
  filter(SURVEYDATE == '1996-04-30')

tmp3 = spt_loc %>% 
  filter(DOW == '77018100')

spt_loc = effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T)

spt_loc %>% 
  filter(DOW == '77018100')

fish_dat %>% 
  filter(DOW == '77018100') %>% 
  filter(COMMON_NAME == 'yellow perch')


range(fish_dat$SURVEYDATE)

length(unique(fish_dat$DOW))






by_fish_vals = fish_dat %>% 
  select(all_of(mean_covs), COMMON_NAME, TOTAL_CATCH) %>% 
  rename('depth' = 'MAX_DEPTH_FEET',
         'area' = 'LAKE_AREA_GIS_ACRES') %>% 
  filter(TOTAL_CATCH > 0) %>% 
  select(-TOTAL_CATCH) %>% 
  group_by(COMMON_NAME) %>% 
  summarise_all(list(min = ~min(.),
                     max = ~max(.), 
                     mean = ~mean(.),
                     median = ~median(.),
                     sd = ~sd(.))) %>% 
  pivot_longer(depth_min:secchi_sd, names_to = c('Variable', 'Summary'), names_sep = '_', values_to = 'Value') %>% 
  pivot_wider(names_from = 'Summary', values_from = 'Value')

by_fish_vals %>% 
  filter(Variable == 'secchi')

by_fish_vals %>% 
  filter(Variable == 'DD5')


b_hat = apply(beta[,,-c(1:5000)], c(1,2), mean)

round(b_hat, 3) %>%
  as_tibble(.name_repair = ~b_names) %>%
  mutate(Species = pars$fish_names) %>%
  relocate(Species)




# relative abundance vs. cpue unitless ------------------------------------

b_hat = apply(beta, c(1,2), mean)
colnames(b_hat) = b_names

# construct omega
K = length(pars$fish_names)
# omega_hat = pars$P %*% matrix(apply(omega, 1, mean), ncol = K, byrow = F)
omega_hat = pars$P %*% apply(omega, c(1,2), mean)
ind_array = tibble(id = pars$lake_id, as_tibble(omega_hat, .name_repair = ~LETTERS[1:K]))
lake_array = tibble(id = pars$lake_index)
omega_hat = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))


lam_hat = list()
for(k in 1:pars$K){
  lam_hat[[k]] = exp(pars$X[[k]] %*% b_hat[k,] + omega_hat[,k])
}


rel_abun = rbind(fish_dat %>% filter(COMMON_NAME == 'black crappie') %>% mutate(Abundance = c(lam_hat[[1]])),
                 fish_dat %>% filter(COMMON_NAME == 'bluegill') %>% mutate(Abundance = c(lam_hat[[2]])),
                 fish_dat %>% filter(COMMON_NAME == 'largemouth bass') %>% mutate(Abundance = c(lam_hat[[3]])),
                 fish_dat %>% filter(COMMON_NAME == 'northern pike') %>% mutate(Abundance = c(lam_hat[[4]])),
                 fish_dat %>% filter(COMMON_NAME == 'walleye') %>% mutate(Abundance = c(lam_hat[[5]])),
                 fish_dat %>% filter(COMMON_NAME == 'yellow perch') %>% mutate(Abundance = c(lam_hat[[6]])))

rel_abun = rel_abun %>% 
  select(DOW, COMMON_NAME, CPUE, Abundance) %>% 
  group_by(DOW, COMMON_NAME) %>% 
  summarise_at(vars(CPUE, Abundance), median) %>% 
  ungroup() %>% 
  mutate(Fish = str_replace_all(COMMON_NAME, " ", "_")) %>% 
  left_join(effort %>% 
              left_join(static, by = 'DOW') %>% 
              select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
              mutate(DOW = factor(DOW)) %>% 
              distinct(DOW, .keep_all = T), by = 'DOW') %>% 
  select(DOW, COMMON_NAME, CPUE, Abundance, Fish, LAKE_CENTER_LONG_DD5, LAKE_CENTER_LAT_DD5) %>% group_by(COMMON_NAME) %>% 
  mutate_at(vars(CPUE, Abundance), ~ (. - mean(.))/sd(.)) %>% 
  ungroup() %>% 
  pivot_longer(c(Abundance, CPUE), names_to = 'Est', values_to = 'Abun')

tmp = rel_abun %>% 
  filter(COMMON_NAME == 'northern pike') %>% 
  filter(LAKE_CENTER_LONG_DD5 < -95.5 & LAKE_CENTER_LONG_DD5 > -96.5) %>% 
  filter(LAKE_CENTER_LAT_DD5 < 45 & LAKE_CENTER_LAT_DD5 > 44)

lats = range(rel_abun$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(rel_abun$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))


create_spatial_abun <- function(fish, abun, title, save_fig = FALSE){
  
  p = ggplot() +
    geom_sf(data = usa) +
    coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
    geom_point(data = rel_abun %>% 
                 filter(Est == abun) %>% 
                 filter(COMMON_NAME == fish), aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = Abun), size = 1.5) +
    # scale_color_gradientn(colors = c('red', 'white', 'blue')) +
    scale_color_fermenter(breaks = c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3), palette = 'RdBu', limits = c(-4, 4)) + #
    xlab("Longitude") +
    ylab("Latitude") +
    theme(legend.title=element_blank(),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = 'bold'),
          title = element_text(size = 16, face = 'bold')) +
    ggtitle(title)
  
  if(save_fig == TRUE){
    path = paste0('/Users/joshuanorth/Desktop/erin_talk_plots/', str_replace_all(fish, ' ', '_'), '_', abun, '.png')
    ggsave(path, p, width = 8, height = 8)
  }else{
    return(p)
  }
  
}

create_spatial_abun('black crappie', 'CPUE', 'Black Crappie - CPUE', save_fig = F)
create_spatial_abun('black crappie', 'Abundance', 'Black Crappie - Abundance', save_fig = F)
create_spatial_abun('bluegill', 'CPUE', 'Bluegill - CPUE', save_fig = F)
create_spatial_abun('bluegill', 'Abundance', 'Bluegill - Abundance', save_fig = F)
create_spatial_abun('largemouth bass', 'CPUE', 'Largemouth Bass - CPUE', save_fig = F)
create_spatial_abun('largemouth bass', 'Abundance', 'Largemouth Bass - Abundance', save_fig = F)
create_spatial_abun('northern pike', 'CPUE', 'Northern Pike - CPUE', save_fig = F)
create_spatial_abun('northern pike', 'Abundance', 'Northern Pike - Abundance', save_fig = F)
create_spatial_abun('walleye', 'CPUE', 'Walleye - CPUE', save_fig = F)
create_spatial_abun('walleye', 'Abundance', 'Walleye - Abundance', save_fig = F)
create_spatial_abun('yellow perch', 'CPUE', 'Yellow Perch - CPUE', save_fig = F)
create_spatial_abun('yellow perch', 'Abundance', 'Yellow Perch - Abundance', save_fig = F)

create_spatial_abun('black crappie', 'CPUE', 'Black Crappie - CPUE', save_fig = T)
create_spatial_abun('black crappie', 'Abundance', 'Black Crappie - Abundance', save_fig = T)
create_spatial_abun('bluegill', 'CPUE', 'Bluegill - CPUE', save_fig = T)
create_spatial_abun('bluegill', 'Abundance', 'Bluegill - Abundance', save_fig = T)
create_spatial_abun('largemouth bass', 'CPUE', 'Largemouth Bass - CPUE', save_fig = T)
create_spatial_abun('largemouth bass', 'Abundance', 'Largemouth Bass - Abundance', save_fig = T)
create_spatial_abun('northern pike', 'CPUE', 'Northern Pike - CPUE', save_fig = T)
create_spatial_abun('northern pike', 'Abundance', 'Northern Pike - Abundance', save_fig = T)
create_spatial_abun('walleye', 'CPUE', 'Walleye - CPUE', save_fig = T)
create_spatial_abun('walleye', 'Abundance', 'Walleye - Abundance', save_fig = T)
create_spatial_abun('yellow perch', 'CPUE', 'Yellow Perch - CPUE', save_fig = T)
create_spatial_abun('yellow perch', 'Abundance', 'Yellow Perch - Abundance', save_fig = T)

# difference in catchability cold/hot -------------------------------------



d_plot_lake = function(phis, fish_names, yr, fish_dat){
  
  fish_names = str_replace_all(fish_names, " ", "_")
  
  tC = temp %>%
    filter(year(SURVEYDATE) == yr) %>%
    group_by(SURVEYDATE) %>% 
    summarise(temp_0 = mean(temp_0))
  
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
  GN = apply(post_stores_GN, c(1,2), mean)
  
  # TN = apply(post_stores_TN, c(1,2), median)
  # GN = apply(post_stores_GN, c(1,2), median)
  
  colnames(TN) = paste0("TN_",fish_names)
  colnames(GN) = paste0("GN_",fish_names)
  
  
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
  
  return(mean)
}

yr_1 = 2000
yr_2 = 2019
fnames = c('crappie', 'bluegill', 'bass', 'pike', 'walleye', 'perch')
plt_dat_1 = d_plot_lake(phi, fnames, yr_1, fish_dat)
plt_dat_2 = d_plot_lake(phi, fnames, yr_2, fish_dat)


catch_plt = plt_dat_1 %>% 
  select(-c(GN_survey:TN_survey)) %>% 
  pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "catch") %>% 
  left_join(plt_dat_1 %>%
              select(SURVEYDATE:TN_survey) %>% 
              pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "Date"), by = c('Gear', 'SURVEYDATE')) %>% 
  mutate(SURVEYDATE=ymd(format(SURVEYDATE, paste0(yr_2, "-%m-%d"))),
         Date=ymd(format(Date, paste0(yr_2, "-%m-%d")))) %>% 
  select(-sample_day) %>% 
  mutate(Year = yr_1) %>% 
  rbind(plt_dat_2 %>% 
          select(-c(GN_survey:TN_survey)) %>% 
          pivot_longer(cols = -c(SURVEYDATE), names_to = c("Gear", "Fish"), names_sep = "_", values_to = "catch") %>% 
          left_join(plt_dat_2 %>%
                      select(SURVEYDATE:TN_survey) %>% 
                      pivot_longer(GN_survey:TN_survey, names_to = c('Gear', 'sample_day'), names_sep = "_", values_to = "Date"), by = c('Gear', 'SURVEYDATE')) %>% 
          select(-sample_day) %>% 
          mutate(Year = yr_2)) %>% 
  mutate(Year = as.factor(Year)) %>% 
  filter(SURVEYDATE > ymd(paste0(yr_2, '-06-01')) & SURVEYDATE < ymd(paste0(yr_2, '-09-15'))) %>% 
  group_by(Fish) %>% 
  mutate(catch = (catch - min(catch))/(max(catch) - min(catch)))
  


ggplot(catch_plt, aes(x = SURVEYDATE, y = catch, color = Year, fill = Year)) +
  geom_line(aes(linetype = Gear), size = 0.9) +
  facet_wrap(~ Fish) +
  geom_rug(aes(x = Date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr_2, '-06-01')), ymd(paste0(yr_2, '-09-15')))) +
  scale_linetype_manual(values=c("solid", 'longdash')) +
  ggtitle(paste0('Gear Type Effectiveness by Region')) +
  xlab('') +
  ylab('Relative Effectiveness') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold')) +
  ylim(c(-0.25, 1.25))

# ggsave('/Users/joshuanorth/Desktop/erin_talk_plots/all_catch.png', width = 14, height = 8)


f_plots = unique(catch_plt$Fish)
f_plot_title = str_to_title(pars$fish_names)

for(i in seq_along(f_plots)){
  
  curr_df = catch_plt %>% 
    filter(Fish == f_plots[i])
  
  ggplot(curr_df, aes(x = SURVEYDATE, y = catch, color = Year, fill = Year)) +
    geom_line(aes(linetype = Gear), size = 1.1) +
    geom_rug(aes(x = Date), sides="b") +
    scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr_2, '-06-01')), ymd(paste0(yr_2, '-09-15')))) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_color_manual(values = c('blue', 'black')) +
    ggtitle(paste0('Gear Type Effectiveness by Region - ', f_plot_title[i])) +
    xlab('') +
    ylab('Relative Effectiveness') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = 'bold'),
          title = element_text(size = 16, face = 'bold'),
          legend.position = "bottom",
          legend.key.width = unit(1.5, 'cm'),
          strip.text = element_text(size = 14),
          legend.spacing = unit(5, 'cm'))
  
  fpath = paste0('/Users/joshuanorth/Desktop/erin_talk_plots/catchability_lake_comp_', f_plots[i], '.png')
  ggsave(fpath, width = 14, height = 8)
  
}

#################3
ch_plot = cold %>% 
  rbind(hot) %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = Effectiveness, color = Year, fill = Year)) +
  geom_line(aes(linetype = Lake), size = 0.9) +
  facet_wrap(~ Fish + Gear) +
  geom_rug(aes(x = Date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr_2, '-06-01')), ymd(paste0(yr_2, '-10-01')))) +
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

# ggsave(paste0('results/spatial_results/catchability_lake_comp_species_', 'cold_hot', '.png'), ch_plot, width = 14, height = 10)


cold %>% 
  rbind(hot) %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(., aes(x = SURVEYDATE, y = exp(Effectiveness), color = Year, fill = Year)) +
  geom_line(aes(linetype = Lake), size = 0.9) +
  facet_wrap(~ Fish + Gear) +
  geom_rug(aes(x = Date), sides="b") +
  scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr_2, '-06-01')), ymd(paste0(yr_2, '-10-01')))) +
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
  ylim(c(0, 180000))

tmp_dat = cold %>% 
  rbind(hot) %>% 
  mutate(Eff = exp(Effectiveness))


catch_df = cold %>% 
  mutate(SURVEYDATE=ymd(format(SURVEYDATE, "2019-%m-%d")),
         Date=ymd(format(Date, "2019-%m-%d"))) %>% 
  rbind(hot) %>% 
  mutate(Year = factor(Year),
         Lake = factor(Lake, levels = c('North', 'Center', 'South')))
f_plots = unique(catch_df$Fish)
f_plot_title = str_to_title(pars$fish_names)

for(i in seq_along(f_plots)){
  
  curr_df = catch_df %>% 
    filter(Fish == f_plots[i])
  
  ggplot(curr_df, aes(x = SURVEYDATE, y = Effectiveness, color = Lake)) +
    geom_line(aes(linetype = Year), size = 1.3) +
    facet_wrap(~ Gear, ncol = 2) +
    geom_rug(aes(x = Date), sides="b") +
    scale_x_date(date_breaks = 'month', date_labels = "%b %d", limits = c(ymd(paste0(yr_2, '-06-01')), ymd(paste0(yr_2, '-10-01')))) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_color_manual(values = c('blue', 'red', 'black')) +
    ggtitle(paste0('Gear Type Effectiveness by Region - ', f_plot_title[i])) +
    xlab('') +
    ylab('Relative Effectiveness') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = 'bold'),
          title = element_text(size = 16, face = 'bold'),
          legend.position = "bottom",
          legend.key.width = unit(1.5, 'cm'),
          strip.text = element_text(size = 14),
          legend.spacing = unit(5, 'cm'))
  
  fpath = paste0('results/spatial_results/catchability_lake_comp_', f_plots[i], '.png')
  ggsave(fpath, width = 14, height = 8)
  
}

lake_dow_south = 53002800 # south
lake_dow_north = 16001900 # north
lake_dow_center = 18005000 # center

lake_loc = tibble(DOW = c(lake_dow_south, lake_dow_north, lake_dow_center), lake = c('South', 'North', 'Center')) %>% 
  mutate(DOW = as.factor(DOW))

lake_loc = effort %>% left_join(static, by = 'DOW') %>% 
  select(DOW, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5) %>% 
  mutate(DOW = factor(DOW)) %>% 
  distinct(DOW, .keep_all = T) %>% 
  right_join(lake_loc)


lats = c(43.50852, 48.72190)
lons = c(-96.75060, -90.06664)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

ggplot(data = lake_loc) +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5), size = 3) +
  geom_text(aes(label=lake, x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5), hjust=0.5, vjust=-1.1, size = 8) + 
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Lake Location") +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, 'cm'))

# ggsave(paste0('results/spatial_results/lake_location.png'), width = 8, height = 8)

# infor for lake for erin -------------------------------------------------

lake_num = '42009600'
fish_dat %>% 
  filter(DOW == lake_num) %>% 
  filter(COMMON_NAME =='northern pike')

tmp = fish_dat %>% 
  filter(DOW == lake_num)

write_csv(tmp, '/Users/joshuanorth/Desktop/erin_talk_plots/question_lake.csv')

tmp = fish_dat %>% 
  filter(DOW == lake_num) %>% 
  filter(COMMON_NAME =='northern pike')

write_csv(tmp, '/Users/joshuanorth/Desktop/erin_talk_plots/question_lake_pike.csv')