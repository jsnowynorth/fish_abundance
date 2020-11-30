##----------------------------------------------------
## Name: Joshua North
##
## Date: 10/19/2020
##
## Project: Fish Abundance
##
## Objective: Simple Poisson model for northern pike, no time
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
library(progress)



# load in data ------------------------------------------------------------

static = read_csv('data/Static_MN_Data_JE.csv')
effort = read_csv('data/MN_fish_catch_effort_all_lakes_years_surveys_JE.csv')



# join data ---------------------------------------------------------------


static = static %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))
effort = effort %>% mutate(DOW = (str_pad(DOW, 8, side = "left", pad = "0")))

# select only northern pike
# fish_dat = effort %>% 
#   left_join(static, by = 'DOW') %>%
#   filter(COMMON_NAME == 'northern pike') %>%
#   select(-c(SURVEY_ID, LKNAME, CTY:SURVEYDATE, COMMON_NAME, EFFORT, TOTAL_CATCH, CTY:SURVEYTYPE, FISHERIES_WATERBODY_ID.x:ALT_LAKE_NAME, LAKE_CENTER_LAT_DD5:LAKE_CENTER_UTM_NORTHING, FETCH_ORIENTATION_DEGREES)) %>%
#   mutate(DOW = as.factor(DOW),
#          GEAR = as.factor(GEAR)) %>%
#   group_by(.dots=c('DOW', 'GEAR')) %>%
#   summarise_all( ~ mean(., na.rm = T)) %>% 
#   ungroup() %>% 
#   filter(complete.cases(.)) %>% 
#   mutate(Gear_ind = as.integer(1)) %>% 
#   pivot_wider(names_from = GEAR, values_from = Gear_ind, values_fill = 0)
# 
# rm(effort, static)


# select only northern pike
# fish_dat = effort %>% 
#   left_join(static, by = 'DOW') %>%
#   filter(COMMON_NAME == 'northern pike') %>%
#   select(-c(SURVEY_ID, LKNAME, CTY:SURVEYDATE, COMMON_NAME, CPUE, EFFORT, CTY:SURVEYTYPE, AVG_MEASURED_WT_LBS:ALT_LAKE_NAME, 
#             LAKE_CENTER_LAT_DD5:LAKE_CENTER_UTM_NORTHING, FETCH_ORIENTATION_DEGREES, PDIST_cum, littoral.zone)) %>%
#   mutate(DOW = as.factor(DOW),
#          GEAR = as.factor(GEAR)) %>%
#   group_by(.dots=c('DOW', 'GEAR')) %>%
#   summarise_all( ~ floor(mean(., na.rm = T))) %>% 
#   ungroup() %>% 
#   filter(complete.cases(.)) %>% 
#   mutate(Gear_ind = as.integer(1)) %>% 
#   pivot_wider(names_from = GEAR, values_from = Gear_ind, values_fill = 0) %>% 
#   mutate(DOW_ind = as.integer(1)) %>% 
#   pivot_wider(names_from = DOW, values_from = DOW_ind, values_fill = 0)
# 
# rm(effort, static)

dat_tmp = effort %>% 
  left_join(static, by = 'DOW')


# select only northern pike
fish_dat = effort %>% 
  left_join(static, by = 'DOW') %>%
  filter(COMMON_NAME == 'northern pike') %>%
  select(-c(SURVEY_ID, LKNAME, CTY:SURVEYTYPE, COMMON_NAME, CPUE, CTY:SURVEYTYPE, AVG_MEASURED_WT_LBS:ALT_LAKE_NAME, 
            LAKE_CENTER_LAT_DD5:LAKE_CENTER_UTM_NORTHING, FETCH_ORIENTATION_DEGREES, PDIST_cum, littoral.zone)) %>%
  mutate(DOW = as.factor(DOW),
         GEAR = as.factor(GEAR),
         DOY = yday(mdy(SURVEYDATE)),
         DOY = (DOY - mean(DOY))/sd(DOY),
         DOY_squared = DOY^2) %>%
  filter(complete.cases(.)) %>% 
  mutate(row = row_number(),
         Gear_ind = as.integer(1)) %>% 
  pivot_wider(names_from = GEAR, values_from = Gear_ind, values_fill = 0) %>% 
  select(-c(SURVEYDATE, row))
# %>% 
#   mutate(DOW_ind = as.integer(1)) %>% 
#   pivot_wider(names_from = DOW, values_from = DOW_ind, values_fill = 0) %>% 
#   select(-row)

# rm(effort, static)

hist(log(fish_dat$AVG_WATER_CLARITY))
hist(log(fish_dat$mean.gdd))
hist(log(fish_dat$area.hectares))
hist(log(fish_dat$LAKE_AREA_DOW_ACRES))
hist(log(fish_dat$LAKE_AREA_GIS_ACRES))
hist(log(fish_dat$LAKE_AREA_PLANIMETERED_ACRES))
hist(log(fish_dat$MAX_DEPTH_FEET))
hist(log(fish_dat$SHORE_LENGTH_MILES))
hist(log(fish_dat$MAX_FETCH_MILES))
hist(log(fish_dat$MEAN_DEPTH_FEET))
hist(log(fish_dat$LITTORAL_AREA_ACRES))
hist(log(fish_dat$secchi.m))
hist(fish_dat$DOY)


# simple model ------------------------------------------------------------

create_pars <- function(fish_dat){
  
  pars = list()
  
  # data
  pars$Y = fish_dat$TOTAL_CATCH
  # X = fish_dat %>% 
  #   select(mean.gdd, MAX_DEPTH_FEET, LAKE_AREA_GIS_ACRES, secchi.m:DOY_squared) %>% #select(AVG_WATER_CLARITY, mean.gdd:DOY_squared) %>% 
  #   mutate_at(vars(-c(DOY:DOY_squared)), ~ ifelse(. == 0, . + 0.001, .)) %>% 
  #   mutate_at(vars(-c(DOY:DOY_squared)), ~ log(.)) %>% 
  #   mutate(Int = 1)
  
  X = fish_dat %>% 
    select(LAKE_AREA_GIS_ACRES, secchi.m:DOY_squared) %>% #select(AVG_WATER_CLARITY, mean.gdd:DOY_squared) %>% 
    mutate_at(vars(-c(DOY:DOY_squared)), ~ ifelse(. == 0, . + 0.001, .)) %>% 
    mutate_at(vars(-c(DOY:DOY_squared)), ~ log(.)) %>% 
    select(-LAKE_AREA_GIS_ACRES)
  
  Z = fish_dat %>% 
    select(GN:GSH) %>% 
    mutate(Int = 1) %>% 
    select(-GN)
  
  pars$X = as.matrix(cbind(X, Z))
  pars$effort = fish_dat$EFFORT
  
  # parameters
  pars$n = length(pars$Y)
  pars$p = ncol(pars$X)
  
  pars$beta = rep(0, pars$p)
  pars$beta_accept = rep(0, pars$p)
  
  # hyperpriors
  pars$sig_prop = rep(0.5, pars$p)
  
  return(pars)
  
}

update_beta <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  effort = pars$effort
  
  # parameters
  n = pars$n
  p = pars$p
  
  beta_accept =  rep(0, pars$p)
  
  beta_curr = pars$beta
  
  sig_prop = pars$sig_prop

  # proposal
  # sig = sig_prop * diag(length(beta_curr))
  # beta_prop = t(rmvnorm(1, mean = beta_curr, sigma = sig))
  # 
  # like_curr = sum(dpois(Y, lambda = exp(X %*% beta_curr), log = T))
  # like_prop = sum(dpois(Y, lambda = exp(X %*% beta_prop), log = T))
  # 
  # 
  for(i in 1:p){
    
    b_prop = rnorm(1, beta_curr[i], sig_prop[i])
    beta_prop = beta_curr
    beta_prop[i] = b_prop
    
    like_curr = sum(dpois(Y, lambda = effort*exp(X %*% beta_curr), log = T))
    like_prop = sum(dpois(Y, lambda = effort*exp(X %*% beta_prop), log = T))
    
    if((like_prop - like_curr) > log(runif(1))){
      beta_curr[i] = b_prop
      beta_accept[i] = 1
    }
    
    
  }
  
  pars$beta = beta_curr
  pars$beta_accept = beta_accept
  
  return(pars)
  
  
}

update_proposal_var <- function(pars, beta_accept_post, i){
  
  sig_prop = pars$sig_prop
  
  bp = beta_accept_post[,(i-99):i]
  accept_rate = apply(bp, 1, mean)
  
  sig_prop = ifelse(accept_rate < 0.25, sig_prop*0.7, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.7, sig_prop)
  
  pars$sig_prop = sig_prop
  
  return(pars)
  
  
}

sampler <- function(fish_dat, nits){
  
  pars = create_pars(fish_dat)
  
  p = pars$p
  
  beta_post = array(NA, dim = c(p, nits))
  beta_accept_post = array(NA, dim = c(p, nits))
  
  pb <- progress_bar$new(
    format = "  Running [:bar] :percent eta: :eta",
    total = nits, clear = FALSE, width = 60)
  
  for(i in seq(1, nits)){
    
    pars <- update_beta(pars)

    beta_post[,i] = pars$beta
    beta_accept_post[,i] = pars$beta_accept

    
    if(i %in% seq(0, nits-1, by = 100)){
      pars <- update_proposal_var(pars, beta_accept_post, i)
    }
    
    
    pb$tick()
    
  }
  
  return(list(beta = beta_post,
              beta_accept = beta_accept_post,
              sig_prop = pars$sig_prop))
  
}

nits = 1000
burnin = 1:500

run = sampler(fish_dat, nits)

pars = create_pars(fish_dat)
nms = colnames(pars$X)


plot(run$beta[1,-c(burnin)], type = 'l', main = nms[1])
plot(run$beta[2,-c(burnin)], type = 'l', main = nms[2])
plot(run$beta[3,-c(burnin)], type = 'l', main = nms[3])
plot(run$beta[4,-c(burnin)], type = 'l', main = nms[4])
plot(run$beta[5,-c(burnin)], type = 'l', main = nms[5])
plot(run$beta[6,-c(burnin)], type = 'l', main = nms[6])
plot(run$beta[7,-c(burnin)], type = 'l', main = nms[7])
plot(run$beta[8,-c(burnin)], type = 'l', main = nms[8])
plot(run$beta[9,-c(burnin)], type = 'l', main = nms[9])
plot(run$beta[10,-c(burnin)], type = 'l', main = nms[10])
plot(run$beta[11,-c(burnin)], type = 'l', main = nms[11])
plot(run$beta[12,-c(burnin)], type = 'l', main = nms[12])
plot(run$beta[13,-c(burnin)], type = 'l', main = nms[13])
plot(run$beta[14,-c(burnin)], type = 'l', main = nms[14])
plot(run$beta[15,-c(burnin)], type = 'l', main = nms[15])
plot(run$beta[16,-c(burnin)], type = 'l', main = nms[16])
plot(run$beta[17,-c(burnin)], type = 'l', main = nms[17])
plot(run$beta[18,-c(burnin)], type = 'l', main = nms[18])
plot(run$beta[19,-c(burnin)], type = 'l', main = nms[19])
plot(run$beta[20,-c(burnin)], type = 'l', main = nms[20])

apply(run$beta_accept, 1, mean)
round(run$sig_prop, 3)


cbind(nms, apply(run$beta, 1, mean))



fish_dat_sd = effort %>% 
  left_join(static, by = 'DOW') %>%
  filter(COMMON_NAME == 'northern pike') %>%
  select(SURVEYDATE) %>%
  mutate(DOY = yday(mdy(SURVEYDATE))) %>% 
  summarise_at(vars(DOY), list(mean = ~mean(.), 
                                      sd = ~sd(.)))



d_plot = function(time){
  return(time*(-0.213704826181231) + (time^2)*(0.0665421242328541))
}

plot(d_plot((1:365 - fish_dat_sd$mean)/fish_dat_sd$sd), type = 'l')

# relative abundance ------------------------------------------------------

pars = create_pars(fish_dat)

b_hat = apply(run$beta[,-c(burnin)], 1, mean)
lam_hat = exp(pars$X[,-c(6:11)] %*% b_hat[-c(6:11)])
lam_hat_full = pars$effort*exp(pars$X %*% b_hat)
sum(pars$Y - lam_hat_full)/length(pars$Y)

plot(pars$Y - lam_hat)

res_dat = fish_dat %>% 
  mutate(rel_inten = lam_hat,
         rel_abun = lam_hat_full) %>% 
  left_join(static, by = 'DOW') %>% 
  select(c(DOW, rel_abun, rel_inten, TOTAL_CATCH, LAKE_CENTER_LAT_DD5, LAKE_CENTER_LONG_DD5)) %>% 
  mutate(e_hat = scale(TOTAL_CATCH - rel_abun))

lats = range(res_dat$LAKE_CENTER_LAT_DD5, na.rm = T)
lons = range(res_dat$LAKE_CENTER_LONG_DD5, na.rm = T)

usa = st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

p1 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = res_dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = rel_inten)) +
  scale_color_gradient(low = 'yellow', high = 'red') +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Intensity")

# ggsave('results/model_plots/intensity.png', p1)

p2 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = res_dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = e_hat)) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red', limits = c(-3,3)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Residuals")

# ggsave('results/model_plots/residuals.png', p2)

p3 <- ggplot() +
  geom_sf(data = usa) +
  coord_sf(xlim = c(lons[1] - 1, lons[2] + 1), ylim = c(lats[1] - 1, lats[2] + 1), expand = FALSE) +
  geom_point(data = res_dat, aes(x = LAKE_CENTER_LONG_DD5, y = LAKE_CENTER_LAT_DD5, color = log(rel_abun))) +
  scale_color_gradient(low = 'yellow', high = 'red', name = "Log Relative Abundance") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.title=element_blank()) +
  ggtitle("Log Relative Abundance")

# ggsave('results/model_plots/log_relative_abundance.png', p3)

p1
p2
p3


# effectiveness of gear ---------------------------------------------------

pars = create_pars(fish_dat)

b_hat = apply(run$beta[,-c(burnin)], 1, mean)
lam_hat = pars$effort * exp(pars$X %*% b_hat)
sum(pars$Y - lam_hat)/length(pars$Y)


gear_effect = as.matrix(fish_dat %>% select(GN:GSH))

# gear_ind = 
# fish_dat %>%
#   select(GN:GSH) %>% 
#   pivot_longer(GN:GSH, names_to = 'Gear', values_to = 'Ind') %>% 
#   filter(Ind == 1) %>% 
#   mutate(effect = c(pars$effort * exp(gear_effect %*% b_hat[15:20]))) %>% 
#   select(-Ind) %>% 
#   group_by(Gear) %>% 
#   summarize_all(mean)

fish_dat %>%
  select(GN:GSH) %>% 
  pivot_longer(GN:GSH, names_to = 'Gear', values_to = 'Ind') %>% 
  filter(Ind == 1) %>% 
  mutate(effect = c(exp(gear_effect %*% b_hat[6:11]))) %>% 
  select(-Ind) %>% 
  group_by(Gear) %>% 
  summarize_all(mean)



fish_dat %>%
  select(TOTAL_CATCH, EFFORT, GN:GSH) %>% 
  mutate(CPUE = TOTAL_CATCH/EFFORT) %>% 
  select(-c(TOTAL_CATCH, EFFORT)) %>% 
  pivot_longer(GN:GSH, names_to = 'Gear', values_to = 'Ind') %>% 
  filter(Ind == 1) %>% 
  select(-Ind) %>% 
  group_by(Gear) %>% 
  summarise_all(mean)










