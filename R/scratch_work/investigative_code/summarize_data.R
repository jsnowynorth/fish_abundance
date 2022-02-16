library(tidyverse)
library(lubridate)
library(viridis)
library(sparklyr)
library(stringr)
library(readxl)


# load in fish data
fish_dat = read_csv('data/fish_dat.csv')

fish_dat = fish_dat %>% 
  mutate(DOW = as.factor(DOW),
         COMMON_NAME = as.factor(COMMON_NAME)) %>% 
  arrange(DOW, year, COMMON_NAME)

# run this to only include data after 2001
fish_dat = fish_dat %>%
  filter(year >= 2000)


# average number of lakes surveyed per year
fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  summarize(number = n()) %>% 
  summarize(min = min(number), max = max(number), mean = mean(number), median = median(number), sd = sd(number))

# histogram in numbers
fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  summarize(number = n()) %>% 
  mutate(int = 1) %>% 
  group_by(number) %>% 
  summarize(n_lakes = n())

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

# average survey date day-month
fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  mutate(SURVEYDATE=ymd(format(SURVEYDATE, "2016-%m-%d"))) %>% 
  ungroup() %>% 
  summarize(min = min(SURVEYDATE), max = max(SURVEYDATE), mean = mean(SURVEYDATE), median = median(SURVEYDATE), sd = sd(SURVEYDATE)) %>% 
  mutate_all(~format(., format="%m-%d"))

# average survey date doth  
fish_dat %>% 
  select(DOW, SURVEYDATE) %>% 
  group_by(DOW, SURVEYDATE) %>% 
  unique() %>% 
  ungroup(SURVEYDATE) %>% 
  mutate(day=yday(SURVEYDATE)) %>% 
  ungroup() %>% 
  summarize(min = min(day), max = max(day), mean = mean(day), median = median(day), sd = sd(day))

# summary of covariates
fish_dat %>% 
  select(all_of(colnames(fish_dat)[c(7, 9, 23, 13:15,17, 25)])) %>% 
  rename('depth' = 'MAX_DEPTH_FEET',
         'area' = 'LAKE_AREA_GIS_ACRES') %>% 
  summarise_all(list(min = ~min(.),
                     max = ~max(.), 
                     mean = ~mean(.),
                     median = ~median(.),
                     sd = ~sd(.))) %>% 
  pivot_longer(depth_min:secchi_sd, names_to = c('Variable', 'Summary'), names_sep = '_', values_to = 'Value') %>% 
  pivot_wider(names_from = 'Summary', values_from = 'Value')

# summary of total catch and cpue
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

# more info about the max values
fish_dat %>% 
  select(DOW:CPUE) %>% 
  group_by(COMMON_NAME) %>% 
  filter(TOTAL_CATCH == max(TOTAL_CATCH)) %>% 
  ungroup()

# summary of total catch and cpue with 0.01 and 0.99 quantile
fish_dat %>% 
  select(DOW:CPUE) %>% 
  rename('fish' = 'COMMON_NAME',
         'total' = 'TOTAL_CATCH') %>% 
  rename_all(~str_to_lower(.)) %>% 
  group_by(fish) %>% 
  summarise_at(vars(total, cpue), list(lower = ~quantile(., probs = 0.01),
                                       upper = ~quantile(., probs = 0.99), 
                                       mean = ~mean(.),
                                       median = ~median(.),
                                       sd = ~sd(.))) %>% 
  pivot_longer(total_lower:cpue_sd, names_to = c('Variable', 'Summary'), names_sep = '_', values_to = 'Value') %>% 
  pivot_wider(names_from = 'Summary', values_from = 'Value')


# more info about the values above the 999 quantile
fish_dat %>% 
  select(DOW:CPUE) %>% 
  group_by(COMMON_NAME) %>% 
  filter(TOTAL_CATCH > quantile(TOTAL_CATCH, probs = 0.999)) %>% 
  ungroup()


