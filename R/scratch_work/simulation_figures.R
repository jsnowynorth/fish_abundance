


# load libraries ----------------------------------------------------------

library(tidyverse)
library(rstan)
library(fields)
library(Matrix)
library(lubridate)
library(stringr)
library(LaplacesDemon)
library(xtable)



# load data ---------------------------------------------------------------

full = read_rds("data/stan_jsdm_simulation.rds")
chainsFull = extract(full, permuted = T)


omit = read_rds("data/stan_jsdm_without_effort_simulation.rds")
chainsOmit = extract(omit, permuted = T)

fnames = dat$fish_names
b_names = colnames(dat$X)

b_names = b_names %>% 
  str_replace_all("_", " ") %>% 
  str_to_title() %>% 
  str_replace_all('Dd5', 'GDD') %>% 
  str_replace_all('Gis', 'GIS')
b_names[1] = "Max Depth"
b_names[2] = "Lake Area"

b_names = c("Int", b_names)


# significance of no effort -----------------------------------------------


beta_0_no = chainsOmit$beta_0
beta_no = chainsOmit$beta

b_no_mean = t(apply(beta_no, c(2,3), mean))
b_no_lower = t(apply(beta_no, c(2,3), quantile, probs = 0.025))
b_no_upper = t(apply(beta_no, c(2,3), quantile, probs = 0.975))

b0_no_mean = apply(beta_0_no, c(2), mean)
b0_no_lower = apply(beta_0_no, c(2), quantile, probs = 0.025)
b0_no_upper = apply(beta_0_no, c(2), quantile, probs = 0.975)

b_no_mean = cbind(b0_no_mean, b_no_mean)
b_no_lower = cbind(b0_no_lower, b_no_lower)
b_no_upper = cbind(b0_no_upper, b_no_upper)

ind_array = array(NA, dim = dim(b_no_mean))
for(i in 1:nrow(b_no_mean)){
  for(j in 1:ncol(b_no_mean)){
    
    range1 = c(b_lower[i,j], b_upper[i,j])
    range2 = c(b_no_lower[i,j], b_no_upper[i,j])
    
    check1 = range1[1] < range2[1] & range1[2] > range2[1]
    check2 = range2[1] < range1[1] & range2[2] > range1[1]
    
    # check1 = range2[1] > range1[1] & range2[1] < range1[2]
    # check2 = range2[2] > range1[1] & range2[2] < range1[2]
    # check3 = sign(range1[1]) == sign(range2[1])
    # check4 = sign(range1[2]) == sign(range2[2])
    
    # ind_array[i,j] = ifelse((check1 | check2) & check3 & check4, 0, 1)
    ind_array[i,j] = ifelse((check1 | check2), 1, 0)
  }
}

b_no_sig = ifelse(ind_array == 1, 0, 1)

colnames(b_no_mean) = b_names
rownames(b_no_mean) = fnames

b_no_mean_sig = round(b_no_mean, 3)

b_no_col = ifelse(b_no_sig, paste0("\\B", round(b_no_mean, 3)),  round(b_no_mean, 3))
colnames(b_no_col) = b_names
rownames(b_no_col) = fnames
print(xtable(t(b_no_col)), sanitize.text.function = identity)


sign_ind_red = b_mean*b_no_sig > b_no_mean*b_no_sig
sign_ind_blue = b_mean*b_no_sig < b_no_mean*b_no_sig

b_no_color = ifelse(sign_ind_red, paste0("\\textcolor{red}{", b_no_mean_sig, "}"), b_no_mean_sig)
b_no_color = ifelse(sign_ind_blue, paste0("\\textcolor{blue}{", b_no_color, "}"), b_no_color)

b_no_col = ifelse(b_no_sig, paste0("\\B", b_no_color),  b_no_color)
colnames(b_no_col) = b_names
rownames(b_no_col) = fnames
print(xtable(t(b_no_col)), sanitize.text.function = identity)


rm(chains_no_catch)



# plot --------------------------------------------------------------------



b_mean = t(apply(chainsFull$beta, c(2,3), mean))
b_lower = t(apply(chainsFull$beta, c(2,3), quantile, probs = 0.025))
b_upper = t(apply(chainsFull$beta, c(2,3), quantile, probs = 0.975))

b0_mean = apply(chainsFull$beta_0, c(2), mean)
b0_lower = apply(chainsFull$beta_0, c(2), quantile, probs = 0.025)
b0_upper = apply(chainsFull$beta_0, c(2), quantile, probs = 0.975)

b_mean = cbind(b0_mean, b_mean)
b_lower = cbind(b0_lower, b_lower)
b_upper = cbind(b0_upper, b_upper)

colnames(b_mean) = b_names
colnames(b_lower) = b_names
colnames(b_upper) = b_names
# rownames(b_mean) = fnames

b_mean = as_tibble(b_mean) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Mean")

b_lower = as_tibble(b_lower) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Lower")

b_upper = as_tibble(b_upper) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Upper")

b_plt = b_mean %>% 
  left_join(b_lower) %>% 
  left_join(b_upper) %>% 
  mutate(sig = as.factor(if_else(sign(Lower) == sign(Upper), 1, 0))) %>% 
  mutate(Variable = factor(Variable, 
                           levels = c('Int', 'Max Depth', 'Lake Area', 'Ag', 'Urban', 'Wetlands', 'GDD', 'Secchi'))) %>% 
  mutate(species = factor(species,
                          levels = c('black crappie', 'bluegill', 'largemouth bass', 'northern pike', 'walleye', 'yellow perch')))


ggplot(b_plt, aes(y = reorder(species, desc(species)))) +
  geom_point(aes(x = Mean, shape = sig, color = sig), size = 1.5) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper, color = sig), width = 0.3) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  facet_wrap(~Variable, scales = 'free_x', ncol = 4) +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual(values = c('darkgray', 'blue')) +
  xlab('Parameter Estimate') +
  ylab('') +
  guides(shape = 'none', color = 'none') 

# ggsave('results/spatial_results/parameter_estimate_plot.png', width = 12, height = 6)
# ggsave('results/eps_figures/parameter_estimate_plot.eps', width = 12, height = 6, device="eps")


b_mean = t(apply(chainsOmit$beta, c(2,3), mean))
b_lower = t(apply(chainsOmit$beta, c(2,3), quantile, probs = 0.025))
b_upper = t(apply(chainsOmit$beta, c(2,3), quantile, probs = 0.975))

b0_mean = apply(chainsOmit$beta_0, c(2), mean)
b0_lower = apply(chainsOmit$beta_0, c(2), quantile, probs = 0.025)
b0_upper = apply(chainsOmit$beta_0, c(2), quantile, probs = 0.975)

b_mean = cbind(b0_mean, b_mean)
b_lower = cbind(b0_lower, b_lower)
b_upper = cbind(b0_upper, b_upper)

colnames(b_mean) = b_names
colnames(b_lower) = b_names
colnames(b_upper) = b_names
# rownames(b_mean) = fnames

b_mean_no = as_tibble(b_mean) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Mean")

b_lower_no = as_tibble(b_lower) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Lower")

b_upper_no = as_tibble(b_upper) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "Upper")

colnames(b_no_sig) = b_names
b_no_sig = as_tibble(b_no_sig) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "sig")


b_plt_no = b_mean_no %>% 
  left_join(b_lower_no) %>% 
  left_join(b_upper_no) %>% 
  left_join(b_no_sig) %>% 
  # mutate(sig = as.factor(if_else(sign(Lower) == sign(Upper), 1, 0))) %>% 
  mutate(Variable = factor(Variable, 
                           levels = c('Int', 'Max Depth', 'Lake Area', 'Ag', 'Urban', 'Wetlands', 'GDD', 'Secchi'))) %>% 
  mutate(species = factor(species,
                          levels = c('black crappie', 'bluegill', 'largemouth bass', 'northern pike', 'walleye', 'yellow perch')))
# %>% 
#   mutate(sig = factor(sig))


for(i in 1:nrow(b_plt_no)){
  
  if(between(b_plt_no$Lower[i], b_plt$Lower[i], b_plt$Upper[i])){
    b_plt_no$sig[i] = 0
  }else if(between(b_plt_no$Upper[i], b_plt$Lower[i], b_plt$Upper[i])){
    b_plt_no$sig[i] = 0
  }else if(between(b_plt$Lower[i], b_plt_no$Lower[i], b_plt_no$Upper[i])){
    b_plt_no$sig[i] = 0
  }else if(between(b_plt$Upper[i], b_plt_no$Lower[i], b_plt_no$Upper[i])){
    b_plt_no$sig[i] = 0
  }else{
    b_plt_no$sig[i] = 1
  }
  
}


betas_plt = rbind(b_plt %>% mutate(model = 'with effort scaling'), b_plt_no %>% mutate(model = 'without effort scaling'))

betas_plt = betas_plt %>% 
  mutate(Variable = str_to_lower(Variable)) %>% 
  mutate(Variable = factor(Variable, levels = c('int', 'max depth', 'lake area', 'ag', 'urban', 'wetlands', 'gdd', 'secchi')))

betas_plt = betas_plt %>% 
  mutate(model = factor(model, levels = c('with effort scaling', 'without effort scaling')))


# https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2

library(gtable)
library(cowplot)
library(grid)
library(gridExtra)
library(lemon)

shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}

p1 = ggplot(betas_plt, aes(y = reorder(species, desc(species)))) +
  geom_point(aes(x = Mean, color = model), size = 1.5) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper, color = model), width = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  scale_x_symmetric() +
  facet_wrap(~Variable, scales = 'free_x', ncol = 2,
             labeller = labeller(Variable = c('int' = 'intercept',
                                              'max depth' = 'maximum depth',
                                              'lake area' = 'lake area',
                                              'urban' = 'percent urban',
                                              'wetlands' = 'percent wetlands',
                                              'gdd' = 'growing degree days',
                                              'secchi' = 'secchi disk depth'))) +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual(values = c('blue', 'red')) +
  xlab('') +
  ylab('') +
  guides(shape = 'none') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.height = unit(1, 'cm'),
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_text(angle = 45, vjust = -1, hjust = 1))
grid.draw(shift_legend(p1))


# width = 1700, height = 750
# ggsave('results/spatial_results/parameter_estimate_plot.png', p1, width = 16, height = 8)



p = ggplot(betas_plt %>% filter(Variable != 'int'), aes(y = reorder(species, desc(species)))) +
  geom_point(aes(x = Mean, color = model), size = 1.5) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper, color = model), width = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  scale_x_symmetric() +
  facet_wrap(~Variable, scales = 'free_x', ncol = 2,
             labeller = labeller(Variable = c('max depth' = 'maximum depth',
                                              'lake area' = 'lake area',
                                              'urban' = 'percent urban',
                                              'wetlands' = 'percent wetlands',
                                              'gdd' = 'growing degree days',
                                              'secchi' = 'secchi disk depth'))) +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual(values = c('blue', 'red')) +
  xlab('') +
  ylab('') +
  guides(shape = 'none') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.key.height = unit(1, 'cm'),
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_text(angle = 45, vjust = -1, hjust = 1),
        legend.position="bottom",
        plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))

# ggsave('results/spatial_results/parameter_estimate_plot_no_int.png', width = 10, height = 10)
# ggsave('results/eps_figures/parameter_estimate_plot_no_int.eps', dpi = 600, width = 10, height = 10, device=grDevices::cairo_ps)

beta = t(beta)

colnames(beta) = b_names[2:7]

beta_vals = as_tibble(beta) %>% 
  mutate(species = fnames) %>% 
  pivot_longer(-species, names_to = "Variable", values_to = "truth") %>% 
  mutate(Variable = str_to_lower(Variable)) %>% 
  mutate(Variable = factor(Variable, levels = c('max depth', 'lake area', 'urban', 'wetlands', 'gdd', 'secchi')))

betas_plt_2 = betas_plt %>% 
  left_join(beta_vals, by = c("species", "Variable")) %>% 
  filter(Variable != 'int')

# p = 1
ggplot(betas_plt_2, aes(y = reorder(species, desc(species)))) +
  geom_point(aes(x = Mean, color = model), size = 1.5) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper, color = model), width = 0.2) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_point(aes(x = truth), shape = 8) + 
  scale_x_symmetric() +
  facet_wrap(~Variable, scales = 'free_x', ncol = 2,
             labeller = labeller(Variable = c('max depth' = 'maximum depth',
                                              'lake area' = 'lake area',
                                              'urban' = 'percent urban',
                                              'wetlands' = 'percent wetlands',
                                              'gdd' = 'growing degree days',
                                              'secchi' = 'secchi disk depth'))) +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual(values = c('blue', 'red')) +
  xlab('') +
  ylab('') +
  guides(shape = 'none') +
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.key.height = unit(1, 'cm'),
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_text(angle = 45, vjust = -1, hjust = 1),
        legend.position="bottom",
        plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))



# ggsave('results/parameter_estimate_plot_simulation.png', width = 10, height = 10)
# ggsave('results/parameter_estimate_plot_simulation.eps', dpi = 600, width = 10, height = 10, device=grDevices::cairo_ps)












