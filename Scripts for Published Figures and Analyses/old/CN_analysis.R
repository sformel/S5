#Script for S5

#Description: Compare C:N among treatments

#Last updated Jan 20, 2020 by Steve Formel

#import data
source("../ITS/import_and_clean_S5_ITS.R")

#make subtitle for metadata on all figures-----
subtitle <- "Exploratory Analysis, Script: CN_analysis.R"

#as of Jan 20, 2020 I'm still waiting on CN for soil from Pardue lab

#look at CN over time----
library(ggplot2)

#Aesthetics

#color palette
cPAL <- c("#E69F00", "#0072B2")


#subset to CN measurements only

S5 <- prune_samples(sample_data(S5)$tissue %in% c("Leaf", "Root") &
                      sample_data(S5)$sampling_period!="NA", S5)

#subset to traits that have values for all four sampling periods
df <- sample_data(S5) %>% 
  data.frame(df)

df.traits.all <- df %>%
  select(c(sampleID_ITS, sampling_period, plantID, oil_added, orig_soil, tissue:C.N)) %>% 
  gather(key = Trait, value = value, percent_N:C.N)

#make avg diversity
library(Rmisc)

sum.val <- summarySE(df.traits.all, measurevar = "value", groupvars = c("Trait", "tissue", "oil_added", "orig_soil", "sampling_period"))

ggplot(data = na.omit(sum.val), 
       aes(x = sampling_period,
           y = value,
           shape = oil_added,
           fill = orig_soil)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = value-2*se,ymax = value+2*se), width= 0.2) +
  facet_grid(cols = vars(tissue), 
             rows = vars(Trait),
             scales = "free") +
  labs(title = "S5 CN leaves and roots",
       caption = "Error Bars represent ? 2 SE") +
       theme(legend.position = "right") +
  guides(fill = guide_legend(override.aes=list(shape=21)),
         color = guide_legend(override.aes=list(shape=21)),
         shape = guide_legend(override.aes=list(size = 4))) +
  scale_fill_manual(values = cPAL) +
  scale_shape_manual(values = c(21,24)) +
  scale_color_manual(values = cPAL)

ggsave("figs/exploratory/CN/S5_CN_exploration_all.png", width = 11, height = 7, units = "in")

## Bayesian MM

library(brms)

df.traits.all <- df %>%
  select(c(sampleID_ITS, sampling_period, plantID, oil_added, orig_soil, tissue:C.N))

#percent C

m <- brm(percent_C ~ 1 + (1|plantID) + sampling_period*tissue*orig_soil*oil_added, 
         data =  df.traits.all %>% 
           filter(tissue=="Root"),
         family = "skew_normal",
         control = list(adapt_delta = 0.99, 
                        max_treedepth = 10),
         cores = 4)

summary(m)

pp_check(object = m, type = "dens_overlay", nsamples = 100)
pp_check(m, type = "stat", stat = 'median', nsamples = 100)
pp_check(m, type = "stat", stat = 'mean', nsamples = 100)
pp_check(m ,type = 'intervals', nsamples = 100)

conditional_effects(m)

#percent N

m <- brm(percent_N ~ 1 + (1|plantID) + sampling_period*orig_soil*oil_added, 
         data = df.traits.all %>% 
           filter(tissue=="Root"),
         family = "skew_normal",
         control = list(adapt_delta = 0.99, 
                        max_treedepth = 15),
         cores = 4)

summary(m)

pp_check(object = m, type = "dens_overlay", nsamples = 100)
pp_check(m, type = "stat", stat = 'median', nsamples = 100)
pp_check(m, type = "stat", stat = 'mean', nsamples = 100)
pp_check(m ,type = 'intervals', nsamples = 100)

conditional_effects(m)