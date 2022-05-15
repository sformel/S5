#Plant biomass and stem count over time
#Last updated: April 3, 2022
#By Steve Formel

#load libraries----

library(cowplot)
library(readxl)
library(tidyverse)
library(compositions)
library(vegan)
library(plotrix)
library(brms)

#Break out some elements for easy adjustment of plots----
annotation_size <- 6/ggplot2::.pt
point_size <- 5/ggplot2::.pt
stroke_size <- 0.5/ggplot2::.pt

S5_theme <- theme_bw(base_size = 7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(2,4,2,2, "mm"),
        axis.title.x = element_text(margin = margin(t = 2, 
                                                  r = 0, 
                                                  b = 0, 
                                                  l = 0,
                                  unit = "mm")),
        axis.title.y = element_text(margin = margin(t = 0, 
                                                  r = 2, 
                                                  b = 0, 
                                                  l = 0,
                                  unit = "mm"))) 


#color palette
cPAL <- c("#F0E442", "#0072B2")

#Read in data and clean

data.FP <- "../Relevant Data/data_available_from_GRIIDC/GOMRI_R5x2860000002_data.xlsx"

morph <- read_excel(data.FP,
                    sheet = "plant_morphology", na = "NA")

sample_key <- read_excel(data.FP,
                         sheet = "sample_key", na = "NA")

df <- right_join(x = sample_key, y = morph, by = 'sampleID_stem')

#fix factors
df <- df %>%
  mutate(across(where(is.character), as.factor)) %>% 
  dplyr::mutate(across(c(plantID, table, block), factor)) %>%
  dplyr::mutate(across(c(num_stems_live, num_nodes), as.integer)) %>%
  mutate(sampling_period = factor(sampling_period, 
                                  levels = c("N16", "J17", "N17", "J18"))) %>% 
  filter(sampling_period!="M16")


#Gussy up the factor levels
df$oil_added <- plyr::revalue(df$oil_added, c("Y" = "Oil Added", "N" = "No Oil Added"))
df$orig_soil <- plyr::revalue(df$orig_soil, c("Y" = "Prev-Oiled Inoc.", "N" = "Not Prev-Oiled Inoc"))

  
## Model Num of live stems

m <- brm(num_stems_live ~ 1 + (1|plantID) + sampling_period*orig_soil*oil_added, 
         data = df,
         family = "negbinomial",
         control = list(adapt_delta = 0.99, 
                        max_treedepth = 10),
         cores = 4)

summary(m)

#These are commented out so they don't hinder plot generation when run as a script
#pp_check(object = m, type = "dens_overlay", nsamples = 100)
#pp_check(m, type = "stat", stat = 'median', nsamples = 100)
#pp_check(m, type = "stat", stat = 'mean', nsamples = 100)
#pp_check(m ,type = 'intervals', nsamples = 100)

#conditional_effects(m)

## plot
p <- ggplot(df,
       aes(x = unclass(sampling_period),
           y = num_stems_live,
           fill = orig_soil,
           shape = oil_added)) +
  geom_point(size = point_size,
             alpha = 0.5,
             position = position_dodge(width = 0.4)) +
  S5_theme +
  labs(y = "Number of Live Stems",
       x = "Time",
       fill = "Soil Inoculum",
       color = "Soil Inoculum",
       shape = "Oil Addition") +
  theme(legend.position = "right") +
  theme(legend.position = "bottom", 
        legend.box="horizontal",
        legend.margin=margin(t = 1, r = 1, b = 1, l = 1,
                             unit='mm'),
        legend.title.align=0.5,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "lightgray",
                                             size = stroke_size),
        legend.key.width = unit(5, "mm")) +
  geom_smooth(aes(color = orig_soil),
              formula = y ~ log(x), 
              method = "glm",
              size = stroke_size + 0.5,
              se = FALSE) +
  scale_fill_manual(values = cPAL) +
  scale_shape_manual(values = c(21,24)) +
  scale_color_manual(values = cPAL) +
  scale_x_continuous(breaks = c(1:5),
                     labels = c("Start", "6 Months", "1 Year", "1.5 Years", "2 Years")) +
    guides(fill=guide_legend(title.position = "top",
                           nrow = 2,
                           override.aes=list(shape=21,
                                             size = stroke_size + 0.5)),
         shape=guide_legend(title.position = "top",
                            nrow = 2,
                            override.aes=list(size = point_size)))

# #This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

p

ggsave(filename = "S5_figure3.pdf",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#For draft and co-author review
ggsave(filename = "S5_figure3.png",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#convert pdf to TIFF for publisher
library(pdftools)

pdf_convert(pdf = "figures/S5_figure3.pdf",
format = "tiff",
filenames = "figures/S5_figure3.tif",
dpi = 300)


## Poisson on final stem count (J18)

#Following this example: https://www.theanalysisfactor.com/generalized-linear-models-in-r-part-6-poisson-regression-count-variables/
#https://www.theanalysisfactor.com/glm-r-overdispersion-count-regression/
#Also read this about the scale parameter in quasi-poisson. 
#http://dept.stat.lsa.umich.edu/~kshedden/Courses/Stat504/posts/glm/

hist(df$num_stems_live[df$sampling_period=="J18"])
  
m <- glm(num_stems_live ~ orig_soil*oil_added, 
         data = df %>% 
           filter(sampling_period=="J18"),
         family = poisson)

summary(m) 

#Overdispersed because Residual deviance >> degrees of freedom

#Effect size
df %>% 
  filter(sampling_period=="J18") %>% 
  group_by(oil_added) %>% 
  dplyr::summarise(mean = mean(num_stems_live),
                   se = std.error(num_stems_live))

### Quasipoisson

#See this for interpretation: https://stats.stackexchange.com/questions/470640/interpreting-interaction-among-2-categorical-iv-in-quasi-poisson-regression

m <- glm(num_stems_live ~ orig_soil*oil_added, 
         data = df %>% 
           filter(sampling_period=="J18"),
         family = quasipoisson(link = "log"))

summary(m) 

#F-test is most appropriate for quasipoisson according to the help for anova.glm
anova(m, test = "F")


#Not Significant Traits------

#Number of nodes

ggplot(df,
       aes(x = unclass(sampling_period),
           y = num_nodes,
           fill = orig_soil,
           shape = oil_added)) +
  geom_point(size = 2,
             alpha = 0.5,
             position = position_dodge(width = 0.4)) +
  theme_bw() +
  labs(y = "Mean Nodes per Stem",
       x = "Time",
       fill = "Inoculum",
       shape = "Oil Addition") +
  theme(legend.position = "right") +
  guides(fill = guide_legend(override.aes=list(shape=21)),
         color = guide_legend(override.aes=list(shape=21)),
         shape = guide_legend(override.aes=list(size = 4))) +
  scale_fill_manual(values = cPAL) +
  scale_shape_manual(values = c(21,24)) +
  scale_color_manual(values = cPAL)

#It's clear I did a bad job communicating how to count nodes.  Can't use this data.  Also clear that there isn't that big of a difference within any time point.


#Stem Height
ggplot(df,
       aes(x = unclass(sampling_period),
           y = stem_ht,
           fill = orig_soil,
           shape = oil_added)) +
  geom_point(size = 2,
             alpha = 0.5,
             position = position_dodge(width = 0.4)) +
  theme_bw() +
  labs(y = "Mean Stem Height (cm)",
       x = "Time",
       fill = "Inoculum",
       shape = "Oil Addition") +
  theme(legend.position = "right") +
  guides(fill = guide_legend(override.aes=list(shape=21)),
         color = guide_legend(override.aes=list(shape=21)),
         shape = guide_legend(override.aes=list(size = 4))) +
  scale_fill_manual(values = cPAL) +
  scale_shape_manual(values = c(21,24)) +
  scale_color_manual(values = cPAL)

## Model
m <- brm(stem_ht ~ 1 + (1|plantID) + unclass(sampling_period) + orig_soil*oil_added, 
         data = df,
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

#Periods 1 + 3 (November) are taller than 2 + 4, so that makes sense.  It doesn't look like there is anything worth pursuing here.  Stems may have gotten shorter over time, but it doesn't look like it's unique to any treatment.

#Stem Diameter
ggplot(df,
       aes(x = unclass(sampling_period),
           y = stem_diam,
           fill = orig_soil,
           shape = oil_added)) +
  geom_point(size = 2,
             alpha = 0.5,
             position = position_dodge(width = 0.4)) +
  theme_bw() +
  labs(y = "Mean Stem Height (mm)",
       x = "Time",
       fill = "Inoculum",
       shape = "Oil Addition") +
  theme(legend.position = "right") +
  guides(fill = guide_legend(override.aes=list(shape=21)),
         color = guide_legend(override.aes=list(shape=21)),
         shape = guide_legend(override.aes=list(size = 4))) +
  scale_fill_manual(values = cPAL) +
  scale_shape_manual(values = c(21,24)) +
  scale_color_manual(values = cPAL)

## Model
m <- brm(stem_diam ~ 1 + (1|plantID) + unclass(sampling_period) + orig_soil*oil_added, 
         data = df,
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

#Like height, it doesn't look like there is any effect of the treatments.  Diameters may get slightly smaller in J18, but that could also be a function of who was measuring them (if they squeezed the stems, or measured height poorly).  Still I'll plot the means of this one to see if there is a difference that corresponds to the growth chamber.

## plot stem diam means
df %>%
  group_by(sampling_period, oil_added, orig_soil) %>% 
  summarise_at("stem_diam",
               funs(mean,
                     sd, 
                     se=sd(.)/sqrt(n()))) %>% 
  ggplot(aes(x = sampling_period,
           y = mean,
           fill = oil_added,
           shape = orig_soil))  +
  geom_point(size = 4,
             alpha = 0.5,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(min = mean - 2*se,
                    max = mean + 2*se),
                width = 0.2,
                position = position_dodge(width = 0.4)) +
  labs(y = "Mean Stem Diameter (mm)",
       x = NULL,
       fill = "Oil Addition", 
       caption = "Error bars represent \U00B1 2 SE") +
  scale_fill_manual(values = c("white", "black")) +
  scale_shape_manual(values = c(21,24)) +
  theme_bw()  +
  theme(legend.position = "right") +
  guides(fill=guide_legend(override.aes=list(shape=21)))