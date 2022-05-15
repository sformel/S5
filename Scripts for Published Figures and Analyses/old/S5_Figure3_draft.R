#PAH Decomposition over time
#Last updated: March 24, 2022
#By Steve Formel

#load libraries----

library(dplyr)
library(cowplot)
library(readxl)
library(tidyverse)
library(compositions)
library(plyr)
library(vegan)
library(plotrix)

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

#CLR vs Total PAH for N16 and J18 only, full set of samples
#remove S5_60_J18, which had something funky happen in analysis per Vijai

#Figure 3: Full set of N16 and J18 samples only

oil <- read_excel("../data/GOMRI_R5x2860000002_supplementary_oil_data.xlsx", 
                 sheet = "data", na = "NA")

sample_key <- read_excel("../data/uploaded_to_GRIIDC/GOMRI_R5x2860000002_data.xlsx", 
                                     sheet = "sample_key", na = "NA")

oil <- right_join(x = sample_key, 
                  y = oil, by = 'sampleID_stem')

#Clean Data
oil <- oil %>%
  select(-plantID.y) %>% 
  rename(c("plantID.x" = "plantID")) %>%
  dplyr::mutate(across(where(is.character), factor)) %>% 
  dplyr::mutate(across(c(plantID, table, block), factor)) %>% 
  mutate(sampling_period = factor(sampling_period, 
                                  levels = c("N16", "J17", "N17", "J18")))

#Select and summarize PAHs
oil <- oil %>% 
  mutate(tot_naph = `C1-naphthalenes`+`C2-naphthalenes`+`C3-naphthalenes`+`C4-naphthalenes`,
         tot_phen = `C1-phenanthrenes`+`C2-phenanthrenes`+`C3-phenanthrenes`+`C4-phenanthrenes`,
         tot_dibenz = `C1-dibenzothiophenes`+`C2-dibenzothiophenes`+`C3-dibenzothiophenes`,
         tot_chrys = `C1-chrysenes`+`C2-chrysenes`+`C3-chrysenes`,
         tot_PAH = tot_naph + tot_phen + tot_dibenz + tot_chrys) %>% 
  select(sampleID_stem, 
         plantID, 
         plant_trt, 
         orig_soil, 
         oil_added, 
         sampling_period, 
         tot_PAH, 
         tot_naph, 
         tot_phen, 
         tot_dibenz, 
         tot_chrys)  

#average if there were multiple samples per time period
oil <- oil %>% 
  group_by(sampleID_stem, 
         plantID, 
         plant_trt, 
         orig_soil, 
         oil_added, 
         sampling_period) %>% 
  dplyr::summarise(across(tot_PAH:tot_chrys, mean),
                   .groups = 'drop')

#Gussy up the factor levels
oil$oil_added <- revalue(oil$oil_added, c("Y" = "Oiled", "N" = "No Oil Added")) %>% 
  as.factor()
oil$orig_soil <- revalue(oil$orig_soil, c("Y" = "Prev-Oiled Inoculum", "N" = "Not Prev-Oiled Inoculum")) %>% 
  as.factor()
oil$plant_trt <- revalue(oil$plant_trt, c("Y" = "Plant", "N" = "No Plant")) %>% 
  as.factor()

#note, can't use mutate(across()) here because CLR needs to work on the columns as a matrix, not indivudal columns

clr.PAH <- oil %>% 
  filter(oil_added=="Oiled",
         sampling_period %in% c("N16","J18"), 
         !(plantID==60 & sampling_period=="J18")) %>% 
  select(tot_naph, 
         tot_dibenz, 
         tot_phen, 
         tot_chrys) %>% 
  na.omit() %>% 
  clr() %>% 
  as.data.frame()

colnames(clr.PAH) <- paste("CLR", colnames(clr.PAH), sep = "_")

#bind
df.clr <- cbind(oil %>%
                  filter(oil_added=="Oiled",
                         sampling_period %in% c("N16","J18"),
                         !(plantID==60 & sampling_period=="J18")), 
                clr.PAH)

#plot

p <- df.clr %>%
  gather(key = "PAH", 
         value = "clr_val", 
         CLR_tot_naph: CLR_tot_chrys) %>%
  ggplot(aes(x = clr_val,
             y = tot_PAH,
             fill = PAH)) +
  geom_jitter(aes(shape = plant_trt),
              size = point_size,
              stroke = stroke_size,
              alpha = 0.5,
              color = "black",
              width = point_size) +
  scale_shape_manual(values = c(21,
                                24)) +
  scale_fill_manual(values = c("red", "gray", "black", "white"),
                    labels = c("chrysenes", 
                               "dibenzothiophenes", 
                               "naphthalenes", 
                               "phenanthrenes")) +
  labs(shape = "Plant Treatment",
       fill = "PAH",
       x = "Centered Log-Ratio",
       y = "Total PAHs (ug/g)") +
  S5_theme  +
  theme(legend.position = "bottom", 
        legend.box="horizontal",
        legend.margin=margin(t = 1, r = 1, b = 1, l = 1,
                             unit='mm'),
        legend.title.align=0.5,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "lightgray",
                                             size = stroke_size)) +
  guides(fill=guide_legend(title.position = "top",
                           nrow = 2,
                           override.aes=list(shape=21,
                                             size = annotation_size + 1)),
         shape=guide_legend(title.position = "top",
                            nrow = 2,
                            override.aes=list(size = annotation_size + 1))) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(sampling_period = c("N16" = "6 months", "J18" = "2 years"))) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  scale_y_log10()

p

# #This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

ggsave(filename = "S5_Figure3.pdf",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#For draft and co-author review
ggsave(filename = "S5_Figure3.png",
       path = "figures/", 
       width = 85,
       height = 85,
       units = "mm", 
       dpi = 300,
       scale = 1)


#Plot points as averages with error bars

oil.avg <- df.clr %>% 
  gather(key = "PAH", 
         value = "Value", 
         c(CLR_tot_naph:CLR_tot_chrys)) %>% 
  group_by(plant_trt, 
         sampling_period,
         PAH) %>% 
  dplyr::summarise(across(c(Value), 
                          list(mean = mean,
                               se = std.error)),
                   .groups = 'drop') %>% 
  right_join(., df.clr %>%
               group_by(plant_trt,
                        sampling_period) %>%
               dplyr::summarise(across(c(tot_PAH),
                                       list(mean = mean,
                                            se = std.error))), 
             .groups = 'drop')

#plot

p <- oil.avg %>%
  ggplot(aes(x = Value_mean,
             y = tot_PAH_mean,
             fill = PAH)) +
  geom_point(aes(shape = plant_trt),
              size = point_size,
              stroke = stroke_size,
              alpha = 0.5,
              color = "black") +
  geom_errorbar(aes(ymin = tot_PAH_mean - tot_PAH_se,
                    ymax = tot_PAH_mean + tot_PAH_se),
                size = stroke_size) +
  geom_errorbarh(aes(xmin = Value_mean - Value_se,
                    xmax = Value_mean + Value_se),
                 size = stroke_size) +
  scale_shape_manual(values = c(21,
                                24)) +
  scale_fill_manual(values = c("red", "gray", "black", "white"),
                    labels = c("chrysenes", 
                               "dibenzothiophenes", 
                               "naphthalenes", 
                               "phenanthrenes")) +
  labs(shape = "Plant Treatment",
       fill = "PAH",
       x = "Centered Log-Ratio",
       y = "Total PAHs (ug/g)") +
  S5_theme  +
  theme(legend.position = "bottom", 
        legend.box="horizontal",
        legend.margin=margin(t = 1, r = 1, b = 1, l = 1,
                             unit='mm'),
        legend.title.align=0.5,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "lightgray",
                                             size = stroke_size)) +
  guides(fill=guide_legend(title.position = "top",
                           nrow = 2,
                           override.aes=list(shape=21,
                                             size = annotation_size + 1)),
         shape=guide_legend(title.position = "top",
                            nrow = 2,
                            override.aes=list(size = annotation_size + 1))) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(sampling_period = c("N16" = "6 months", "J18" = "2 years"))) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  scale_y_log10()

p


# #This was extremely helpful for tweaking the plot
nflplotR::ggpreview(plot = p,
                    width = 85,
                    height = 85,
                    units = "mm",
                    dpi = 300,
                    scale = 1,
                    device = "pdf")


#point cloud with mean as bright point

p <- ggplot() +
  geom_jitter(data = df.clr %>%
                gather(key = "PAH",
                       value = "clr_val",
                       CLR_tot_naph: CLR_tot_chrys),
              inherit.aes = FALSE,
              aes(x = clr_val,
             y = tot_PAH,
             fill = PAH,
             shape = plant_trt),
              size = point_size,
              stroke = stroke_size,
              alpha = 0.2,
              color = "black",
              width = point_size) +
  geom_point(data = oil.avg,
             inherit.aes = FALSE,
             aes(x = Value_mean,
             y = tot_PAH_mean,
             fill = PAH,
             shape = plant_trt),
              size = point_size + 1,
              stroke = stroke_size,
              alpha = 1,
              color = "black") +
  scale_shape_manual(values = c(21,
                                24)) +
  scale_fill_manual(values = c("red", "gray", "black", "white"),
                    labels = c("chrysenes", 
                               "dibenzothiophenes", 
                               "naphthalenes", 
                               "phenanthrenes")) +
  labs(shape = "Plant Treatment",
       fill = "PAH",
       x = "Centered Log-Ratio",
       y = "Total PAHs (ug/g)") +
  S5_theme  +
  theme(legend.position = "bottom", 
        legend.box="horizontal",
        legend.margin=margin(t = 1, r = 1, b = 1, l = 1,
                             unit='mm'),
        legend.title.align=0.5,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "lightgray",
                                             size = stroke_size)) +
  guides(fill=guide_legend(title.position = "top",
                           nrow = 2,
                           override.aes=list(shape=21,
                                             size = annotation_size + 1)),
         shape=guide_legend(title.position = "top",
                            nrow = 2,
                            override.aes=list(size = annotation_size + 1))) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(sampling_period = c("N16" = "6 months", "J18" = "2 years"))) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  scale_y_log10()

p

# #This was extremely helpful for tweaking the plot
nflplotR::ggpreview(plot = p,
                    width = 85,
                    height = 85,
                    units = "mm",
                    dpi = 300,
                    scale = 1,
                    device = "pdf")

#Load in a clean data of balanced set of samples-----

oil <- read_excel("../data/uploaded_to_GRIIDC/GOMRI_R5x2860000002_data.xlsx", 
                 sheet = "oil", na = "NA")

sample_key <- read_excel("../data/uploaded_to_GRIIDC/GOMRI_R5x2860000002_data.xlsx", 
                                     sheet = "sample_key", na = "NA")

df <- right_join(x = sample_key, y = oil, by = 'sampleID_stem')

#fix factors
df <- df %>%
  mutate_if(is.character, as.factor)

df$plantID <- as.factor(df$plantID)
df$table <- as.factor(df$table)
df$block <- as.factor(df$block)

df$sampling_period <- factor(df$sampling_period, levels = c("N16", "J17", "N17", "J18"))

#Gussy up the factor levels
df$oil_added <- revalue(df$oil_added, c("Y" = "Oiled", "N" = "No Oil Added")) %>% 
  as.factor()
df$orig_soil <- revalue(df$orig_soil, c("Y" = "Prev-Oiled Inoculum", "N" = "Not Prev-Oiled Inoculum")) %>% 
  as.factor()
df$plant_trt <- revalue(df$plant_trt, c("Y" = "Plant", "N" = "No Plant")) %>% 
  as.factor()

# PAHs in Soil----

df.PAHcomp <- df %>% 
  select(plantID,
         plant_trt,
         orig_soil,
         oil_added,
         sampling_period,
         tissue,
         table,
         block,
         total_methy_naph,
         total_methy_dibenz,
         total_methy_phen,
         total_methy_chry,
         total_relevant_PAHs) %>% 
  filter(tissue=="Soil" & oil_added=="Oiled")
         
## CLR Transformation

clr.PAH <- clr(na.omit(df.PAHcomp[,c(9:12)])) %>% 
  data.frame()

#PCA as visual companion to PERMANOVA
pc <- princomp(x = clr.PAH)

df.pca <- pc$scores

cbind(na.omit(df.PAHcomp), df.pca) %>%
  ggplot(aes(x = Comp.1,
             y = Comp.2,
             fill = orig_soil,
             shape = plant_trt)) +
  geom_point(size = 4,
             stroke = 1,
             alpha = 0.5,
             color = "black") +
  scale_fill_manual(values = cPAL) +
  scale_shape_manual(values = c(21,24)) +
  labs(fill = "Inoculum",
       shape = "Plant Presence") +
  theme_bw() +
  theme(legend.position = "right") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  facet_wrap(~ sampling_period)

## PERMANOVA

clr.PAH <- clr(na.omit(df.PAHcomp[,c(9:12)])) %>% 
  as.data.frame()

df.clr <- cbind(df.PAHcomp[,1:8], clr.PAH)

set.seed(1)
adonis(formula = clr.PAH ~ plant_trt*orig_soil, 
       data = df.clr, 
       permutations = 9999, 
       method = "euclidean",
       strata = df.clr$sampling_period)

adonis2(formula = clr.PAH ~ plant_trt*orig_soil, 
       data = df.clr, 
       permutations = 9999, 
       method = "euclidean",
       strata = df.clr$sampling_period,
       by = "margin")

# Inoculum had an effect.  Not a surprise given that PAHs came in with the inocula. Plant trt didn't have an effect. Time didn't interact with either. In general the PCA makes it obvious that the PAH data got noisier over time regardless of treatment.  More evidence that this is difficult to measure.

df.clr %>%
  gather(key = "PAH", value = "clr_val", c(9:12)) %>%
  ggplot(aes(x = PAH,
         y = clr_val),
         color = "black") +
  geom_jitter(aes(fill = orig_soil,
                  shape = plant_trt),
              size = 4,
              stroke = 1,
              alpha = 0.5,
              width = 0.2) +
  scale_fill_manual(values = cPAL) +
  scale_shape_manual(values = c(21,24)) +
  labs(x = NULL,
       fill = "Inoculum",
       shape = "Plant Presence") +
  theme_bw() +
  theme(legend.position = "bottom", legend.box="vertical") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  facet_wrap(~ sampling_period,
             labeller = labeller(sampling_period = c("N16" = "6 Months", "J17" = "1 Year", "N17" = "1.5 Years", "J18" = "2 Years"))) +
  scale_x_discrete(labels = c("Chrysenes", "Dibenzothiophenes", "Naphthalenes", "Phenanthrenes")) +
  labs(y = "Centered Log-Ratio") +
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0,
             linetype = "dashed")


### CLR and total PAH on same plot

df.clr$total_PAH <- df.PAHcomp[df.PAHcomp$oil_added=="Oiled",]$total_relevant_PAHs

## GLM on total PAHs

library(lme4)
m <- lmer(scale(total_PAH) ~ orig_soil*plant_trt + (1|sampling_period), 
         data = df.clr,
         REML = F)

summary(m) 

confint(m)


set.seed(1)
m <- adonis(formula = clr.PAH ~ orig_soil*plant_trt, data = df.clr, permutations = 9999, method = "euclidean",
            strata = df.clr$sampling_period)

m


### Test change in PAH over time

C <- B %>% 
  filter(sampleID_stem!="S5_60_J18") %>% 
  select(plantID, plant_trt, orig_soil, oil_added, sampling_period, tot_PAH) %>% 
  group_by(plantID,sampling_period) %>% filter(duplicated(plantID) | n()==1)

df.delta <- C %>% 
  filter(oil_added=="Oiled") %>% 
  spread(key = sampling_period, value = tot_PAH) %>% 
  mutate("deltaPAH" = J18 - N16) %>% 
  na.omit()

m <- lm(deltaPAH ~ plant_trt*orig_soil, data = df.delta)

summary(m)
hist(resid(m))
shapiro.test(resid(m))

kruskal.test(x = df.delta$deltaPAH, g = df.delta$plant_trt)


#What about just at J18?
D <- C %>% 
  filter(sampling_period=="J18") %>% 
  filter(oil_added=="Oiled")

kruskal.test(x = D$tot_PAH, g = D$plant_trt)


m <- lm(log(tot_PAH) ~ sampling_period*plant_trt*orig_soil, data = C %>% 
          filter(oil_added=="Oiled"))
summary(m)
hist(resid(m))
shapiro.test(resid(m))

car::Anova(m, type = "III")

```

### Plot

```{r}

ggplot(df.delta,
       aes(x = orig_soil,
           y = -deltaPAH,
           fill = orig_soil)) +
  geom_jitter(alpha = 0.5,
             color = "black",
             width = 0.1,
             size = 4,
             shape = 21) +
  scale_fill_manual(values = c("white", "black")) +
  labs(y = "\u0394 Total PAHs (ug/g)",
       fill = "Oil Addition") +
  theme_bw()  +
  theme(legend.position = "right") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  facet_grid(cols = vars(plant_trt)) +
  scale_y_log10()
```

```{r}

# Just J18

C %>% 
  filter(sampling_period=="J18") %>% 
  filter(oil_added=="Oiled") %>% 
  ggplot(aes(x = orig_soil,
           y = tot_PAH,
           fill = orig_soil)) +
  geom_jitter(alpha = 0.5,
             color = "black",
             width = 0.1,
             size = 4,
             shape = 21) +
  scale_fill_manual(values = c("white", "black")) +
  labs(y = "\u0394 Total PAHs (ug/g)",
       fill = "Oil Addition") +
  theme_bw()  +
  theme(legend.position = "right") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  facet_grid(cols = vars(plant_trt))
```

```{r}
#What about N16?

C %>% 
  filter(sampling_period=="N16") %>% 
  filter(oil_added=="Oiled") %>% 
  ggplot(aes(x = orig_soil,
           y = tot_PAH,
           fill = orig_soil)) +
  geom_jitter(alpha = 0.5,
             color = "black",
             width = 0.1,
             size = 4,
             shape = 21) +
  scale_fill_manual(values = c("white", "black")) +
  labs(y = "\u0394 Total PAHs (ug/g)",
       fill = "Oil Addition") +
  theme_bw()  +
  theme(legend.position = "right") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  facet_grid(cols = vars(plant_trt))

```

```{r}
#Together

C %>% 
  filter(oil_added=="Oiled") %>% 
  ggplot(aes(x = plant_trt,
           y = tot_PAH,
           fill = orig_soil)) +
  geom_jitter(alpha = 0.5,
             color = "black",
             width = 0.1,
             size = 4,
             shape = 21) +
  scale_fill_manual(values = cPAL) +
  labs(x = NULL,
       y = "Total PAHs (ug/g)",
       fill = "Inoculum") +
  theme_bw()  +
  theme(legend.position = "right") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(sampling_period = c("N16" = "6 months", "J18" = "2 years"))) +
  scale_y_log10()

ggsave("S5_total_PAH_extra_samples.png", width = 6, height = 5)
```

```{r}
D <- C %>% 
  filter(oil_added=="Oiled") %>% 
  filter(sampling_period=="J18") 

hist(D$tot_PAH)

#probably a GLM with gamma dist is more appropriate

m <- glm(tot_PAH ~ plant_trt*orig_soil*sampling_period,
        data = C %>% 
  filter(oil_added=="Oiled"),
        family = Gamma(link = "inverse"))

anova(m, test = "F")
summary(m)
hist(resid(m))
shapiro.test(resid(m))



kruskal.test(D$tot_PAH, D$plant_trt)
kruskal.test(D$tot_PAH, D$orig_soil)

#What about effect size?
library(Rmisc)
summarySE(D, measurevar = "tot_PAH", groupvars = c("plant_trt"))
```
