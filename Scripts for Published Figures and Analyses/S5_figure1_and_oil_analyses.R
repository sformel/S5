#PAH Decomposition over time
#Last updated: May 14, 2022
#By Steve Formel

#load libraries----

library(cowplot)
library(readxl)
library(tidyverse)
library(compositions)
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

oil <- read_excel("../Relevant Data/data_available_from_GRIIDC/GOMRI_R5x2860000002_supplementary_oil_data.xlsx", 
                 sheet = "data", na = "NA")

sample_key <- read_excel("../Relevant Data/data_available_from_GRIIDC/GOMRI_R5x2860000002_data.xlsx", 
                                     sheet = "sample_key", na = "NA")

oil <- right_join(x = sample_key, 
                  y = oil, by = 'sampleID_stem')

#Clean Data
oil <- oil %>%
  rename(c("plantID" = "plantID.x")) %>%
  select(-plantID.y) %>% 
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
oil$oil_added <- plyr::revalue(oil$oil_added, c("Y" = "Oiled", "N" = "No Oil Added")) %>% 
  as.factor()
oil$orig_soil <- plyr::revalue(oil$orig_soil, c("Y" = "Prev-Oiled Inoculum", "N" = "Not Prev-Oiled Inoculum")) %>% 
  as.factor()
oil$plant_trt <- plyr::revalue(oil$plant_trt, c("Y" = "Plant", "N" = "No Plant")) %>% 
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

ggsave(filename = "S5_figure1.pdf",
       path = "figures/",
       width = 85,
       height = 85,
       units = "mm",
       dpi = 300,
       scale = 1)

#For draft and co-author review
ggsave(filename = "S5_figure1.png",
       path = "figures/", 
       width = 85,
       height = 85,
       units = "mm", 
       dpi = 300,
       scale = 1)

#convert pdf to TIFF for publisher
library(pdftools)

pdf_convert(pdf = "figures/S5_figure1.pdf",
format = "tiff",
filenames = "figures/S5_figure1.tif",
dpi = 300)



##Stats on beginning and end

#PCA as visual companion to PERMANOVA
pc <- princomp(x = clr.PAH)

df.pca <- pc$scores

cbind(na.omit(df.clr), df.pca) %>%
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
#Strata argument
#https://github.com/vegandevs/vegan/issues/427

perm <- how(nperm = 9999)
setBlocks(perm) <- with(df.clr, sampling_period)

set.seed(1)
adonis(formula = clr.PAH ~ plant_trt*orig_soil, 
       data = df.clr, 
       permutations = perm, 
       method = "euclidean")

adonis2(formula = clr.PAH ~ plant_trt*orig_soil, 
       data = df.clr, 
       method = "euclidean",
       permutations = perm,
       by = "margin")


### ANOVA on difference in PAHs from N16 to J18 in PAHs----

oil.anova <- oil %>% 
          filter(oil_added=="Oiled", 
                 sampleID_stem!="S5_60_J18",
                 sampling_period=="J18") 

m <- lm(sqrt(tot_PAH) ~ plant_trt*orig_soil, data = oil.anova)
summary(m)
hist(resid(m))
shapiro.test(resid(m))

car::Anova(m, type = "III")


#mean and se
oil.anova %>% 
  group_by(plant_trt, sampling_period) %>% 
  summarise(across(tot_PAH, list(mean = mean,
                          se = std.error)))
