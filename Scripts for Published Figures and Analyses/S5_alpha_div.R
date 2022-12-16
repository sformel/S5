#S5 Project Alpha Diversity
#Last updated: 2022-12-16

## Description
#Assess alpha diversity of 16S and ITS communities

## Load Libraries

library(readxl)
library(plyr)
library(dplyr)
library(phyloseq)
library(data.table)
library(reshape2)
library(car)
library(cowplot)
library(vegan)

## Set up aesthetic for plotting

#Set up plotting objects

#Break out some elements for easy adjustment of plots----
annotation_size <- 6 / ggplot2::.pt
point_size <- 5 / ggplot2::.pt
stroke_size <- 0.25 / ggplot2::.pt

#color palette
cPAL <- c("#F0E442", "#0072B2")

S5_theme <- theme_bw(base_size = 7) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(7, 4, 2, 2, "mm"),
    axis.title.x = element_text(margin = margin(
      t = 2,
      r = 0,
      b = 0,
      l = 0,
      unit = "mm"
    )),
    axis.title.y = element_text(margin = margin(
      t = 0,
      r = 2,
      b = 0,
      l = 0,
      unit = "mm"
    ))
  )


## Import 16S phyloseq object and clean

#read in data

S5 <- readRDS("../Relevant Data/feature_tables_as_phyloseq_objects/S5_16S_pseq.rds")

taxa_names(S5) <- paste0("Seq", seq(ntaxa(S5)))

S5.noCL.noMT <- S5 %>%
  subset_taxa(Family != "Mitochondria" & Order != "Chloroplast")

S5.16S <- S5.noCL.noMT

sort(sample_sums(S5.16S))

S5 <- S5.16S

sample_data(S5)$seq_count <- sample_sums(S5)

#fix factors
sample_data(S5) <- sample_data(S5) %>%
  data.frame() %>%
  mutate(across(where(is.character), as.factor)) %>%
  dplyr::mutate(across(c(plantID, table, block), factor))

#Decided to discard all 16S but soil from analyses. Not enough samples remained.
S5 <- subset_samples(sample_data(S5)$tissue=="Soil", physeq = S5) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

#Find minimal slope for rarecurve without dropping too many samples

A <- lapply(sort(sample_sums(S5))[c(1:10)], function(x){
  format(median(rareslope(otu_table(S5), sample = x)), scientific = FALSE)
})

B <- lapply(sort(sample_sums(S5))[c(1:10)], function(x){
  format(mean(rareslope(otu_table(S5), sample = x)), scientific = FALSE)
})

C <- lapply(sort(sample_sums(S5))[c(1:10)], function(x){
  sum(sort(sample_sums(S5))[c(1:10)] < x)
})

D <- lapply(sort(sample_sums(S5))[c(1:10)], function(x){
  print(x)
  })

dropdf <- data.frame(median = unlist(A),
                     mean = unlist(B),
                     num_dropped = unlist(C),
                     num_seqs = unlist(D))  

dropdf

#Plot Shannon diversity for each tissue type by time point.
S5.rare <- rarefy_even_depth(physeq = S5,
                              sample.size = 10441)
S5 <- S5.rare
sort(sample_sums(S5))

#Soil
p <- plot_richness(S5, measures = "Shannon")

#All points
p$data %>%
  ggplot(aes(x = sampling_period,
           y = value,
           shape = oil_added))  +
  geom_point(size = 4,
             alpha = 0.5,
             position = position_dodge(width = 0.4)) +
  labs(title = "S5 16S Shannon's Diversity",
       y = "Mean Shannon's Diversity",
       x = NULL,
       fill = "Oil Addition", 
       caption = "Error bars represent \U00B1 2 SE") +
  scale_shape_manual(values = c(21,24,23)) +
  theme_bw()  +
  theme(legend.position = "right") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  facet_grid(cols = vars(orig_soil),
             rows = vars(plant_trt),
             labeller = labeller(plant_trt = c("Y" = "Planted", "N" = "No Plant")))


#Means and SE
df <- p$data

prok.alpha <- df %>%
  group_by(tissue, sampling_period, plant_trt, oil_added, orig_soil) %>% 
  summarise_at("value",
               funs(mean,
                     sd, 
                     se=sd(.)/sqrt(n()))) %>% 
  ggplot(aes(x = sampling_period,
           y = mean,
           shape = oil_added,
           fill = tissue))  +
  geom_point(stroke = stroke_size,
             size = point_size,
             alpha = 1,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(min = mean - 2*se,
                    max = mean + 2*se),
                size = stroke_size,
                position = position_dodge(width = 0.4)) +
  labs(title = "Shannon's Diversity - Prokaryote Soil Communities",
       y = "Mean Shannon's Diversity",
       x = NULL,
       fill = "Tissue",
       shape = "Oil Addition",
       caption = "Error bars represent \U00B1 2 SE") +
  scale_shape_manual(values = c(21,24,23)) +
  scale_fill_manual(values = "#619CFF") +
  S5_theme  +
  theme(legend.position = "right") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  facet_grid(cols = vars(orig_soil),
             rows = vars(plant_trt),
             labeller = labeller(plant_trt = c("Y" = "Planted", "N" = "No Plant")))

prok.alpha

# #This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = prok.alpha,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(plot = prok.alpha,
#        filename = "S5_supp_fig4.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
# 
# #For draft and co-author review
# ggsave(plot = prok.alpha,
#        filename = "S5_supp_fig4.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)

## Model

## Can't use plantID or sampling period as random effect because of too few levels
# Validating gamma models: https://stats.stackexchange.com/questions/45401/how-to-validate-diagnose-a-gamma-glm-in-r

#soil
m <- glm(value ~ oil_added*orig_soil*plant_trt, data = p$data %>% #not quite normal
          filter(tissue=="Soil"),
         family = "Gamma")
summary(m)
car::Anova(m, type = "III", test.statistic = "F")
hist(resid(m))
shapiro.test(resid(m))

#No changes when rarefied to 10441 seqs

## Import ITS phyloseq object and clean ----

source("../Scripts for Published Figures and Analyses/import_and_clean_S5_ITS.R")

S5.exp <- prune_samples(sample_data(S5)$special_char=="Experimental Sample", S5)
S5.all <- S5

S5 <- S5.exp

#Break out some elements for easy adjustment of plots----
annotation_size <- 6 / ggplot2::.pt
point_size <- 5 / ggplot2::.pt
stroke_size <- 0.25 / ggplot2::.pt

#color palette
cPAL <- c("#F0E442", "#0072B2")

S5_theme <- theme_bw(base_size = 7) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(7, 4, 2, 2, "mm"),
    axis.title.x = element_text(margin = margin(
      t = 2,
      r = 0,
      b = 0,
      l = 0,
      unit = "mm"
    )),
    axis.title.y = element_text(margin = margin(
      t = 0,
      r = 2,
      b = 0,
      l = 0,
      unit = "mm"
    ))
  )

#Find minimal slope for rarecurve without dropping too many samples

A <- lapply(sort(sample_sums(S5))[c(1:10)], function(x){
  format(median(rareslope(otu_table(S5), sample = x)), scientific = FALSE)
})

B <- lapply(sort(sample_sums(S5))[c(1:10)], function(x){
  format(mean(rareslope(otu_table(S5), sample = x)), scientific = FALSE)
})

C <- lapply(sort(sample_sums(S5))[c(1:10)], function(x){
  sum(sort(sample_sums(S5))[c(1:10)] < x)
})

D <- lapply(sort(sample_sums(S5))[c(1:10)], function(x){
  print(x)
  })

dropdf <- data.frame(median = unlist(A),
                     mean = unlist(B),
                     num_dropped = unlist(C),
                     num_seqs = unlist(D))  

dropdf

#Plot Shannon diversity for each tissue type by time point.
S5.rare <- rarefy_even_depth(physeq = S5,
                              sample.size = 12158)
S5 <- S5.rare
sort(sample_sums(S5))

#Plot Shannon diversity for each tissue type by time point.

p <- plot_richness(S5, measures = "Shannon")

df <- p$data

ITS.alpha <- df %>%
  group_by(tissue, sampling_period, plant_trt, oil_added, orig_soil) %>% 
  summarise_at("value",
               funs(mean,
                     sd, 
                     se=sd(.)/sqrt(n()))) %>% 
  ggplot(aes(x = sampling_period,
           y = mean,
           fill = tissue,
           shape = oil_added))  +
  geom_point(stroke = stroke_size,
             size = point_size,
             alpha = 1,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(min = mean - 2*se,
                    max = mean + 2*se),
                size = stroke_size,
                position = position_dodge(width = 0.4)) +
  labs(title = "Shannon's Diversity - Fungal Communities",
       y = "Mean Shannon's Diversity",
       x = NULL,
       fill = "Tissue",
       shape = "Oil Addition",
       caption = "Error bars represent \U00B1 2 SE") +
  scale_shape_manual(values = c(21,24,23)) +
  S5_theme  +
  theme(legend.position = "right") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  facet_grid(cols = vars(orig_soil),
             rows = vars(plant_trt),
             labeller = labeller(plant_trt = c("Plant" = "Planted", "No Plant" = "No Plant")))

ITS.alpha
# #This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = ITS.alpha,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(plot = ITS.alpha,
#        filename = "S5_supp_fig6.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
# 
# #For draft and co-author review
# ggsave(plot = ITS.alpha,
#        filename = "S5_supp_fig6.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)

## Model
## Can't use plantID or sampling period as random effect because of too few levels
# Validating gamma models: https://stats.stackexchange.com/questions/45401/how-to-validate-diagnose-a-gamma-glm-in-r

m <- glm(value ~ tissue, data = df %>% 
          filter(value!=0),
         family = "Gamma")

summary(m)
car::Anova(m, type = "III", test.statistic = "F")
hist(resid(m))
shapiro.test(resid(m))

#Soil
m <- glm(value ~ sampling_period*oil_added*orig_soil*plant_trt, data = df %>% #not quite normal
          filter(tissue=="Soil", value!=0),
         family = "Gamma")
summary(m)
car::Anova(m, type = "III", test.statistic = "F")
hist(resid(m))
shapiro.test(resid(m))

#Roots
m <- glm(value ~ sampling_period*oil_added*orig_soil, data = df %>% 
          filter(tissue=="Root", value!=0),
         family = "Gamma")
summary(m)
car::Anova(m, type = "III", test.statistic = "F")

hist(resid(m))
shapiro.test(resid(m))

#Leaves
m <- glm(value ~ sampling_period*oil_added*orig_soil, data = df %>% 
          filter(tissue=="Leaf", value!=0),
         family = "Gamma")

summary(m)
car::Anova(m, type = "III", test.statistic = "F")

hist(resid(m))
shapiro.test(resid(m))

#ASV Count and overlap
S5.soil <- subset_samples(sample_data(S5)$tissue=="Soil", physeq = S5) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

S5.root <- subset_samples(sample_data(S5)$tissue=="Root", physeq = S5) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

S5.leaf <- subset_samples(sample_data(S5)$tissue=="Leaf", physeq = S5) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

#Count shared taxa
sum(taxa_names(S5.soil) %in% taxa_names(S5.root))
sum(taxa_names(S5.soil) %in% taxa_names(S5.leaf))
sum(taxa_names(S5.root) %in% taxa_names(S5.leaf))

#all three
sum(taxa_names(S5.soil) %in% taxa_names(S5.leaf) & taxa_names(S5.soil) %in% taxa_names(S5.root))

# What phyla/families are these?

#Soil-Roots
soil_root.taxa <- tax_table(S5.soil)[row.names(tax_table(S5.soil))[taxa_names(S5.soil) %in% taxa_names(S5.root)],]
get_taxa_unique(soil_root.taxa, taxonomic.rank = "Phylum")
get_taxa_unique(soil_root.taxa, taxonomic.rank = "Family")

#Soil-leaf
soil_leaf.taxa <- tax_table(S5.soil)[row.names(tax_table(S5.soil))[taxa_names(S5.soil) %in% taxa_names(S5.leaf)],]
get_taxa_unique(soil_leaf.taxa, taxonomic.rank = "Phylum")
get_taxa_unique(soil_leaf.taxa, taxonomic.rank = "Family")

#Root-leaf
root_leaf.taxa <- tax_table(S5.root)[row.names(tax_table(S5.root))[taxa_names(S5.root) %in% taxa_names(S5.leaf)],]
get_taxa_unique(root_leaf.taxa, taxonomic.rank = "Phylum")
get_taxa_unique(root_leaf.taxa, taxonomic.rank = "Family")

#All three
allthree.taxa <- tax_table(soil_root.taxa)[row.names(tax_table(soil_root.taxa))[taxa_names(soil_root.taxa) %in% taxa_names(soil_leaf.taxa)],]
get_taxa_unique(allthree.taxa, taxonomic.rank = "Phylum")
get_taxa_unique(allthree.taxa, taxonomic.rank = "Family")

#What shows up in the soil when a plant is present?
soil_plant <- subset_samples(S5.soil, sample_data(S5.soil)$plant_trt=="Plant") %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, .)
get_taxa_unique(soil_plant, taxonomic.rank = "Phylum") %>% 
  sort()
get_taxa_unique(soil_plant, taxonomic.rank = "Family") %>% 
  sort()

soil_no_plant <- subset_samples(S5.soil, sample_data(S5.soil)$plant_trt=="No Plant") %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, .)
get_taxa_unique(soil_no_plant, taxonomic.rank = "Phylum") %>% 
  sort()
get_taxa_unique(soil_no_plant, taxonomic.rank = "Family") %>% 
  sort()

#What about the different inocula?
soil_oiled <- subset_samples(S5.soil, sample_data(S5.soil)$orig_soil=="Prev-Oiled Inoc.") %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, .)
get_taxa_unique(soil_oiled, taxonomic.rank = "Phylum") %>% 
  sort()

soil_not_oiled <- subset_samples(S5.soil, sample_data(S5.soil)$orig_soil=="Not Prev-Oiled Inoc") %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, .)
get_taxa_unique(soil_not_oiled, taxonomic.rank = "Phylum") %>% 
  sort()