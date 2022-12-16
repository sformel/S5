#Graphic representation of community composition by phyla
#Last updated: 2022-12-16

## Description

#load libraries----

library(cowplot)
library(readxl)
library(dplyr)
library(compositions)
library(vegan)
library(plotrix)
library(phyloseq)
library(data.table)
library(reshape2)
library(car)
library(compositions)
library(ggcorrplot)
library(ggforce)

#Read in data and clean-----

#Import ITS data and clean

#read in data
source("./import_and_clean_S5_ITS.R")

S5.exp <-
  prune_samples(sample_data(S5)$special_char == "Experimental Sample", S5)
S5.all <- S5

S5.ITS <- S5.exp

S5.rare <- rarefy_even_depth(physeq = S5.ITS,
                              sample.size = 3914)
S5.ITS <- S5.rare

## Import 16S phyloseq object and clean
#read in data

S5 <- readRDS("../Relevant Data/feature_tables_as_phyloseq_objects/S5_16S_pseq.rds")

taxa_names(S5) <- paste0("Seq", seq(ntaxa(S5)))

S5.noCL.noMT <- S5 %>%
  subset_taxa(Family != "Mitochondria" & Order != "Chloroplast")

S5.16S <- S5.noCL.noMT

S5 <- S5.16S

S5.rare <- rarefy_even_depth(physeq = S5,
                              sample.size = 1039)
S5 <- S5.rare

sort(sample_sums(S5))

#fix factors
sample_data(S5) <- sample_data(S5) %>%
  data.frame() %>%
  mutate(across(where(is.character), as.factor)) %>%
  dplyr::mutate(across(c(plantID, table, block), factor))

#Set up plotting objects

#Break out some elements for easy adjustment of plots----
annotation_size <- 5 / ggplot2::.pt
point_size <- 5 / ggplot2::.pt
stroke_size <- 0.5 / ggplot2::.pt

#color palette
cPAL <- c("#F0E442", "#0072B2")

S5_theme <- theme_bw(base_size = 7) +
  theme(
    legend.spacing.y = unit(-5, units = "pt"),
    legend.title = element_text(margin = margin(0,0,3,0, unit = "pt")),
    legend.key = element_rect(fill = "transparent"),
    panel.grid.major = element_blank(),
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


SF1_aes <- aes(x = PC1,
               y = PC2,
               shape = orig_soil,
               fill = oil_added)

my.scale_aes <- list(geom_point(size = 2),
                     scale_shape_manual(values = c(21, 24)),
                     scale_fill_manual(values = c(cPAL, "black", "gray")))

### Fungal Soil composition

S5 <- subset_samples(S5.ITS, sample_data(S5.ITS)$tissue=="Soil") %>% 
  prune_taxa(taxa_sums(.) > 0, .) 

#Sum by phylum
S5 <- tax_glom(S5, taxrank="Phylum")

otu <- otu_table(S5) %>%
  as.data.frame() %>% 
  t() %>% 
  as.data.frame()

meta <- sample_data(S5) %>%
  data.frame()

clr.df <- clr(otu)
comm.df <- clr.df %>%
  as.data.frame()

otu_table(S5) <- otu_table(comm.df, taxa_are_rows = TRUE)

sample_data(S5)$plantID <- sample_data(S5)$plantID %>% 
  as.factor()

p <- plot_bar(S5, x = "plantID", 
         y = "Abundance", 
         fill = "Phylum")

#strip prefix from phylum

p$data$Phylum <- sub(pattern = "p__", replacement = "", x = p$data$Phylum)

#make interaction for x-axis grouping
p$data$meso_by_plant <- fct_cross(p$data$plantID, p$data$plant_trt)
levels(p$data$meso_by_plant) <- sub(x = levels(p$data$meso_by_plant), pattern = ":.*",replacement = "")

#Reverse plant Trt Levels so shapes are consistent for other plots
p$data$plant_trt <- fct_rev(p$data$plant_trt)

#sampling_period levels
levels(p$data$sampling_period) <- c("6 months", "1 year", "1.5 years", "2 years")

#Treatment levels
levels(p$data$Treatment) <- gsub(x = levels(p$data$Treatment),pattern = " ; ", replacement = "\n")
  
q <- p$data %>% 
  ggplot(aes(x = meso_by_plant,
             y = Abundance,
             color = Phylum,
             shape = plant_trt)) +
  facet_grid(rows = vars(sampling_period),
             cols = vars(Treatment),
             scales = "free_x") + 
  geom_jitter(height = NULL,
              width = 0.15,
              size = point_size, 
             alpha = 0.7) +
    geom_hline(yintercept = 0, 
             linetype = 3, 
             size = stroke_size) +
  S5_theme +
  guides(shape = guide_legend(byrow = TRUE),
         color = guide_legend(byrow = TRUE)) +
  labs(shape = "Plant Presence",
       x = "Meoscosm ID",
       y = "centered log ratio of taxon abundance") 

q

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = q,
#                     width = 170,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(plot = q,
#        "S5_fungal_soil_composition_supp.png", 
#        path = "figures/",
#        width = 170,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)

### Fungal Root composition

S5 <- subset_samples(S5.ITS, sample_data(S5.ITS)$tissue=="Root") %>% 
  prune_taxa(taxa_sums(.) > 0, .) 

#Sum by phylum
S5 <- tax_glom(S5, taxrank="Phylum")

otu <- otu_table(S5) %>%
  as.data.frame() %>% 
  t() %>% 
  as.data.frame()

meta <- sample_data(S5) %>%
  data.frame()

clr.df <- clr(otu)
comm.df <- clr.df %>%
  as.data.frame()

otu_table(S5) <- otu_table(comm.df, taxa_are_rows = TRUE)

sample_data(S5)$plantID <- sample_data(S5)$plantID %>% 
  as.factor()

p <- plot_bar(S5, x = "plantID", 
         y = "Abundance", 
         fill = "Phylum")

#strip prefix from phylum
p$data$Phylum <- sub(pattern = "p__", replacement = "", x = p$data$Phylum)

#make interaction for x-axis grouping
p$data$meso_by_plant <- fct_cross(p$data$plantID, p$data$plant_trt)
levels(p$data$meso_by_plant) <- sub(x = levels(p$data$meso_by_plant), pattern = ":.*",replacement = "")

#Reverse plant Trt Levels so shapes are consistent for other plots
p$data$plant_trt <- fct_rev(p$data$plant_trt)

#sampling_period levels
levels(p$data$sampling_period) <- c("6 months", "1 year", "1.5 years", "2 years")

#Treatment levels
levels(p$data$Treatment) <- gsub(x = levels(p$data$Treatment),pattern = " ; ", replacement = "\n")

q <- p$data %>% 
  ggplot(aes(x = meso_by_plant,
             y = Abundance,
             color = Phylum)) +
  facet_grid(rows = vars(sampling_period),
             cols = vars(Treatment),
             scales = "free_x") + 
  geom_jitter(height = NULL,
              width = 0.15,
              size = point_size, 
             alpha = 0.7) +
    geom_hline(yintercept = 0, 
             linetype = 3, 
             size = stroke_size) +
  S5_theme +
  guides(shape = guide_legend(byrow = TRUE),
         color = guide_legend(byrow = TRUE)) +
  labs(shape = "Plant Presence",
       x = "Meoscosm ID",
       y = "centered log ratio of taxon abundance") 

q

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = q,
#                     width = 170,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")
# 
# ggsave(plot = q,
#        "S5_fungal_root_composition_supp.png", 
#        path = "figures/",
#        width = 170,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)


### Fungal Leave composition
S5 <- subset_samples(S5.ITS, sample_data(S5.ITS)$tissue=="Leaf") %>% 
  prune_taxa(taxa_sums(.) > 0, .) 

#Sum by phylum
S5 <- tax_glom(S5, taxrank="Phylum")

otu <- otu_table(S5) %>%
  as.data.frame() %>% 
  t() %>% 
  as.data.frame()

meta <- sample_data(S5) %>%
  data.frame()

clr.df <- clr(otu)
comm.df <- clr.df %>%
  as.data.frame()

otu_table(S5) <- otu_table(comm.df, taxa_are_rows = TRUE)

sample_data(S5)$plantID <- sample_data(S5)$plantID %>% 
  as.factor()

p <- plot_bar(S5, x = "plantID", 
         y = "Abundance", 
         fill = "Phylum")

#strip prefix from phylum
p$data$Phylum <- sub(pattern = "p__", replacement = "", x = p$data$Phylum)

#make interaction for x-axis grouping
p$data$meso_by_plant <- fct_cross(p$data$plantID, p$data$plant_trt)
levels(p$data$meso_by_plant) <- sub(x = levels(p$data$meso_by_plant), pattern = ":.*",replacement = "")

#Reverse plant Trt Levels so shapes are consistent for other plots
p$data$plant_trt <- fct_rev(p$data$plant_trt)

#sampling_period levels
levels(p$data$sampling_period) <- c("6 months", "1 year", "1.5 years", "2 years")

#Treatment levels
levels(p$data$Treatment) <- gsub(x = levels(p$data$Treatment),pattern = " ; ", replacement = "\n")

q <- p$data %>% 
  ggplot(aes(x = meso_by_plant,
             y = Abundance,
             color = Phylum)) +
  facet_grid(rows = vars(sampling_period),
             cols = vars(Treatment),
             scales = "free_x") + 
  geom_jitter(height = NULL,
              width = 0.15,
              size = point_size, 
             alpha = 0.7) +
    geom_hline(yintercept = 0, 
             linetype = 3, 
             size = stroke_size) +
  S5_theme +
  guides(shape = guide_legend(byrow = TRUE),
         color = guide_legend(byrow = TRUE)) +
  labs(shape = "Plant Presence",
       x = "Meoscosm ID",
       y = "centered log ratio of taxon abundance") 

q

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = q,
#                     width = 170,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")
# 
# ggsave(plot = q,
#        "S5_fungal_leaf_composition_supp.png", 
#        path = "figures/",
#        width = 170,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)


### Prokaryote Soil composition

S5 <- subset_samples(S5.16S, sample_data(S5.16S)$tissue=="Soil") %>% 
  prune_taxa(taxa_sums(.) > 0, .) 

#Sum by phylum
S5 <- tax_glom(S5, taxrank="Phylum")

otu <- otu_table(S5) %>%
  as.data.frame() %>% 
  t() %>% 
  as.data.frame()

meta <- sample_data(S5) %>%
  data.frame()

clr.df <- clr(otu)
comm.df <- clr.df %>%
  as.data.frame()

otu_table(S5) <- otu_table(comm.df, taxa_are_rows = TRUE)

sample_data(S5)$plantID <- sample_data(S5)$plantID %>% 
  as.factor()

p <- plot_bar(S5, x = "plantID", 
         y = "Abundance", 
         fill = "Phylum")

#make interaction for x-axis grouping
p$data$meso_by_plant <- fct_cross(p$data$plantID, p$data$plant_trt)
levels(p$data$meso_by_plant) <- sub(x = levels(p$data$meso_by_plant), pattern = ":.*",replacement = "")

#Reverse plant Trt Levels so shapes are consistent for other plots
p$data$plant_trt <- fct_rev(p$data$plant_trt)

#sampling_period levels
levels(p$data$sampling_period) <- c("6 months", "1 year", "1.5 years", "2 years")

#Treatment levels
levels(p$data$Treatment) <- gsub(x = levels(p$data$Treatment),pattern = " ; ", replacement = "\n")

q <- p$data %>% 
  ggplot(aes(x = meso_by_plant,
             y = Abundance,
             shape = plant_trt,
             color = Phylum)) +
  facet_grid(rows = vars(sampling_period),
             cols = vars(Treatment),
             scales = "free_x") + 
  geom_jitter(height = NULL,
              width = 0.15,
              size = point_size, 
             alpha = 0.7) +
    geom_hline(yintercept = 0, 
             linetype = 3, 
             size = stroke_size) +
  S5_theme +
  guides(shape = guide_legend(byrow = TRUE),
         color = guide_legend(byrow = TRUE)) +
  labs(shape = "Plant Presence",
       x = "Meoscosm ID",
       y = "centered log ratio of taxon abundance") +
  theme(legend.position = "none")

q

# #This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = q,
#                     width = 170,
#                     height = 100,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")
# 
# ggsave(plot = q,
#        "S5_prok_soil_composition_supp_no_legend.png", 
#        path = "figures/",
#        width = 170,
#        height = 100,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
# 
# #plot legend only
# q.leg <- get_legend(q + theme(legend.position = "right"))
# 
# ggsave(plot = q.leg,
#        "S5_prok_soil_composition_supp_legend.png", 
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
