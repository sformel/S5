#Plant biomass and stem count over time
#Last updated: April 3, 2022
#By Steve Formel

## Description
#Run PERMANOVA and make PCA of 16S based on Aitchison

#load libraries----

library(cowplot)
library(readxl)
library(tidyverse)
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

## Import 16S phyloseq object and clean
#read in data

S5 <-
  readRDS("../Relevant Data/feature_tables_as_phyloseq_objects/S5_16S_pseq.rds")

taxa_names(S5) <- paste0("Seq", seq(ntaxa(S5)))

S5.noCL.noMT <- S5 %>%
  subset_taxa(Family != "Mitochondria" & Order != "Chloroplast")

S5.16S <- S5.noCL.noMT

#fix factors
sample_data(S5) <- sample_data(S5) %>%
  data.frame() %>%
  mutate(across(where(is.character), as.factor)) %>%
  dplyr::mutate(across(c(plantID, table, block), factor)) %>%
  filter(sampling_period == "J18", plant_trt == "Y")


#Set up plotting objects

#Break out some elements for easy adjustment of plots----
annotation_size <- 6 / ggplot2::.pt
point_size <- 5 / ggplot2::.pt
stroke_size <- 0.5 / ggplot2::.pt

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


SF1_aes <- aes(x = PC1,
               y = PC2,
               shape = orig_soil,
               fill = oil_added)

my.scale_aes <- list(geom_point(size = 2),
                     scale_shape_manual(values = c(21, 24)),
                     scale_fill_manual(values = c(cPAL, "black", "gray")))


#16S analysis ----

S5 <- S5.16S
### Plot and PERMANOVA (adonis) of 16S - Aitchison

#By Tissue through time

#takes about 5 min on SF macbook pro
df.all <- otu_table(S5) %>%
  as.data.frame()

df.meta.all <- sample_data(S5) %>%
  data.frame()

#make generic for downstream analysis
df <- df.all
df.meta <- df.meta.all

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.16S <- comm.df

#### Plot PCA - Aitchison
#All tissues

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         shape = orig_soil,
         fill = tissue
       )) +
  my.scale_aes +
  labs(
    title = "PCA on 16S Metagenome - Aitchison",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(cols = vars(sampling_period),
             rows = vars(oil_added)) +
  ggforce::geom_mark_ellipse(aes(linetype = orig_soil))


#NMDS without Soil 3- N17 as outlier
comm.df1 <- comm.df[-132, ]

C <- metaMDS(comm.df1, distance = "euclidean")

D <- as.data.frame(C$points)

ggplot(data = cbind(df.meta[-132, ], D),
       aes(
         x = MDS1,
         y = MDS2,
         shape = orig_soil,
         fill = tissue
       )) +
  my.scale_aes +
  facet_grid(cols = vars(sampling_period),
             rows = vars(oil_added)) +
  ggforce::geom_mark_ellipse(aes(linetype = orig_soil))

## Permanova - Aitchison

set.seed(1)
mod1 <- adonis(
  comm.df ~ tissue,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)
mod1

B <- vegdist(comm.df1, method = "euclidean")

anova(betadisper(B, group = df.meta$tissue[-132]))

#Everything is significant...

## Soil - Did the inoculation work for 16S?

S5.soil <- prune_samples(sample_data(S5)$tissue == "Soil", S5)
S5.soil <-
  prune_samples((
    sample_data(S5.soil)$plantID != "3" |
      sample_data(S5.soil)$sampling_period != "N17"
  ),
  S5.soil)


#takes about 5 min on SF macbook pro
df.soil <- otu_table(S5.soil) %>%
  as.data.frame()

df.meta.soil <- sample_data(S5.soil) %>%
  data.frame()

#make generic for downstream analysis
df <- df.soil
df.meta <- df.meta.soil

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.16S <- comm.df

#### Plot PCA - Aitchison

#Just SOIL

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

#get rid of plant 3 at N17, crazy outlier

ggplot(data = df.PCA %>%
         filter(!(plantID == "3" & sampling_period == "N17")),
       aes(
         x = PC1,
         y = PC2,
         fill = orig_soil,
         shape = oil_added
       )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
                             alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(
    cols = vars(sampling_period),
    labeller = labeller(
      plant_trt = c("Y" = "Planted", "N" = "No Plant"),
      sampling_period = c(
        "N16" = "6 Months",
        "J17" = "1 Year",
        "N17" = "1.5 Years",
        "J18" = "2 Years"
      )
    )
  ) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)


### NMDS

#plantID S3-N17 is a clear outlier

A <- comm.df[rownames(comm.df) != "S3-N17", ]
A <- metaMDS(comm = A, distance = "euclidean")

cbind(df.meta[rownames(df.meta) != "S3-N17", ], A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    fill = orig_soil,
    shape = oil_added
  )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
                             alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment"
  ) +
  facet_grid(
    cols = vars(sampling_period),
    labeller = labeller(
      plant_trt = c("Y" = "Planted", "N" = "No Plant"),
      sampling_period = c(
        "N16" = "6 Months",
        "J17" = "1 Year",
        "N17" = "1.5 Years",
        "J18" = "2 Years"
      )
    )
  ) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL) +
  geom_text(aes(label = plantID))

## NMDS Plot for publication

#The panels I made above were too much and there were some critiques of the PCA, so I've made something that is visually simpler and is an NMDS.  It's hard to see the yellow though, so I'm going to change the colors to something darker.


p.16S.soil <- cbind(df.meta, A$points) %>%
  group_by(plant_trt, oil_added, orig_soil, sampling_period) %>%
  dplyr::summarise(MDS1.mean = mean(MDS1),
                   MDS2.mean = mean(MDS2)) %>%
  ggplot(aes(
    x = MDS1.mean,
    y = MDS2.mean,
    group = interaction(plant_trt, oil_added, orig_soil)
  )) +
  geom_path(
    aes(color = orig_soil,
        linetype = fct_rev(plant_trt)),
    arrow = arrow(type = "open",
                  length = unit(2, "mm")),
    size = 0.5,
    show.legend = TRUE
  ) +
  geom_point(
    data = . %>%
      filter(sampling_period == "N16"),
    aes(shape = oil_added,
        fill = orig_soil),
    size = 3
  ) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(cPAL, "black", "gray")) +
  scale_color_manual(values = c(cPAL, "black", "gray")) +
  S5_theme +
  theme(
    legend.position = "bottom" ,
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    color = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    shape = guide_legend(override.aes = list(
      shape = c(21, 24),
      size = 4,
      linetype = NULL
    ))
  ) +
  labs(
    title = NULL,
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    color = "Inoculum",
    fill = "Inoculum",
    linetype = "Plant Presence",
    x = "MDS1",
    y = "MDS2"
  )

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.16S.soil,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(plot = p.16S.soil,
#        filename = "S5_Figure6.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#
# #For draft and co-author review
# ggsave(plot = p.16S.soil,
#        filename = "S5_Figure6.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)

#below I will combine this plot with the ITS analog

# Permanova - Aitchison

#Just SOIL

set.seed(1)
mod1 <- adonis(
  comm.df ~ plant_trt * orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)
mod1

#check with dbRDA (actually RDA because it's euclidean, but the syntax was more obvious to me in dbRDA than in RDA)

mod1 <-
  dbrda(
    comm.df ~ plant_trt * orig_soil * oil_added + Condition(sampling_period),
    data = df.meta,
    distance = "euclidean"
  )

anova(mod1,
      by = "margin",
      permutations = 9999)


#Everything is significant except for the interaction of all three. But dbRDA doesn't produce a good fit.

## Roots - Did the inoculation work?

S5.roots <- prune_samples(sample_data(S5)$tissue == "Root", S5)

#S5.roots <- prune_samples(sample_data(S5.roots)$season=="Growing", S5.roots)

#takes about 5 min on SF macbook pro
df.roots <- otu_table(S5.roots) %>%
  as.data.frame()

df.meta.roots <- sample_data(S5.roots) %>%
  data.frame()

#make generic for downstream analysis
df <- df.roots
df.meta <- df.meta.roots

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.16S <- comm.df

#### Plot PCA - Aitchison

#Just ROOTS

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         fill = orig_soil,
         shape = oil_added
       )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(color = orig_soil),
                             alpha = 0.2) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    color = "Inoculum",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(cols = vars(season),
             labeller = labeller(season = c(
               "Growing" = "November", "Non-Growing" = "June"
             ))) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

# A <- p +
#   transition_states(sampling_period,
#                     transition_length = 2,
#                     state_length = 1)

## NMDS

A <- comm.df
A <- metaMDS(comm = A, distance = "euclidean")

cbind(df.meta, A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    fill = orig_soil,
    shape = oil_added
  )) +
  my.scale_aes +
  # ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
  #                            alpha = 0.2) +
  # scale_linetype_manual(labels = c("No Plant",
  #                                  "Planted"),
  #                       values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment"
  ) +
  facet_grid(cols = vars(season),
             labeller = labeller(season = c(
               "Growing" = "November", "Non-Growing" = "June"
             ))) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL) +
  scale_fill_manual(values = cPAL)

## Alternative NMDS Plot

#The panels I made above were too much and there were some critiques of the PCA, so I've made something that is visually simpler and is an NMDS.  It's hard to see the yellow though, so I'm going to change the colors to something darker.

p.16S.roots <- cbind(df.meta, A$points) %>%
  group_by(oil_added, orig_soil, sampling_period) %>%
  dplyr::summarise(MDS1.mean = mean(MDS1),
                   MDS2.mean = mean(MDS2)) %>%
  ggplot(aes(
    x = MDS1.mean,
    y = MDS2.mean,
    group = interaction(oil_added, orig_soil)
  )) +
  # geom_link2(aes(alpha = sampling_period),
  #            size = 3) +
  geom_path(
    aes(color = orig_soil),
    arrow = arrow(type = "open",
                  length = unit(2, "mm")),
    size = 0.5,
    show.legend = TRUE
  ) +
  geom_point(
    data = . %>%
      filter(sampling_period == "N16"),
    aes(shape = oil_added,
        fill = orig_soil),
    size = 3
  ) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(cPAL, "black", "gray")) +
  scale_color_manual(values = c(cPAL, "black", "gray")) +
  S5_theme +
  theme(
    legend.position = "bottom" ,
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    color = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    shape = guide_legend(override.aes = list(
      shape = c(21, 24),
      size = 4,
      linetype = NULL
    ))
  ) +
  labs(
    title = NULL,
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    color = "Inoculum",
    fill = "Inoculum",
    x = "MDS1",
    y = "MDS2"
  )

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.16S.roots,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(filename = "S5_Figure7.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#
# #For draft and co-author review
# ggsave(filename = "S5_Figure7.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)


#stats

set.seed(1)
mod1 <- adonis(
  comm.df ~ orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)
mod1

## Permanova - Aitchison

#Just ROOTS

set.seed(1)
mod1 <- adonis(
  comm.df ~ orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)
mod1

#check with dbRDA (actually RDA because it's euclidean, but the syntax was more obvious to me in dbRDA than in RDA)

mod1 <-
  dbrda(
    comm.df ~ orig_soil + oil_added + Condition(sampling_period),
    data = df.meta,
    distance = "euclidean"
  )

anova(mod1,
      by = "margin",
      permutations = 9999)

#Main effects are significant, not the interaction. dbRDA supports this.

### What taxa are responsible?

#quick approach based off a phyloseq tutorial
#https://joey711.github.io/phyloseq/preprocess.html

S5.roots.2 <- S5.roots

transform_sample_counts(S5.roots.2, function(OTU)
  OTU / sum(OTU))
#S5.roots = filter_taxa(S5.roots, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

total = median(sample_sums(S5.roots.2))
standf = function(x, t = total)
  round(t * (x / sum(x)))
gps = transform_sample_counts(S5.roots.2, standf)

#Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
#gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)

#Subset the data to Proteobacteria
gpsfb = subset_taxa(gps, Phylum == "Proteobacteria")

#Now let’s summarize this slice of the data with some graphics.

plot_bar(gpsfb, "sampling_period", "Abundance")

# Other phyla
gpsfb = subset_taxa(gps, Phylum == "Bacteroidetes")
plot_bar(gpsfb, "sampling_period", "Abundance")

# Other phyla
gpsfb = subset_taxa(gps, Phylum == "Cyanobacteria")
plot_bar(gpsfb, "sampling_period", "Abundance")

# Other phyla
gpsfb = subset_taxa(gps, Phylum == "Actinobacteria")
plot_bar(gpsfb, "sampling_period", "Abundance")

# Other phyla
gpsfb = subset_taxa(gps, Phylum == "Firmicutes")
plot_bar(gpsfb, "sampling_period", "Abundance")

# Other phyla
gpsfb = subset_taxa(gps, Phylum == "Epsilonbacteraeota")
plot_bar(gpsfb, "sampling_period", "Abundance")


## Leaves - Did the inoculation work?

#Just LEAVES

S5.leaf <- prune_samples(sample_data(S5)$tissue == "Leaf", S5)

#takes about 5 min on SF macbook pro
df.leaf <- otu_table(S5.leaf) %>%
  as.data.frame()

df.meta.leaf <- sample_data(S5.leaf) %>%
  data.frame()

#make generic for downstream analysis
df <- df.leaf
df.meta <- df.meta.leaf

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.16S <- comm.df

#### Plot PCA - Aitchison

#Just LEAVES

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         fill = orig_soil,
         shape = oil_added
       )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(color = orig_soil),
                             alpha = 0.2) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    color = "Inoculum",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(
               sampling_period = c(
                 "N16" = "6 Months",
                 "J17" = "1 Year",
                 "N17" = "1.5 Years",
                 "J18" = "2 Years"
               )
             )) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## NMDS

A <- comm.df
A <- metaMDS(comm = A, distance = "euclidean")

cbind(df.meta, A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    fill = orig_soil,
    shape = oil_added
  )) +
  my.scale_aes +
  # ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
  #                            alpha = 0.2) +
  # scale_linetype_manual(labels = c("No Plant",
  #                                  "Planted"),
  #                       values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment"
  ) +
  # facet_grid(cols = vars(sampling_period),
  # labeller = labeller(sampling_period = c("N16" = "6 Months", "J17" = "1 Year", "N17" = "1.5 Years", "J18" = "2 Years"))) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## Alternative NMDS Plot

#The panels I made above were too much and there were some critiques of the PCA, so I've made something that is visually simpler and is an NMDS.  It's hard to see the yellow though, so I'm going to change the colors to something darker.

p.16S.leaves <- cbind(df.meta, A$points) %>%
  group_by(oil_added, orig_soil, sampling_period) %>%
  dplyr::summarise(MDS1.mean = mean(MDS1),
                   MDS2.mean = mean(MDS2)) %>%
  ggplot(aes(
    x = MDS1.mean,
    y = MDS2.mean,
    group = interaction(oil_added, orig_soil)
  )) +
  # geom_link2(aes(alpha = sampling_period),
  #            size = 3) +
  geom_path(
    aes(color = orig_soil),
    arrow = arrow(type = "open",
                  length = unit(2, "mm")),
    size = 0.5,
    show.legend = TRUE
  ) +
  geom_point(
    data = . %>%
      filter(sampling_period == "N16"),
    aes(shape = oil_added,
        fill = orig_soil),
    size = 3
  ) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(cPAL, "black", "gray")) +
  scale_color_manual(values = c(cPAL, "black", "gray")) +
  S5_theme +
  theme(
    legend.position = "bottom" ,
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    color = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    shape = guide_legend(override.aes = list(
      shape = c(21, 24),
      size = 4,
      linetype = NULL
    ))
  ) +
  labs(
    title = NULL,
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    color = "Inoculum",
    fill = "Inoculum",
    x = "MDS1",
    y = "MDS1"
  )


#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.16S.leaves,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(filename = "S5_Figure5.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#
# #For draft and co-author review
# ggsave(filename = "S5_Figure5.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)


## Permanova - Aitchison

#Just Leaves

set.seed(1)
mod1 <- adonis(
  comm.df ~ orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)
mod1

A <- vegdist(x = comm.df, method = "euclidean")
B <- betadisper(d = A, group = df.meta$orig_soil)

anova(B)

B <- betadisper(d = A, group = df.meta$oil_added)
anova(B)


set.seed(1)
mod1 <-
  dbrda(
    comm.df ~ orig_soil * oil_added + Condition(sampling_period),
    distance = "euclidean",
    data = df.meta
  )

anova(mod1, by = "margin", permutations = 9999)

#Sampling_period + sampling_period:orig_soil are the only levels that are significant.  If you remove sampling period, nothing else is significant.  dbRDA isn't a good fit.

### What taxa are responsible?

#quick approach based off a phyloseq tutorial
#https://joey711.github.io/phyloseq/preprocess.html

S5.leaf.2 <- S5.leaf

transform_sample_counts(S5.leaf.2, function(OTU)
  OTU / sum(OTU))
#S5.roots = filter_taxa(S5.roots, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

total = median(sample_sums(S5.leaf.2))
standf = function(x, t = total)
  round(t * (x / sum(x)))
gps = transform_sample_counts(S5.leaf.2, standf)

#Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
#gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)

#Subset the data to Proteobacteria
gpsfb = subset_taxa(gps, Phylum == "Proteobacteria")

#Now let’s summarize this slice of the data with some graphics.

plot_bar(gpsfb, "sampling_period", "Abundance")

# Other phyla
gpsfb = subset_taxa(gps, Phylum == "Bacteroidetes")
plot_bar(gpsfb, "sampling_period", "Abundance")

# Other phyla
gpsfb = subset_taxa(gps, Phylum == "Cyanobacteria")
plot_bar(gpsfb, "sampling_period", "Abundance")

# Other phyla
gpsfb = subset_taxa(gps, Phylum == "Actinobacteria")
plot_bar(gpsfb, "sampling_period", "Abundance")

# Other phyla
gpsfb = subset_taxa(gps, Phylum == "Firmicutes")
plot_bar(gpsfb, "sampling_period", "Abundance")

# Other phyla
gpsfb = subset_taxa(gps, Phylum == "Epsilonbacteraeota")
plot_bar(gpsfb, "sampling_period", "Abundance")

# Roots continued
set.seed(1)

S5.ord <- ordinate(S5.roots, distance = "bray")

B <-
  plot_ordination(S5.roots, S5.ord, type = "samples", color = "season")

#What if you get rid of Cyanobactera?

S5.roots_noCB <- subset_taxa(S5.roots, Phylum != "Cyanobacteria")
S5.ord <- ordinate(S5.roots_noCB, distance = "bray")

B <-
  plot_ordination(S5.roots_noCB, S5.ord, type = "samples", color = "season")

#delve into it
S5.roots_noCB <-
  subset_taxa(S5.roots, Phylum == "Proteobacteria" &
                taxa_sums(S5.roots) > 100)
plot_bar(
  S5.roots_noCB,
  x = "oil_added",
  fill = "Phylum",
  facet_grid = ~ season
)

#ITS analsysis----

S5 <- S5.ITS
### Plot of Inocula and Soil together
#takes about 5 min on SF macbook pro

df.all <- otu_table(S5) %>%
  as.data.frame()

df.meta.all <- sample_data(S5) %>%
  data.frame()

#make generic for downstream analysis
df <- df.all
df.meta <- df.meta.all

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.ITS <- comm.df

#### Plot PCA - Aitchison
#All tissues + Inocula

#make variable for inocula

df.meta$orig_soil[str_detect(df.meta$sampleID_stem, pattern = "O$")] <-
  "Prev-Oiled Inoc."
df.meta$orig_soil[str_detect(df.meta$sampleID_stem, pattern = "U$")] <-
  "Not Prev-Oiled Inoc"


#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         shape = orig_soil,
         fill = tissue
       )) +
  geom_point(size = 2) +
  scale_shape_manual(values = c(21, 24)) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "PCA on ITS Metagenome - Aitchison",
    caption = "Ellipses are visual guides and have no statistical meaning",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)"),
    shape = "Inoculum"
  ) +
  ggforce::geom_mark_ellipse(aes(color = tissue),
                             alpha = 0.1)


### Plot and PERMANOVA (adonis) of ITS - Aitchison

#By Tissue through time

#takes about 5 min on SF macbook pro
df.all <- otu_table(S5) %>%
  as.data.frame()

df.meta.all <- sample_data(S5) %>%
  data.frame()

#make generic for downstream analysis
df <- df.all
df.meta <- df.meta.all

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.ITS <- comm.df

#### Plot PCA - Aitchison
#All tissues


#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         shape = orig_soil,
         fill = tissue
       )) +
  geom_point(size = 2) +
  scale_shape_manual(values = c(21, 24)) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "PCA on ITS Metagenome - Aitchison",
    caption = "Ellipses are visual guides and have no statistical meaning",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)"),
    shape = "Inoculum"
  ) +
  facet_grid(cols = vars(sampling_period),
             rows = vars(oil_added)) +
  ggforce::geom_mark_ellipse(aes(linetype = orig_soil))


## Permanova - Aitchison

#Control for sampling period.  This means permutations will only happen within a sampling period.  It isn't exactly repeated-measures, but it's in th right d

set.seed(1)
mod1 <- adonis(
  comm.df ~ tissue,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  df.meta$sampling_period
)
mod1

## Soil - Did the inoculation work?

S5.soil <- prune_samples(sample_data(S5)$tissue == "Soil", S5)

#takes about 5 min on SF macbook pro
df.soil <- otu_table(S5.soil) %>%
  as.data.frame()

df.meta.soil <- sample_data(S5.soil) %>%
  data.frame()

#make generic for downstream analysis
df <- df.soil
df.meta <- df.meta.soil

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.ITS <- comm.df

#### Plot PCA - Aitchison

#Just SOIL

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         fill = orig_soil,
         shape = oil_added
       )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
                             alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(
    cols = vars(sampling_period),
    labeller = labeller(
      plant_trt = c("Y" = "Planted", "N" = "No Plant"),
      sampling_period = c(
        "N16" = "6 Months",
        "J17" = "1 Year",
        "N17" = "1.5 Years",
        "J18" = "2 Years"
      )
    )
  ) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

### NMDS

A <- comm.df
A <- metaMDS(comm = A, distance = "euclidean")

cbind(df.meta, A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    fill = orig_soil,
    shape = oil_added
  )) +
  my.scale_aes +
  #ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
  #                           alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment"
  ) +
  #   facet_grid(cols = vars(sampling_period),
  #              labeller = labeller(plant_trt = c("Y" = "Planted", "N" = "No Plant"),
  #                                  sampling_period = c("N16" = "6 Months", "J17" = "1 Year", "N17" = "1.5 Years", "J18" = "2 Years"))) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## Alternative NMDS Plot

#The panels I made above were too much and there were some critiques of the PCA, so I've made something that is visually simpler and is an NMDS.  It's hard to see the yellow though, so I'm going to change the colors to something darker.

p.ITS.soil <- cbind(df.meta, A$points) %>%
  group_by(plant_trt, oil_added, orig_soil, sampling_period) %>%
  dplyr::summarise(MDS1.mean = mean(MDS1),
                   MDS2.mean = mean(MDS2)) %>%
  ggplot(aes(
    x = MDS1.mean,
    y = MDS2.mean,
    group = interaction(plant_trt, oil_added, orig_soil)
  )) +
  # geom_link2(aes(alpha = sampling_period),
  #            size = 3) +
  geom_path(
    aes(color = orig_soil,
        linetype = fct_rev(plant_trt)),
    arrow = arrow(type = "open",
                  length = unit(2, "mm")),
    size = 0.5,
    show.legend = TRUE
  ) +
  geom_point(
    data = . %>%
      filter(sampling_period == "N16"),
    aes(shape = oil_added,
        fill = orig_soil),
    size = 3
  ) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(cPAL, "black", "gray")) +
  scale_color_manual(values = c(cPAL, "black", "gray")) +
  S5_theme +
  theme(
    legend.position = "bottom" ,
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    color = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    shape = guide_legend(override.aes = list(
      shape = c(21, 24),
      size = 4,
      linetype = NULL
    ))
  ) +
  labs(
    title = NULL,
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    color = "Inoculum",
    fill = "Inoculum",
    linetype = "Plant Presence",
    x = "MDS1",
    y = "MDS1"
  )

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.ITS.soil,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(filename = "S5_Figure5.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#
# #For draft and co-author review
# ggsave(filename = "S5_Figure5.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#


## Permanova - Aitchison

#Just SOIL

set.seed(1)

mod1 <- adonis(
  comm.df ~ plant_trt * orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)

mod1

#check with dbRDA (actually RDA because it's euclidean, but the syntax was more obvious to me in dbRDA than in RDA)

mod1 <-
  dbrda(
    comm.df ~ plant_trt + orig_soil + oil_added + Condition(sampling_period),
    data = df.meta,
    distance = "euclidean"
  )

anova(mod1,
      by = "margin",
      permutations = 9999)


#sampling period, plant_trt, orig_soil and the interaction of plant_trt:orig_soil are significant.  But it's still flimsy.  If I take sampling period out, it's still the same, results otherwise.  So I mine as well keep it in.

## dbRDA with oil components + CN

mod1 <-
  dbrda(
    comm.df ~ pH + conductivity + Hopanes + Total.Naph. + Total.Dibenz. + Total.Phenan. + Total.Chrys. + Condition(sampling_period + table + block),
    data = df.meta,
    distance = "euclidean"
  )

mod2 <- dbrda(comm.df ~ 1,
              data = df.meta,
              distance = "euclidean")

ordistep(mod2, scope = formula(mod1))

anova(mod1,
      by = "margin",
      permutations = 9999)

#Significant for Hopanes, pH, table, block.  After controlling for table and block, Hopanes are significant.  Oil addition wasn't significant in the permanova, so that's not very reassuring tha these results mean much.

## Roots - Did the inoculation work?

S5.roots <- prune_samples(sample_data(S5)$tissue == "Root", S5)

#takes about 5 min on SF macbook pro
df.roots <- otu_table(S5.roots) %>%
  as.data.frame()

df.meta.roots <- sample_data(S5.roots) %>%
  data.frame()

#make generic for downstream analysis
df <- df.roots
df.meta <- df.meta.roots

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.ITS <- comm.df

#### Plot PCA - Aitchison

#Just ROOTS

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

#remove plant 60, J18 as an outlier

ggplot(data = df.PCA %>%
         filter(!(plantID == "60" & sampling_period == "J18")),
       aes(
         x = PC1,
         y = PC2,
         fill = orig_soil,
         shape = oil_added
       )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(color = orig_soil),
                             alpha = 0.2) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    color = "Inoculum",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(
               sampling_period = c(
                 "N16" = "6 Months",
                 "J17" = "1 Year",
                 "N17" = "1.5 Years",
                 "J18" = "2 Years"
               )
             )) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## NMDS

A <- comm.df[rownames(comm.df) != "S5_R60_J18", ]
A <- metaMDS(comm = A, distance = "euclidean")

cbind(df.meta[rownames(df.meta) != "S5_R60_J18", ], A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    fill = orig_soil,
    shape = oil_added
  )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
                             alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment"
  ) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(
               sampling_period = c(
                 "N16" = "6 Months",
                 "J17" = "1 Year",
                 "N17" = "1.5 Years",
                 "J18" = "2 Years"
               )
             )) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## Alternative NMDS Plot

#The panels I made above were too much and there were some critiques of the PCA, so I've made something that is visually simpler and is an NMDS.  It's hard to see the yellow though, so I'm going to change the colors to something darker.


p.ITS.roots <-
  cbind(df.meta[rownames(df.meta) != "S5_R60_J18", ], A$points) %>%
  group_by(oil_added, orig_soil, sampling_period) %>%
  dplyr::summarise(MDS1.mean = mean(MDS1),
                   MDS2.mean = mean(MDS2)) %>%
  ggplot(aes(
    x = MDS1.mean,
    y = MDS2.mean,
    group = interaction(oil_added, orig_soil)
  )) +
  # geom_link2(aes(alpha = sampling_period),
  #            size = 3) +
  geom_path(
    aes(color = orig_soil),
    arrow = arrow(type = "open",
                  length = unit(2, "mm")),
    size = 0.5,
    show.legend = TRUE
  ) +
  geom_point(
    data = . %>%
      filter(sampling_period == "N16"),
    aes(shape = oil_added,
        fill = orig_soil),
    size = 3
  ) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(cPAL, "black", "gray")) +
  scale_color_manual(values = c(cPAL, "black", "gray")) +
  S5_theme +
  theme(
    legend.position = "bottom" ,
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    color = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    shape = guide_legend(override.aes = list(
      shape = c(21, 24),
      size = 4,
      linetype = NULL
    ))
  ) +
  labs(
    title = NULL,
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    color = "Inoculum",
    fill = "Inoculum",
    x = "MDS1",
    y = "MDS1"
  )


#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.ITS.roots,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(filename = "S5_Figure5.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#
# #For draft and co-author review
# ggsave(filename = "S5_Figure5.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#
## Permanova - Aitchison

#Just ROOTS

set.seed(1)
mod1 <- adonis(
  comm.df ~ orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)
mod1

#check with dbRDA (actually RDA because it's euclidean, but the syntax was more obvious to me in dbRDA than in RDA)

mod1 <-
  dbrda(
    comm.df ~ orig_soil + oil_added + Condition(sampling_period),
    data = df.meta,
    distance = "euclidean"
  )

anova(mod1,
      by = "margin",
      permutations = 9999)

#Orig_soil, oil_added, orig_soil:oil_added. dbRDA backs up the main effects idea.

## dbRDA with oil components + CN

comm.df.sub <- subset(comm.df, rownames(comm.df) != "S5_R29_N17")
df.meta.sub <- df.meta %>%
  filter(plantID != 29 | sampling_period != "N17")

mod1 <-
  dbrda(
    comm.df.sub ~ percent_N + percent_C + C.N + table + block++stem_diam + num_nodes + num_stems_live + pH + conductivity + Hopanes + Total.Naph. + Total.Dibenz. + Total.Phenan. + Total.Chrys. + Condition(sampling_period),
    data = df.meta.sub,
    distance = "euclidean"
  )

mod2 <- dbrda(comm.df.sub ~ 1,
              data = df.meta.sub,
              distance = "euclidean")

ordistep(mod2, scope = formula(mod1))

anova(mod1,
      by = "margin",
      permutations = 9999)

#Block is the only significant factor.

## Leaves - Did the inoculation work?

#Just LEAVES

S5.leaf <- prune_samples(sample_data(S5)$tissue == "Leaf", S5)

#takes about 5 min on SF macbook pro
df.leaf <- otu_table(S5.leaf) %>%
  as.data.frame()

df.meta.leaf <- sample_data(S5.leaf) %>%
  data.frame()

#make generic for downstream analysis
df <- df.leaf
df.meta <- df.meta.leaf

clr.df <- clr(df)
pca = prcomp(clr.df, center = FALSE, scale. = FALSE)

PCA.loadings = pca$rotation
PCA.scores = pca$x

#For below analyses
comm.df <- clr.df %>%
  as.data.frame()

#for PERMANOVA with all chemistry
comm.df.ITS <- comm.df


#### Plot PCA - Aitchison

#Just LEAVES

#get summary data
pca.sum <- as.data.frame(summary(pca)$importance)
df.PCA <- cbind(PCA.scores, df.meta)

ggplot(data = df.PCA,
       aes(
         x = PC1,
         y = PC2,
         fill = orig_soil,
         shape = oil_added
       )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(color = orig_soil),
                             alpha = 0.2) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    color = "Inoculum",
    x = paste("PC1 (", round(pca.sum$PC1[2] * 100, digits = 1), "%)"),
    y = paste("PC2 (", round(pca.sum$PC2[2] * 100, digits = 1), "%)")
  ) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(
               sampling_period = c(
                 "N16" = "6 Months",
                 "J17" = "1 Year",
                 "N17" = "1.5 Years",
                 "J18" = "2 Years"
               )
             )) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## NMDS

A <- comm.df
A <- metaMDS(comm = A, distance = "euclidean")

cbind(df.meta, A$points) %>%
  ggplot(aes(
    x = MDS1,
    y = MDS2,
    fill = orig_soil,
    shape = oil_added
  )) +
  my.scale_aes +
  ggforce::geom_mark_ellipse(aes(linetype = plant_trt),
                             alpha = 0.2) +
  scale_linetype_manual(labels = c("No Plant",
                                   "Planted"),
                        values = c("dashed", "solid")) +
  labs(
    title = NULL,
    caption = "Ellipses are visual guides and have no statistical meaning",
    shape = "Oil Addition",
    fill = "Inoculum",
    linetype = "Plant Treatment"
  ) +
  facet_grid(cols = vars(sampling_period),
             labeller = labeller(
               sampling_period = c(
                 "N16" = "6 Months",
                 "J17" = "1 Year",
                 "N17" = "1.5 Years",
                 "J18" = "2 Years"
               )
             )) +
  theme(legend.position = "bottom" , legend.direction = "vertical") +
  scale_color_manual(values = cPAL)

## Alternative NMDS Plot

#The panels I made above were too much and there were some critiques of the PCA, so I've made something that is visually simpler and is an NMDS.  It's hard to see the yellow though, so I'm going to change the colors to something darker.

p.ITS.leaves <- cbind(df.meta, A$points) %>%
  group_by(oil_added, orig_soil, sampling_period) %>%
  dplyr::summarise(MDS1.mean = mean(MDS1),
                   MDS2.mean = mean(MDS2)) %>%
  ggplot(aes(
    x = MDS1.mean,
    y = MDS2.mean,
    group = interaction(oil_added, orig_soil)
  )) +
  geom_path(
    aes(color = orig_soil),
    arrow = arrow(type = "open",
                  length = unit(2, "mm")),
    size = 0.5,
    show.legend = TRUE
  ) +
  geom_point(
    data = . %>%
      filter(sampling_period == "N16"),
    aes(shape = oil_added,
        fill = orig_soil),
    size = 3
  ) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(cPAL, "black", "gray")) +
  scale_color_manual(values = c(cPAL, "black", "gray")) +
  S5_theme +
  theme(
    legend.position = "bottom" ,
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    color = guide_legend(override.aes = list(
      shape = 22,
      size = 4,
      linetype = NULL
    )),
    shape = guide_legend(override.aes = list(
      shape = c(21, 24),
      size = 4,
      linetype = NULL
    ))
  ) +
  labs(
    title = NULL,
    caption = paste0("Stress = ", round(A$stress, 2)),
    shape = "Oil Addition",
    color = "Inoculum",
    fill = "Inoculum",
    x = "MDS1",
    y = "MDS1"
  )

#This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p.ITS.leaves,
#                     width = 85,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

# ggsave(filename = "S5_Figure5.pdf",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#
# #For draft and co-author review
# ggsave(filename = "S5_Figure5.png",
#        path = "figures/",
#        width = 85,
#        height = 85,
#        units = "mm",
#        dpi = 300,
#        scale = 1)
#

## Permanova - Aitchison

#Just Leaves
set.seed(1)
mod1 <- adonis(
  comm.df ~ orig_soil * oil_added,
  data = df.meta,
  permutations = 9999,
  method = "euclidean",
  strata = df.meta$sampling_period
)
mod1

#check with dbRDA (actually RDA because it's euclidean, but the syntax was more obvious to me in dbRDA than in RDA)

mod1 <-
  dbrda(
    comm.df ~ orig_soil + oil_added + Condition(sampling_period),
    data = df.meta,
    distance = "euclidean"
  )

anova(mod1,
      by = "margin",
      permutations = 9999)

#orig_soil, and orig_soil:oil_added are the only levels that are significant. dbRDA supports main effect.

## dbRDA with variables but not treatments

comm.df.sub <- subset(comm.df, rownames(comm.df) != "S5_L08_N16")
df.meta.sub <- df.meta %>%
  filter(plantID != 8 | sampling_period != "N16")

mod1 <-
  dbrda(
    comm.df.sub ~ percent_N + percent_C + C.N + table + block++stem_ht + stem_diam + num_nodes + num_stems_live + pH + conductivity + Condition(sampling_period),
    data = df.meta.sub,
    distance = "euclidean"
  )

mod2 <- dbrda(comm.df.sub ~ 1,
              data = df.meta.sub,
              distance = "euclidean")

ordistep(mod2, scope = formula(mod1))

anova(mod1,
      by = "margin",
      permutations = 9999)

#The only significant thing here is table and block.

#Save plots -----

#two columns width = 180mm single column = 85mm

#soil

soil.legend <- get_legend(
  p.16S.soil +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(
        t = 1,
        r = 1,
        b = 1,
        l = 1
      ),
      legend.title.align = 0.5,
      legend.background = element_blank(),
      legend.box.background = element_rect(colour = "lightgray",
                                           size = stroke_size),
      legend.key.width = unit(5, "mm")
    )
)

top <-
  plot_grid(
    p.16S.soil + S5_theme + theme(legend.position = "none"),
    p.ITS.soil + S5_theme + theme(legend.position = "none"),
    labels = c("Prokaryotes", "Fungi"),
    hjust = c(-0.38, -0.71)
  )

p <- plot_grid(top,
               soil.legend,
               ncol = 1,
               rel_heights = c(1, 0.3))

# #This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p,
#                     width = 180,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

ggsave(
  filename = "S5_figure4.pdf",
  path = "figures/",
  width = 180,
  height = 85,
  units = "mm",
  dpi = 300,
  scale = 1
)

# #For draft and co-author review
ggsave(
  filename = "S5_figure4.png",
  path = "figures/",
  width = 180,
  height = 85,
  units = "mm",
  dpi = 300,
  scale = 1
)

pdf_convert(
  pdf = "figures/S5_figure4.pdf",
  format = "tiff",
  filenames = "figures/S5_figure4.tif",
  dpi = 300
)

#roots

roots.legend <- get_legend(
  p.16S.roots +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(
        t = 1,
        r = 1,
        b = 1,
        l = 1
      ),
      legend.title.align = 0.5,
      legend.background = element_blank(),
      legend.box.background = element_rect(colour = "lightgray",
                                           size = stroke_size),
      legend.key.width = unit(5, "mm")
    )
)

top <-
  plot_grid(
    p.16S.roots + S5_theme + theme(legend.position = "none"),
    p.ITS.roots + S5_theme + theme(legend.position = "none"),
    labels = c("Prokaryotes", "Fungi"),
    hjust = c(-0.38, -0.71)
  )

p <- plot_grid(top,
               roots.legend,
               ncol = 1,
               rel_heights = c(1, 0.3))

# #This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p,
#                     width = 180,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

ggsave(
  filename = "S5_figure5.pdf",
  path = "figures/",
  width = 180,
  height = 85,
  units = "mm",
  dpi = 300,
  scale = 1
)

# #For draft and co-author review
ggsave(
  filename = "S5_figure5.png",
  path = "figures/",
  width = 180,
  height = 85,
  units = "mm",
  dpi = 300,
  scale = 1
)

pdf_convert(
  pdf = "figures/S5_figure5.pdf",
  format = "tiff",
  filenames = "figures/S5_figure5.tif",
  dpi = 300
)

#leaves

leaves.legend <- get_legend(
  p.16S.leaves +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(
        t = 1,
        r = 1,
        b = 1,
        l = 1
      ),
      legend.title.align = 0.5,
      legend.background = element_blank(),
      legend.box.background = element_rect(colour = "lightgray",
                                           size = stroke_size),
      legend.key.width = unit(5, "mm")
    )
)

top <-
  plot_grid(
    p.16S.leaves + S5_theme + theme(legend.position = "none"),
    p.ITS.leaves + S5_theme + theme(legend.position = "none"),
    labels = c("Prokaryotes", "Fungi"),
    hjust = c(-0.38, -0.71)
  )

p <- plot_grid(top,
               leaves.legend,
               ncol = 1,
               rel_heights = c(1, 0.3))

# #This was extremely helpful for tweaking the plot
# nflplotR::ggpreview(plot = p,
#                     width = 180,
#                     height = 85,
#                     units = "mm",
#                     dpi = 300,
#                     scale = 1,
#                     device = "pdf")

ggsave(
  filename = "S5_figure6.pdf",
  path = "figures/",
  width = 180,
  height = 85,
  units = "mm",
  dpi = 300,
  scale = 1
)

# #For draft and co-author review
ggsave(
  filename = "S5_figure6.png",
  path = "figures/",
  width = 180,
  height = 85,
  units = "mm",
  dpi = 300,
  scale = 1
)

#convert pdf to TIFF for publisher
library(pdftools)

pdf_convert(
  pdf = "figures/S5_figure6.pdf",
  format = "tiff",
  filenames = "figures/S5_figure6.tif",
  dpi = 300
)