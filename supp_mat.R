#### Supplementary materials ####

# Code last update: 2026-06-29
# R version 4.4.2 (2024-10-31 ucrt)
# Details of the R session at the end of the script

# Author:
# Elena Ianeva

# Clean R memory
rm(list =ls()) 

#### Load Packages #### 

# Load required packages for data manipulation, statistics, and visualization
library(mixOmics)
library(tidyverse) # Edit the data
library(dplyr) # To modify the data
library(ggplot2) # To make the plots
library(plotly)
library(viridis)
library(vegan)
library(pairwiseAdonis)
library(gt)
library(ggordiplots)
library(patchwork)




##### 1. Pilot study #####

# PLS-DA (or sPLS-DA) of samples from the pilot study (the 3 exposure times 
# + control group)

setwd("C:/Users/Elena Ianeva/Desktop/master_thesis/metabolo/Elena Ianeva")
features <- read.csv("features.csv", sep=";")

colnames(features)

rownames(features) <- features$row.ID
features2 <- as.data.frame(t(features[,-c(1:3)]))

MetaData <- read.csv("MetaData.csv", sep=";")
rownames(MetaData) <- MetaData$Filename

df <- merge(MetaData, features2, by = 0)

df2 <- df[-c(1:5), c(9:10, 16:863)]

df2$Treatment <- as.factor(df2$Treatment)
df2$Time_exposed <- as.factor(df2$Time_exposed)
#df2$Wet_weight_mg


names(df2)[1:20]

###
X <- df2[, c(6:ncol(df2))]
Treat <- df2$Treatment
Time <- df2$Time_exposed

# PLS-DA
plsda <- plsda(X,Time, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda <- perf(plsda, validation = 'Mfold', folds = 3, 
                         progressBar = FALSE,  # Set to TRUE to track progress
                         nrepeat = 10)         # We suggest nrepeat = 50

plot(perf.plsda, sd = TRUE, legend.position = 'horizontal') 
perf.plsda$choice.ncomp # ncomp = 3

# tunning (spls-da) _ optimal nb of comp for the model 

tune <- tune.splsda(
  X,
  Time,
 ncomp = 10,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10
)

# optimal number
tune$choice.ncomp$ncomp  # = 3

# FInal model  PLSDA
final_plsda<- plsda(X,Time, ncomp = 3)

# FInal model sPLS-DA
final_splsda <- splsda(
  X,
  Time,
  ncomp = 3,
  keepX = 3
)


plsda_plot <- plotIndiv(final_plsda, ind.names = FALSE, legend=TRUE,
                  comp=c(1,2), ellipse = TRUE)
splsda_plot <- plotIndiv(final_splsda, ind.names = FALSE, legend=TRUE, ellipse = TRUE)

# 3D plot sPLS-DA
scores <- final_splsda$variates$X

df <- data.frame(Comp1 = scores[,1],
                 Comp2 = scores[,2],
                 Comp3 = scores[,3],
                 group = Time)

plot_ly(df,
        x = ~Comp1,
        y = ~Comp2,
        z = ~Comp3,
        color = ~group,
        colors = c("red", "blue", "green"),
        type = "scatter3d",
        mode = "markers")

# PLS-DA with ggplot (Figure S4)

scores <- as.data.frame(final_plsda$variates$X)
scores$group <- Time

scores$group <- factor(
  scores$group,
  levels = c("0h", "24h", "72h", "120h")
)

FigureS4 <- ggplot(scores, aes(x = comp1, y = comp2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_viridis_d(option = "D") +
  labs(
    x = "PLS-DA comp 1",
    y = "PLS-DA comp 2",
    color = "Exposure Time"
  ) +
  theme_minimal()

FigureS4 = FigureS4+ annotate("text",
                                    x = Inf, y = Inf,
                                    label = "CER = 0.21 \n p = 0.001",
                                    hjust = 1.1, vjust = 1.1,
                                    size = 3) ; FigureS4
                                    # Export figure 700 x 420

#RVAideMemoire::MVA.test(X = X, Y = Time, model = "PLS-DA", nperm = 999)
#Permutation test based on cross-validation

#data:  X and Time
#Model: PLS-DA
#8 components
#999 permutations
#CER = 0.21286, p-value = 0.001

# Permanova
vegan::adonis2(dist(X) ~ Time)

dist.mat <- vegdist(X, method = "euclidean")
#default is 999 permutations
res<-pairwiseAdonis::pairwise.adonis(dist.mat, df2[,"Time_exposed"])
res[,3:5] <- round(res[,3:5],2) 
res <- res %>%
  mutate(sig = if_else(p.adjusted < 0.05, "*", ""))
res

TableS2 = gt(res) ; TableS2
# Export figure 715 x 420

##### 2. Main experiment #####

rm(list=ls())

setwd("C:/Users/Elena Ianeva/Desktop/Treatment")

formatted_peak_table <- read.csv("formatted_peak_table.csv", h = T, sep = ",")

formatted_peak_table = as.data.frame(formatted_peak_table)
row.names(formatted_peak_table) = formatted_peak_table$sample_id

# import metadata
setwd("C:/Users/Elena Ianeva/Desktop/master_thesis/metabolo_26")

MetaData <- read.csv("MetaData_Elena_26_.csv", sep=";")
MetaData$X.1 <- NULL
MetaData$X.2 <- NULL
MetaData$X <- NULL

MetaData$correct_name <- sub(".*sample_|\\.mzML$", "", MetaData$correct_name)
MetaData <- MetaData[MetaData$correct_name != "",]
MetaData <- MetaData[ !(MetaData$Plant_species %in% c("M. verticillatum", "M. heterophyllum")), ]
MetaData$correct_name <- sub("\\.mzML$", "", MetaData$correct_name)
MetaData <- MetaData[-134, ] # it's a QC

row.names(MetaData) = MetaData$correct_name


# merge both 
df <- merge(MetaData, formatted_peak_table, by = 0)

df$Treatment <- as.factor(df$Treatment.x)
df$Plant_species <- as.factor(df$Plant_species)

X <- df[, 22:ncol(df)]  # remove metadata, keep only features
X <- as.data.frame(lapply(X, as.numeric)) 
X <- as.matrix(X)
X_scaled = scale(X)

spp.nmds <- metaMDS(X, distance="bray", k=2)

Plant_species <- df$Plant_species
Treat <- df$Treatment

nmds_points_veg <- spp.nmds$points # the coord. of the samples. Each row = one sample (total = 226)
colnames(nmds_points_veg) <- c("NMDS1", "NMDS2")
nmds_points_veg <- as.data.frame(nmds_points_veg)


nmds_points_veg$Sample <- rownames(nmds_points_veg)

nmds_tab_spp <- data.frame(
  Plant_species = Plant_species,
  Treat = Treat,
  nmds_points_veg
) # table with NMDS coord, Plant species and treatment 



# NMDS with ggplot (~ Plant species ) or ~ treatment
nmds_stress <- spp.nmds$stress
sp.hull <- 
  nmds_tab_spp %>% 
  group_by(Plant_species) %>% 
  slice(chull(NMDS1, NMDS2))
centroids_sp <- nmds_tab_spp %>%
  dplyr::group_by(Plant_species) %>%
  dplyr::summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2)
  ) %>%
  dplyr::ungroup()

# ~ Plant species 
p1 = ggplot() + 
  geom_point(data = nmds_points_veg,
             aes(x = NMDS1, y = NMDS2, 
                 colour = Plant_species),
             size = 3) +
  geom_text(data = centroids_sp, 
            aes(x = NMDS1, y = NMDS2, label = Plant_species),
            color = "black") +
 # geom_text(data = nmds_points_veg, # pr vérifier les outliers
                 # aes(x = NMDS1, y = NMDS2,
                  #    label = Sample),
                  #size = 3) +
  scale_colour_manual(values = c(
    "L. grandiflora" = "#00798c",
    "M. aquaticum" = "#d1495b",
    "L. repens" = "#edae49",
    "P. palustris" = "#66a182",
    "M. hippuroides" = "#8d96a3"
  )) +
  scale_fill_manual(values = c(
    "L. grandiflora" = "#00798c",
    "M. aquaticum" = "lightcoral",
    "L. repens" = "#edae49",
    "P. palustris" = "#66a182",
    "M. hippuroides" = "#8d96a3"
  )) +
  annotate("text", x = 0, y = 1.4, 
           label = paste("2d stress =", round(nmds_stress, 3))) +
  geom_polygon(data = sp.hull, 
               aes(x = NMDS1, y = NMDS2, fill = Plant_species, group = Plant_species), 
               alpha = 0.30) +
  labs(
    colour = "Plant species",
    fill = "Plant species"
  )+
  coord_equal() +
  guides(colour = "none") +
  
  theme_bw()

# NMDS with ggplot (~ Treatment )
nmds_stress <- spp.nmds$stress
nmds_tab_subset <- nmds_tab_spp %>%
  dplyr::filter(Treat %in% c("Control", "Herbivory 5 days"))

p2 = ggplot() + 
  geom_point(data = nmds_tab_subset,
             aes(x = NMDS1, y = NMDS2, 
                 , 
                 colour = Treat),
             size = 3) +
  #geom_text(data = centroids, 
  #  aes(x = NMDS1, y = NMDS2, label = Treat),
  # color = "black") +
  scale_colour_manual(values = c(
    "Control" = "#00798c",
    "Herbivory 5 days" = "#8d96a3"
  )) +
  annotate("text", x = 0, y = 1.4, 
           label = paste("2d stress =", round(nmds_stress, 3))) +
  #geom_polygon(data = sp.hull, 
  # aes(x = NMDS1, y = NMDS2, fill = Treat, group = Treat), 
  # alpha = 0.30) +
  labs(color = "Treatment") +
  coord_equal() +
  theme_bw()

# Figure S6
FigureS6 = p1 + p2 + plot_annotation(tag_levels = "A"); FigureS6
# Export figure 890 x 420


# PLS-DA Mechanical Damage Vs Herbivory 5 day

df_subset <- subset(df, Treatment %in% c("Herbivory 5 days", "Mechanical Damage")
                    & Plant_species %in% c("M. aquaticum"))


X <- df_subset[, 22:ncol(df_subset)]  # remove metadata, keep only features
X <- as.data.frame(lapply(X, as.numeric)) 
Y <- as.factor(df_subset$Treatment.x)

# PLS-DA
plsda <- plsda(X,Y, ncomp = 10)

set.seed(30) # For reproducibility with this handbook, remove otherwise
perf.plsda <- perf(plsda, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  # Set to TRUE to track progress
                   nrepeat = 10)         # We suggest nrepeat = 50

plot(perf.plsda, sd = TRUE, legend.position = 'horizontal') 
perf.plsda$choice.ncomp # ncomp = 2

# tunning (spls-da) _ optimal nb of comp for the model 

tune <- tune.splsda(
 X,
 Y,
 ncomp = 10,
 validation = "Mfold",
 folds = 5,
 nrepeat = 10
)
# optimal number
tune$choice.ncomp$ncomp  # = 1

# FInal model  PLSDA
final_plsda<- plsda(X,Y, ncomp = 2)

# FInal model sPLS-DA
final_splsda <- splsda(
  X,
  Y,
  ncomp = 2,
  keepX = 3
)


plsda_plot <- plotIndiv(final_plsda, ind.names = FALSE, legend=TRUE,
                        comp=c(1,2), ellipse = TRUE)

# PLS-DA with ggplot (Figure S5)

scores <- as.data.frame(final_plsda$variates$X)
scores$group <- Y

FigureS5 <- ggplot(scores, aes(x = comp1, y = comp2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_viridis_d(option = "D") +
  labs(
    x = "PLS-DA comp 1",
    y = "PLS-DA comp 2",
    color = "Treatment"
  ) +
  theme_minimal()
FigureS5  = FigureS5 + annotate("text",
                                      x = 20, y = Inf,
                                      label = "CER = 0.01 \np = 0.001",
                                      hjust = 1.1, vjust = 1.1,
                                      size = 3) ; FigureS5 
                                      # Export Figure 600 x 420

#RVAideMemoire::MVA.test(X = X, Y = Y, model = "PLS-DA", nperm = 999)
#Permutation test based on cross-validation

#data:  X and Y
#Model: PLS-DA
#8 components
#999 permutations
#CER = 0.01, p-value = 0.001

# Permanova
vegan::adonis2(dist(X) ~ Y)

dist.mat <- vegdist(X, method = "euclidean")
#default is 999 permutations
#res<-pairwiseAdonis::pairwise.adonis(dist.mat, Y)
#res[,3:5] <- round(res[,3:5],2)  
#res

# Print session information for reproducibility (R version + package versions)
sessionInfo()

#R version 4.4.2 (2024-10-31 ucrt)
#Running under: Windows 10 x64 (build 19045)

#Matrix products: default


#locale:
#  [1] LC_COLLATE=French_Switzerland.utf8  LC_CTYPE=French_Switzerland.utf8   
#[3] LC_MONETARY=French_Switzerland.utf8 LC_NUMERIC=C                       
#[5] LC_TIME=French_Switzerland.utf8    

#time zone: Europe/Paris
#tzcode source: internal

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] patchwork_1.3.2      gt_1.3.0             pairwiseAdonis_0.4.1 cluster_2.1.8.2     
#[5] vegan_2.7-2          permute_0.9-8        viridis_0.6.5        viridisLite_0.4.3   
#[9] plotly_4.10.4        lubridate_1.9.3      forcats_1.0.0        stringr_1.6.0       
#[13] dplyr_1.1.4          purrr_1.2.2          readr_2.1.5          tidyr_1.3.2         
#[17] tibble_3.2.1         tidyverse_2.0.0      mixOmics_6.30.0      ggplot2_4.0.2       
#[21] lattice_0.22-6       MASS_7.3-61         

#loaded via a namespace (and not attached):
#  [1] gtable_0.3.6        ellipse_0.5.0       htmlwidgets_1.6.4   ggrepel_0.9.6      
#[5] mathjaxr_2.0-0      tzdb_0.4.0          numDeriv_2016.8-1.1 crosstalk_1.2.1    
#[9] vctrs_0.6.5         tools_4.4.2         generics_0.1.4      parallel_4.4.2     
#[13] rARPACK_0.11-0      pkgconfig_2.0.3     Matrix_1.7-1        data.table_1.16.2  
#[17] RColorBrewer_1.1-3  S7_0.2.0            lifecycle_1.0.5     compiler_4.4.2     
#[25] yaml_2.3.10         lazyeval_0.2.2      pillar_1.11.1       BiocParallel_1.40.0
#[29] metadat_1.4-0       nlme_3.1-166        RSpectra_0.16-2     tidyselect_1.2.1   
#[33] digest_0.6.39       stringi_1.8.7       reshape2_1.4.5      labeling_0.4.3     
#[37] splines_4.4.2       fastmap_1.2.0       grid_4.4.2          cli_3.6.5          
#[41] metafor_4.8-0       magrittr_2.0.3      dichromat_2.0-0.1   effsize_0.8.1      
#[45] corpcor_1.6.10      withr_3.0.2         scales_1.4.0        timechange_0.3.0   
#[49] httr_1.4.7          matrixStats_1.4.1   igraph_2.1.4        gridExtra_2.3      
#[53] hms_1.1.3           mgcv_1.9-4          rlang_1.2.0         Rcpp_1.0.13-1      
#[57] glue_1.8.0          xml2_1.3.6          jsonlite_1.8.9      rstudioapi_0.18.0  
#[61] R6_2.6.1            plyr_1.8.9          fs_1.6.6 
