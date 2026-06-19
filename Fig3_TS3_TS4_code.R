####### Code for Chemical Richness, Shannon Index and PLS-DA Analysis ####### 


# Code last update: 2026-06-05
# R version 4.4.2 (2024-10-31 ucrt)
# Details of the R session at the end of the script

# Author:
# Elena Ianeva

# Clean R memory
rm(list =ls()) 

#### Load Packages #### 

# Load required packages for data manipulation, statistics, and visualization

library(dplyr) # To modify the data
library(FSA) # For statistical analysis
library(multcompView) # convert multiple comparison test results into compact 
# letter display
library(mixOmics) # for multivariate analysis of omics data
library(viridis) # color palettes
library(ggplot2) # for visualization 
library(patchwork) # to combine multiple ggplot objects
library(ggpubr) # to simply plotting
library(broom) # convert stat model outputs into tidy data frames
library(gt) # to create tables in HTML
library(vegan) # to calculate ecological diversity indices
library(tidyr) #to reshape data between wide and long formats

#### Import data ####

# Set working directory and import feature table 
setwd("C:/Users/Elena Ianeva/Desktop/Treatment")
features <- read.csv("formatted_peak_table.csv")

rownames(features) = features$sample_id
features <- as.data.frame(features[,-c(1:10)])
features$Sample = rownames(features)
features_numeric <- features  # copie

# convert the abundance into num values
features_numeric[, -ncol(features_numeric)] <- lapply(
  features_numeric[, -ncol(features_numeric)], 
  function(x) {
    if(is.factor(x) || is.character(x)) {
      as.numeric(as.character(x))
    } else {
      as.numeric(x)
    }
  }
)

# Import metadata
setwd("C:/Users/Elena Ianeva/Desktop/master_thesis/metabolo_26")
MetaData <-read.csv("MetaData_Elena_26_.csv", sep = ";")

# Remove unused columns and excluse non-taget species/samples
MetaData$X.1 <- NULL
MetaData$X.2 <- NULL
MetaData$X <- NULL

# keep a simple name for the samples
MetaData$correct_name <- sub(".*sample_|\\.mzML$", "", MetaData$correct_name)
MetaData <- MetaData[MetaData$correct_name != "",]
MetaData$correct_name <- sub("\\.mzML$", "", MetaData$correct_name)

# remove samples that will not be analyzed
MetaData <- MetaData[ !(MetaData$Plant_species %in% c("M. verticillatum", 
                                                      "M. heterophyllum")), ]
MetaData <- MetaData[MetaData$correct_name !=
                       "\\2026\\Rasmann\\Elena\\20260204_Elena_QC_1", ]

# keep only the desired columns
MetaData <- MetaData[, -c(1, 2, 3)]

# transform MetaData in a data frame
MetaData <- as.data.frame(MetaData)       
rownames(MetaData) <- NULL 

# define the levels of the Treatment factor
MetaData$Treatment <- factor(MetaData$Treatment,
                             levels = c("Control", "Herbivory 1 day", 
                                        "Herbivory 5 days",
                                        "Mechanical Damage"))

#### Chemical richness #### 

chemical_cols <- features_numeric[, 1:(ncol(features_numeric)-1)]

# Calculate chemical richness per sample
# (number of features detected metabolites per sample)
# Presence defined as abundance > 0
MetaData$Richness <- rowSums(chemical_cols > 0)

# Test differences in chemical richness among treatments 
# within each plant species
# Apply ANOVA + Tukey (normal data) or Kruskal-Wallis + Dunn (non-normal data)

letters_results_chem_richness <- MetaData %>%
  group_by(Plant_species) %>%
  group_modify(~ {
    sp <- .y$Plant_species
    x <- .x
    
    # Test for normality
    normal <- shapiro.test(x$Richness)$p.value > 0.05
    
    if(normal){
      # ANOVA + Tukey
      aov_res <- aov(Richness ~ Treatment, data = x)
      cat("ANOVA result for ", sp," =\n"); 
      print(shapiro.test(x$Richness)$p.value)        
      
      tukey <- TukeyHSD(aov_res)
      cat("Post hoc for", sp, ": "); print(tukey)
      comps <- rownames(tukey$Treatment)
      pvals <- tukey$Treatment[,"p adj"]
      comp_list <- strsplit(comps, "-")
      groups <- unique(unlist(comp_list))
      
      pmat <- matrix(1, nrow=length(groups), ncol=length(groups),
                     dimnames=list(groups, groups))
      
      for(i in seq_along(comp_list)){
        g1 <- trimws(comp_list[[i]][1])
        g2 <- trimws(comp_list[[i]][2])
        pmat[g1,g2] <- pvals[i]
        pmat[g2,g1] <- pvals[i]
      }
      
      letters <- multcompLetters(pmat, threshold=0.05)$Letters
      
      tibble(
        Treatment = names(letters),
        Letters = letters,
        Species = sp
      )
      
    } else {
      # Kruskal-Wallis + Dunn
      dunn <- dunnTest(Richness ~ Treatment, data = x, 
                       method = "bonferroni")$res
      cat("Dunn result for ", sp," =\n"); print(dunn)     
      
      comp_list <- strsplit(dunn$Comparison, " - ")
      groups <- unique(unlist(comp_list))
      
      pmat <- matrix(1, nrow=length(groups), ncol=length(groups),
                     dimnames=list(groups, groups))
      
      for(i in seq_along(comp_list)){
        g1 <- comp_list[[i]][1]
        g2 <- comp_list[[i]][2]
        pmat[g1,g2] <- dunn$P.adj[i]
        pmat[g2,g1] <- dunn$P.adj[i]
      }
      
      letters <- multcompLetters(pmat, threshold=0.05)$Letters
      
      tibble(
        Treatment = names(letters),
        Letters = letters,
        Species = sp
      )
    }
  }) %>%
  ungroup()


####  Shannon index #### 

# Compute Shannon diversity index based on metabolomic profiles
shannon_values <- diversity(chemical_cols, index = "shannon")

MetaData$Shannon <- shannon_values

# Test differences in Shannon diversity among treatments per species
# Apply TUkey or Dunn depending on normality

letters_shannon <- MetaData %>%
  group_by(Plant_species) %>%
  group_modify(~ {
    sp <- .y$Plant_species
    x <- .x
    
    # Test de normalité
    normal <- shapiro.test(x$Shannon)$p.value > 0.05
    
    if(normal){
      # ANOVA + Tukey
      aov_res <- aov(Shannon ~ Treatment, data = x)
      cat("ANOVA result for ", sp," =\n");
      print(shapiro.test(x$Shannon)$p.value)        
      
      tukey <- TukeyHSD(aov_res)
      cat("Post hoc for", sp, ": "); print(tukey)
      comps <- rownames(tukey$Treatment)
      pvals <- tukey$Treatment[,"p adj"]
      comp_list <- strsplit(comps, "-")
      groups <- unique(unlist(comp_list))
      
      pmat <- matrix(1, nrow=length(groups), ncol=length(groups),
                     dimnames=list(groups, groups))
      
      for(i in seq_along(comp_list)){
        g1 <- trimws(comp_list[[i]][1])
        g2 <- trimws(comp_list[[i]][2])
        pmat[g1,g2] <- pvals[i]
        pmat[g2,g1] <- pvals[i]
      }
      
      letters <- multcompLetters(pmat, threshold=0.05)$Letters
      
      tibble(
        Treatment = names(letters),
        Letters = letters,
        Species = sp
      )
      
    } else {
      # Kruskal-Wallis + Dunn
      dunn <- dunnTest(Shannon ~ Treatment, data = x, method = "bonferroni")$res
      cat("Dunn result for ", sp," =\n"); print(dunn)     
      
      comp_list <- strsplit(dunn$Comparison, " - ")
      groups <- unique(unlist(comp_list))
      
      pmat <- matrix(1, nrow=length(groups), ncol=length(groups),
                     dimnames=list(groups, groups))
      
      for(i in seq_along(comp_list)){
        g1 <- comp_list[[i]][1]
        g2 <- comp_list[[i]][2]
        pmat[g1,g2] <- dunn$P.adj[i]
        pmat[g2,g1] <- dunn$P.adj[i]
      }
      
      letters <- multcompLetters(pmat, threshold=0.05)$Letters
      
      tibble(
        Treatment = names(letters),
        Letters = letters,
        Species = sp
      )
    }
  }) %>%
  ungroup()


df_long <- MetaData %>%
  filter(Treatment %in% c("Control", "Herbivory 5 days")) %>%
  pivot_longer(cols = c(Shannon, Richness),
               names_to = "Index",
               values_to = "Value")
df_long$Plant_species <- factor(df_long$Plant_species,
                                levels = c("M. aquaticum", 
                                           setdiff(unique
                                                   (df_long$Plant_species),
                                                   "M. aquaticum"))
)

df_long <- df_long %>%
  mutate(Index = case_when(
    Index == "Shannon" ~ "Shannon Index",
    Index == "Richness" ~ "Chemical Richness",
    TRUE ~ Index
  ))

letters_results_chem_richness$Index <- "Chemical Richness"

letters_shannon$Index  <- "Shannon Index"

# Create a data frame with the significance letters
letters_all <- bind_rows(letters_results_chem_richness,letters_shannon)
letters_all$Index <- factor(letters_all$Index,
                            levels = c("Shannon Index", "Chemical Richness")
)

letters_all$Plant_species <- factor(letters_all$Plant_species,
                                    levels = levels(df_long$Plant_species)
)

# keep only control and Herbivory 5 days conditions (as there other treatments
# in M. aquaticum which are not analyzed here)

letters_all <- letters_all %>%
  filter(Treatment %in% c("Control", "Herbivory 5 days"))


y_max <- df_long %>%
  dplyr::group_by(Plant_species, Treatment, Index) %>%
  dplyr::summarise(
    y_max = max(Value, na.rm = TRUE),
    .groups = "drop")

y_pos = y_max$y_max * 1.05

letters_all <- letters_all %>%
  left_join(y_max, by = c("Plant_species", "Treatment", "Index")) %>%
  mutate(y_pos = y_max * 1.05)


letters_all <- letters_all %>%
  mutate(
    Letters = ifelse(Plant_species == "M. aquaticum", "a", Letters)
  )

####  PLS-DA #### 

# Perform Partial Least Squares Discriminant Analysis 
# Evaluate separation of treatment grops within each species

setwd("C:/Users/Elena Ianeva/Desktop/Treatment")
features_2 <- read.csv("formatted_peak_table.csv")

features_2$Treatment[features_2$Treatment == "Induced"] <- "Herbivory 5 days"
features_2$Treatment[features_2$Treatment == 
                       "Herbivory_5_days"] <- "Herbivory 5 days"

# number of the column when the metabolites variables starts
col_metadata <- 10
# number of the column when the meadata stops
col_feature_table <- col_metadata + 1
features_2[, col_feature_table:ncol(features_2)] <- 
  lapply(features_2[, col_feature_table:ncol(features_2)], as.numeric)

###### Myriophyllum aquaticum ###### 

metabo_PF_df  <- subset(features_2, ATTTRIBUTE_species == "M._aquaticum")
metabo_PF <- metabo_PF_df[, (col_metadata + 1):ncol(metabo_PF_df)]
metabo_PF <- metabo_PF[, colSums(metabo_PF, na.rm = TRUE) != 0]

keep <- metabo_PF_df$Treatment %in% c("Control", "Herbivory 5 days")

metabo_PF_sub <- metabo_PF[keep, ]
metabo_PF_df_sub <- metabo_PF_df[keep, ]

metabo_PF_df_sub$Treatment <- factor(metabo_PF_df_sub$Treatment)

plsda_PF <- plsda(
  X = as.matrix(metabo_PF_sub),
  Y = metabo_PF_df_sub$Treatment,
  ncomp = 2
)

# test if the metabolomic profiles differ between the treatments 

#RVAideMemoire::MVA.test(X = as.matrix(metabo_PF_sub), 
#Y = (metabo_PF_df_sub$Treatment), model = "PLS-DA", nperm = 999) 

#Permutation test based on cross-validation

#data:  as.matrix(metabo_PF_sub) and (metabo_PF_df_sub$Treatment)
#Model: PLS-DA
#8 components
#999 permutations
#CER = 0.34643, p-value = 0.026

###### Myriophyllum hippuroides  ###### 
metabo_hipp_df  <- subset(features_2, ATTTRIBUTE_species == "M._hippuroides")
metabo_hipp <- metabo_hipp_df[, (col_metadata + 1):ncol(metabo_hipp_df)]
metabo_hipp <- metabo_hipp[, colSums(metabo_hipp, na.rm = TRUE) != 0]

plsda_hipp  <- plsda(metabo_hipp, metabo_hipp_df$Treatment, ncomp = 2)
metabo_hipp_df$Treatment <- as.factor(metabo_hipp_df$Treatment)
#RVAideMemoire::MVA.test(X = as.matrix(metabo_hipp), 
#Y = (metabo_hipp_df$Treatment), model = "PLS-DA", nperm = 999) 

#Permutation test based on cross-validation

#data:  as.matrix(metabo_hipp) and (metabo_hipp_df$Treatment)
#Model: PLS-DA
#8 components
#999 permutations
#CER = 0.3413, p-value = 0.026

###### Ludwigia grandiflora  ###### 

metabo_grandi_df  <- subset(features_2, ATTTRIBUTE_species == "L._grandiflora")
metabo_grandi <- metabo_grandi_df[, (col_metadata + 1):ncol(metabo_grandi_df)]
metabo_grandi <- metabo_grandi[, colSums(metabo_grandi, na.rm = TRUE) != 0]

plsda_grandi  <- plsda(metabo_grandi, metabo_grandi_df$Treatment, ncomp = 2)
metabo_grandi_df$Treatment <- as.factor(metabo_grandi_df$Treatment)

#RVAideMemoire::MVA.test(X = as.matrix(metabo_grandi),
#Y = (metabo_grandi_df$Treatment), model = "PLS-DA", nperm = 999) 

#Permutation test based on cross-validation

#data:  as.matrix(metabo_grandi) and (metabo_grandi_df$Treatment)
#Model: PLS-DA
#8 components
#999 permutations
#CER = 0.34167, p-value = 0.116

###### Laurembergia repens  ###### 

metabo_repens_df <- subset(features_2, ATTTRIBUTE_species == "L._repens")
metabo_repens <- metabo_repens_df[, (col_metadata + 1):ncol(metabo_repens_df)]
metabo_repens <- metabo_repens[, colSums(metabo_repens, na.rm = TRUE) != 0]

plsda_repens  <- plsda(metabo_repens, metabo_repens_df$Treatment, ncomp = 2)
metabo_repens_df$Treatment <- as.factor(metabo_repens_df$Treatment)

#RVAideMemoire::MVA.test(X = as.matrix(metabo_repens), 
#Y = (metabo_repens_df$Treatment), model = "PLS-DA", nperm = 999) 

#Permutation test based on cross-validation

#data:  as.matrix(metabo_repens) and (metabo_repens_df$Treatment)
#Model: PLS-DA
#8 components
#999 permutations
#CER = 0.27708, p-value = 0.033


###### Proserpinaca palustris   ###### 

metabo_palu_df <- subset(features_2, ATTTRIBUTE_species == "P._palustris")
metabo_palu<- metabo_palu_df[, (col_metadata + 1):ncol(metabo_palu_df)]
metabo_palu <- metabo_palu[, colSums(metabo_palu, na.rm = TRUE) != 0]

plsda_palu  <- plsda(metabo_palu, metabo_palu_df$Treatment, ncomp = 2)
metabo_palu_df$Treatment <- as.factor(metabo_palu_df$Treatment)

#RVAideMemoire::MVA.test(X = as.matrix(metabo_palu), 
#Y = (metabo_palu_df$Treatment), model = "PLS-DA", nperm = 999) 

#Permutation test based on cross-validation

#data:  as.matrix(metabo_palu) and (metabo_palu_df$Treatment)
#Model: PLS-DA
#8 components
#999 permutations
#CER = 0.16875, p-value = 0.003


#### Function to plot the plsda with ggplot #### 
plot_plsda <- function(plsda_model, group,sample_names = NULL, title = NULL) {
  
  indiv_df <- data.frame(
    Comp1 = plsda_model$variates$X[, 1],
    Comp2 = plsda_model$variates$X[, 2],
    Group = group
  )
  
  p <- ggplot(indiv_df, aes(x = Comp1, y = Comp2, color = Group, 
                            shape = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(type = "norm", level = 0.95, linewidth = 1) +
    scale_color_viridis_d(option = "D") +
    theme_minimal() +
    theme(legend.position = "left")
  
  return(p)
}

plot_PLSDA_PF <- plot_plsda(
  plsda_model = plsda_PF,
  sample_names = rownames(metabo_PF_sub),
  group = metabo_PF_df_sub$Treatment,
  # title = "PLS-DA M. aquaticum"
) + theme(legend.position = "none")

plot_PLSDA_MH <- plot_plsda(
  plsda_model = plsda_hipp,
  sample_names = rownames(metabo_hipp),
  group = metabo_hipp_df$Treatment,
  # title = "PLS-DA M. hippuroides"
) + theme(legend.position = "none")

plot_PLSDA_LG <- plot_plsda(
  plsda_model = plsda_grandi,
  sample_names = rownames(metabo_grandi),
  group = metabo_grandi_df$Treatment,
  #title = "PLS-DA L. grandiflora"
) + theme(legend.position = "none")

plot_PLSDA_LR <- plot_plsda(
  plsda_model = plsda_repens,
  sample_names = rownames(metabo_repens),
  group = metabo_repens_df$Treatment,
  # title = "PLS-DA L. repens"
) + theme(legend.position = "none")

plot_PLSDA_PP <- plot_plsda(
  plsda_model = plsda_palu,  
  sample_names = rownames(metabo_palu),
  group = metabo_palu_df$Treatment,
  # title = "PLS-DA P. palustris"
) 

# Plot Chemical Richness, Shannon Index and PLS-DA (Figure 2)

species <- c(
  "M. aquaticum",
  "M. hippuroides",
  "L. grandiflora",
  "L. repens",
  "P. palustris"
)


plot_A_fun <- function(sp) {
  ggplot(subset(df_long,
                Index == "Chemical Richness" &
                  Plant_species == sp),
         aes(x = Treatment, y = Value, fill = Treatment)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    scale_fill_viridis_d(option = "D") +
    geom_text(
      data = subset(letters_all,
                    Index == "Chemical Richness" &
                      Plant_species == sp),
      aes(x = Treatment, y = y_pos, label = Letters),
      inherit.aes = FALSE, size = 3
    ) +
    theme_bw() +
    labs(title = NULL, x = NULL, y = "Chemical Richness") +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 8),
      axis.text.y  = element_text(size = 7)
    )
}

plot_B_fun <- function(sp) {
  ggplot(subset(df_long,
                Index == "Shannon Index" &
                  Plant_species == sp),
         aes(x = Treatment, y = Value, fill = Treatment)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    scale_fill_viridis_d(option = "D") +
    geom_text(
      data = subset(letters_all,
                    Index == "Shannon Index" &
                      Plant_species == sp),
      aes(x = Treatment, y = y_pos, label = Letters),
      inherit.aes = FALSE, size = 3
    ) +
    theme_bw() +
    labs(title = NULL, x = NULL, y = "Shannon Index") +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 8),
      axis.text.y  = element_text(size = 7)
    )
}

title_fun <- function(label){
  ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = label,
             fontface = "italic",
             size = 4) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = ggplot2::margin(
        t = 5.5, r = 5.5, b = 5.5, l = 30
      )
    )
}

header_row <- ggarrange(
  plotlist = lapply(species, title_fun),
  ncol = 5
)
row_A <- ggarrange(
  plotlist = lapply(species, plot_A_fun),
  ncol = 5
)

row_B <- ggarrange(
  plotlist = lapply(species, plot_B_fun),
  ncol = 5
)

row_C_plots <- ggarrange(
  plot_PLSDA_PF,
  plot_PLSDA_MH,
  plot_PLSDA_LG,
  plot_PLSDA_LR,
  plot_PLSDA_PP,
  ncol = 5
)
legend <- get_legend(
  plot_PLSDA_PP + theme(legend.position = "right")
)
row_C <- ggarrange(
  plot_PLSDA_PF + theme(legend.position = "none"),
  plot_PLSDA_MH + theme(legend.position = "none"),
  plot_PLSDA_LG + theme(legend.position = "none"),
  plot_PLSDA_LR + theme(legend.position = "none"),
  plot_PLSDA_PP + theme(legend.position = "none"),  
  ncol = 5
)

main_fig <- ggarrange(
  header_row,
  row_A,
  row_B,
  row_C,
  ncol = 1,
  heights = c(0.2, 1, 1, 1.2),
  labels = c("", "A", "B", "C")
)

bottom_row <- ggarrange(
  legend,
  ncol = 1
)

final_fig <- ggarrange(
  main_fig,
  bottom_row,
  ncol = 1,
  heights = c(10, 1.2)
)

final_fig # Export 1200 x 600



#### Summarize statistical tests for supplementary materials #### 
# Include test statistics, p-values, and post-hoc comparisons

stats_summary <- MetaData %>%
  group_by(Plant_species) %>%
  group_modify(~ {
    
    x <- .x
    
    # NORMALITY TEST 
    normal <- shapiro.test(x$Richness)$p.value > 0.05
    
    if(normal){
      
      #ANOVA
      aov_res <- aov(Richness ~ Treatment, data = x)
      aov_tab <- broom::tidy(aov_res)
      
      #  TUKEY 
      tukey <- TukeyHSD(aov_res)$Treatment
      
      sig_pairs <- rownames(tukey)[tukey[, "p adj"] < 0.05]
      sig_pvals_num <- tukey[tukey[, "p adj"] < 0.05, "p adj"]
      
      sig_pvals_text <- ifelse(
        sig_pvals_num < 0.05,
        "< 0.05",
        "> 0.05"
      )
      
      tibble(
        Test = "ANOVA",
        Df = aov_tab$df[1],
        Statistic = round(aov_tab$statistic[1], 3),
        P_value = ifelse(aov_tab$p.value[1] < 0.05,
                         "< 0.05",
                         "> 0.05"),
        Post_hoc = "Tukey",
        Significant_comparisons =
          ifelse(length(sig_pairs) == 0,
                 "None",
                 paste(sig_pairs, collapse = "; ")),
        Posthoc_p_values =
          ifelse(length(sig_pvals_text) == 0,
                 "None",
                 paste(sig_pvals_text, collapse = "; "))
      )
      
    } else {
      
      #  KRUSKAL-WALLIS 
      kw <- kruskal.test(Richness ~ Treatment, data = x)
      
      #  DUNN 
      dunn <- dunnTest(
        Richness ~ Treatment,
        data = x,
        method = "bonferroni"
      )$res
      
      sig_pairs <- dunn$Comparison[dunn$P.adj < 0.05]
      sig_pvals_num <- dunn$P.adj[dunn$P.adj < 0.05]
      
      sig_pvals_text <- ifelse(
        sig_pvals_num < 0.05,
        "< 0.05",
        "> 0.05"
      )
      
      tibble(
        Test = "Kruskal-Wallis",
        Df = kw$parameter,
        Statistic = round(as.numeric(kw$statistic), 3),
        P_value = ifelse(kw$p.value < 0.05,
                         "< 0.05",
                         "> 0.05"),
        Post_hoc = "Dunn",
        Significant_comparisons =
          ifelse(length(sig_pairs) == 0,
                 "None",
                 paste(sig_pairs, collapse = "; ")),
        Posthoc_p_values =
          ifelse(length(sig_pvals_text) == 0,
                 "None",
                 paste(sig_pvals_text, collapse = "; "))
      )
    }
    
  }) %>%
  ungroup()

#  GT TABLE 

stats_summary %>%
  gt() %>%
  
  tab_header(
    title = "Statistical analysis of chemical richness by plant species"
  ) %>%
  
  cols_label(
    Plant_species = "Plant species",
    Test = "Test used",
    Df = "Degrees of freedom",
    Statistic = "F / χ² value",
    P_value = "Test p-value",
    Post_hoc = "Post-hoc test",
    Significant_comparisons = "Significant comparisons",
    Posthoc_p_values = "Post-hoc p-values"
  ) %>%
  
  # Italic species names
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = Plant_species)
  ) %>%
  
  # Bold significant global p-values
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = P_value,
      rows = P_value == "< 0.05"
    )
  ) %>%
  
  # Bold significant posthoc p-values
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Posthoc_p_values,
      rows = Posthoc_p_values == "< 0.05"
    )
  )


stats_summary <- MetaData %>%
  group_by(Plant_species) %>%
  group_modify(~ {
    
    x <- .x
    
    # NORMALITY TEST 
    normal <- shapiro.test(x$Shannon)$p.value > 0.05
    
    if(normal){
      
      #ANOVA
      aov_res <- aov(Shannon ~ Treatment, data = x)
      aov_tab <- broom::tidy(aov_res)
      
      #  TUKEY 
      tukey <- TukeyHSD(aov_res)$Treatment
      
      sig_pairs <- rownames(tukey)[tukey[, "p adj"] < 0.05]
      sig_pvals_num <- tukey[tukey[, "p adj"] < 0.05, "p adj"]
      
      sig_pvals_text <- ifelse(
        sig_pvals_num < 0.05,
        "< 0.05",
        "> 0.05"
      )
      
      tibble(
        Test = "ANOVA",
        Df = aov_tab$df[1],
        Statistic = round(aov_tab$statistic[1], 3),
        P_value = ifelse(aov_tab$p.value[1] < 0.05,
                         "< 0.05",
                         "> 0.05"),
        Post_hoc = "Tukey",
        Significant_comparisons =
          ifelse(length(sig_pairs) == 0,
                 "None",
                 paste(sig_pairs, collapse = "; ")),
        Posthoc_p_values =
          ifelse(length(sig_pvals_text) == 0,
                 "None",
                 paste(sig_pvals_text, collapse = "; "))
      )
      
    } else {
      
      #  KRUSKAL-WALLIS 
      kw <- kruskal.test(Shannon ~ Treatment, data = x)
      
      #  DUNN 
      dunn <- dunnTest(
        Shannon ~ Treatment,
        data = x,
        method = "bonferroni"
      )$res
      
      sig_pairs <- dunn$Comparison[dunn$P.adj < 0.05]
      sig_pvals_num <- dunn$P.adj[dunn$P.adj < 0.05]
      
      sig_pvals_text <- ifelse(
        sig_pvals_num < 0.05,
        "< 0.05",
        "> 0.05"
      )
      
      tibble(
        Test = "Kruskal-Wallis",
        Df = kw$parameter,
        Statistic = round(as.numeric(kw$statistic), 3),
        P_value = ifelse(kw$p.value < 0.05,
                         "< 0.05",
                         "> 0.05"),
        Post_hoc = "Dunn",
        Significant_comparisons =
          ifelse(length(sig_pairs) == 0,
                 "None",
                 paste(sig_pairs, collapse = "; ")),
        Posthoc_p_values =
          ifelse(length(sig_pvals_text) == 0,
                 "None",
                 paste(sig_pvals_text, collapse = "; "))
      )
    }
    
  }) %>%
  ungroup()

#  GT TABLE 

stats_summary %>%
  gt() %>%
  
  tab_header(
    title = "Statistical analysis of Shannon Index by plant species"
  ) %>%
  
  cols_label(
    Plant_species = "Plant species",
    Test = "Test used",
    Df = "Degrees of freedom",
    Statistic = "F / χ² value",
    P_value = "Test p-value",
    Post_hoc = "Post-hoc test",
    Significant_comparisons = "Significant comparisons",
    Posthoc_p_values = "Post-hoc p-values"
  ) %>%
  
  # Italic species names
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = Plant_species)
  ) %>%
  
  # Bold significant global p-values
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = P_value,
      rows = P_value == "< 0.05"
    )
  ) %>%
  
  # Bold significant posthoc p-values
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Posthoc_p_values,
      rows = Posthoc_p_values == "< 0.05"
    )
  )


# Print session information for reproducibility (R version + package versions)
sessionInfo()

#R version 4.4.2 (2024-10-31 ucrt)
#Platform: x86_64-w64-mingw32/x64
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
#  [1] tidyr_1.3.2         vegan_2.7-2         permute_0.9-8       gt_1.3.0           
#[5] broom_1.0.7         ggpubr_0.6.0        patchwork_1.3.2     viridis_0.6.5      
#[9] viridisLite_0.4.3   mixOmics_6.30.0     ggplot2_4.0.2       lattice_0.22-6     
#[13] MASS_7.3-61         multcompView_0.1-10 FSA_0.10.1          dplyr_1.1.4        

#loaded via a namespace (and not attached):
#  [1] gtable_0.3.6        ellipse_0.5.0       ggrepel_0.9.6       rstatix_0.7.2      
#[5] vctrs_0.6.5         tools_4.4.2         generics_0.1.4      parallel_4.4.2     
#[9] tibble_3.2.1        cluster_2.1.8.2     rARPACK_0.11-0      pkgconfig_2.0.3    
#[13] Matrix_1.7-1        RColorBrewer_1.1-3  S7_0.2.0            lifecycle_1.0.5    
#[17] compiler_4.4.2      farver_2.1.2        stringr_1.6.0       codetools_0.2-20   
#[21] carData_3.0-5       sass_0.4.9          htmltools_0.5.8.1   Formula_1.2-5      
#[25] pillar_1.11.1       car_3.1-3           BiocParallel_1.40.0 dunn.test_1.3.7    
#[29] abind_1.4-8         nlme_3.1-166        RSpectra_0.16-2     tidyselect_1.2.1   
#[33] digest_0.6.39       stringi_1.8.7       reshape2_1.4.5      purrr_1.2.2        
#[37] labeling_0.4.3      splines_4.4.2       cowplot_1.1.3       fastmap_1.2.0      
#[41] grid_4.4.2          cli_3.6.5           magrittr_2.0.3      dichromat_2.0-0.1  
#[45] corpcor_1.6.10      withr_3.0.2         scales_1.4.0        backports_1.5.1    
#[49] matrixStats_1.4.1   igraph_2.1.4        gridExtra_2.3       ggsignif_0.6.4     
#[53] mgcv_1.9-4          rlang_1.2.0         Rcpp_1.0.13-1       glue_1.8.0         
#[57] xml2_1.3.6          rstudioapi_0.18.0   R6_2.6.1            plyr_1.8.9         
#[61] fs_1.6.6 
