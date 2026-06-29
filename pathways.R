####### Code for metabolomic analyses (pathway regulation) ####### 


# Code last update: 2026-06-29
# R version 4.4.2 (2024-10-31 ucrt)
# Details of the R session at the end of the script

# Author:
# Elena Ianeva

# Clean R memory
rm(list =ls()) 

#### Load Packages #### 
library(tidyverse) # Edit the data
library(dplyr) # To modify the data
library(ggplot2) # To make the plots
library(tidyr) # To reshape data 
library(effsize) # To calculate effect sizes
library(metafor) # Meta-analysis
library(gt) # To create tables in HTML
library(purrr) # To apply functions repeatedly


#### Import data #### 
setwd("C:/Users/Elena Ianeva/Desktop/Treatment")
metabo <- read.csv("formatted_peak_table.csv", h = T, sep = ",")
#metabo$Treatment  = as.factor(metabo$Treatment)

metabo$Treatment[metabo$Treatment == "Induced"] <- "Herbivory 5 days"
metabo$Treatment[metabo$Treatment == "Herbivory_1_day"] <-
  "Herbivory 1 day (plant material collected after 120h)"
metabo$Treatment[metabo$Treatment == "Herbivory_5_days"] <- "Herbivory 5 days"
metabo$Treatment[metabo$Treatment == "Mechanical_Damage"] <- "Mechanical Damage"

metabo$ATTTRIBUTE_species[metabo$ATTTRIBUTE_species == "L._grandiflora"] <-
  "L. grandiflora"
metabo$ATTTRIBUTE_species[metabo$ATTTRIBUTE_species =="L._repens"] <-
  "L. repens"
metabo$ATTTRIBUTE_species[metabo$ATTTRIBUTE_species == "M._aquaticum"] <- 
  "M. aquaticum"
metabo$ATTTRIBUTE_species[metabo$ATTTRIBUTE_species == "M._hippuroides"] <- 
  "M. hippuroides"
metabo$ATTTRIBUTE_species[metabo$ATTTRIBUTE_species == "P._palustris"] <- 
  "P. palustris"


rownames(metabo) <- NULL
annotations <- read.csv("formatted_annotation.csv", h = T, sep = ";")    
# Sirius annotations

#### Exploration #### 
# Look at the data
head(metabo[,1:30])
head(annotations)

# number of the column when the metabolites variables starts
col_metadata <- 10
# number of the column when the meadata stops
col_feature_table <- col_metadata + 1

# Keep only the feature number as colnames 
names_feature <- colnames(metabo[,col_feature_table:length(metabo[1,])])
names_feature <- gsub("feature", "", names_feature)
names_feature <- separate(as.data.frame(names_feature), 
                          names_feature,
                          into = c("feature", "trash1", "trash2"), 
                          sep = "_")[,1]
names_feature <- c(colnames(metabo[,1:col_metadata]), names_feature)
colnames(metabo) <- names_feature

# Pathways
# Check the annotation quality and keep only feature with a score of > 0.7
annotations$pathways <- annotations$feature_pred_tax_npc_01pat_val
annotations["pathways"][annotations["feature_pred_tax_npc_01pat_score"] < 0.7] <-"Unannotated"
table(annotations["pathways"]) # number of features per pathway

# Barplot
# Count number of features per pathway
pathway_counts <- table(annotations$pathways)  # or annotations["pathways"]
pathway_counts <- as.data.frame(pathway_counts)
pathway_counts$Var1 <- factor(
  pathway_counts$Var1,
  levels = pathway_counts$Var1[order(pathway_counts$Freq, decreasing = TRUE)]
)
colnames(pathway_counts) <- c("Pathway", "Count")

# Plot
ggplot(pathway_counts, aes(x = Pathway, y = Count)) +
  geom_bar(stat = "identity", fill = "#00C19A", color = "black", 
           linewidth = 0.8) +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Pathway", y = "Number of features", 
       title = "Feature counts per pathway", 
       caption = "Only features with a pathway annotation score > 0.7 are shown")+
  scale_y_continuous(breaks = seq(0, 3000, by = 1000), limits = c(0, 3000))
# Export figure 842 x 544

# Superclasses

# Check the annotation quality and keep only feature with a score of > 0.7
annotations$superclasses <- annotations$feature_pred_tax_npc_02sup_val
annotations["superclasses"][annotations["feature_pred_tax_npc_02sup_score"] < 0.7] <-"Unannotated"
table(annotations["superclasses"])

# Barplot
# Count number of features per superclass
superclass_counts <- table(annotations$superclasses)  # or annotations["pathways"]
superclass_counts <- as.data.frame(superclass_counts)

colnames(superclass_counts) <- c("Superclass", "Count")
superclass_counts$Superclass <- factor(
  superclass_counts$Superclass,
  levels = superclass_counts$Superclass[order(superclass_counts$Count, 
                                              decreasing = TRUE)]
)
target_class <- "Unannotated"

ggplot(subset(superclass_counts, Superclass != target_class), 
       # without unannotated because scale pb
       aes(Superclass, Count)) +
  geom_col(fill = "#edae49", color = "black") +
  theme_minimal() +
  labs(x = "Superclass", 
       y = "Number of features",
       title = "Feature counts per superclass", 
       caption = "Only features with annotation score > 0.7 are shown") +
  coord_cartesian(ylim = c(0, 150)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
# Export figure 1379 x 574

# Classes

# Check the annotation quality and keep only feature with a score of > 0.7
annotations$classes <- annotations$feature_pred_tax_npc_03cla_val
annotations["classes"][annotations["feature_pred_tax_npc_03cla_score"] < 0.7] <-"Unannotated"
table(annotations["classes"])

# Order "metabo"
metabo <- metabo[order(metabo$sample_id),]

# Barplot
class_counts <- table(annotations$classes)  # or annotations["pathways"]
class_counts <- as.data.frame(class_counts)
colnames(class_counts) <- c("Class", "Count")
class_counts$Class <- factor(
  class_counts$Class,
  levels = class_counts$Class[order(class_counts$Count, decreasing = TRUE)]
)

target_class <- "Unannotated"

ggplot(subset(class_counts, Class != target_class), 
       # without unannotated because scale pb
       aes(Class, Count)) +
  geom_col(fill = "#00798c", color = "black") +
  theme_minimal() +
  labs(x = "Class", 
       y = "Number of features", 
       title = "Feature counts per class", 
       caption = "Only features with annotation score > 0.7 are shown") +
  coord_cartesian(ylim = c(0, 120)) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
# Export figure 1156 x 550


#### Effect size #### 

col_feature_table <- 11
col_metadata <- 10



df_annotation_treatment <- function(metabo, annotations, col_annoation, 
                                    col_feature_table, col_metadata, 
                                    levels = NULL, my_theme = theme_minimal()) {  
  cohen_tot <- rbind()
  treat_levels <- levels(as.factor(metabo$Treatment))
  
  label1 <- treat_levels[1]
  label2 <- treat_levels[2]
  
  print(label1)
  
  for(i_sp in unique(metabo$ATTTRIBUTE_species)){
    
    metabo_i_sp <- subset(metabo, ATTTRIBUTE_species == i_sp)
    
    # keep only the matrix-like part of the df which are needed 
    matrix_superclass <- metabo_i_sp[,col_feature_table:ncol(metabo_i_sp)]
    matrix_superclass <- matrix_superclass %>% select_if(colSums(.) != 0) 
    # Remove all the columns without any molecules
    matrix_superclass_scaled <- as.data.frame(cbind(metabo_i_sp[,1:col_metadata], 
                                                    scale((matrix_superclass))))
    
    df_scaled <- gather(matrix_superclass_scaled, 
                        key = col_annoation, 
                        value = "Expression", 
                        all_of(colnames(matrix_superclass)))
    head(df_scaled)
    length(df_scaled$filename)
    
    # Initialize an empty vector to store Cohen's d results
    cohen_d_col_annoation <- c()
    
    # Loop through each elevation
    for (el in unique(df_scaled$col_annoation)) {
      df_scaled_subset <- subset(df_scaled, col_annoation == el)
      
      
      # Calculate Cohen's d effect size
      es <- effsize::cohen.d(df_scaled_subset$Expression,
                             df_scaled_subset$Treatment, 
                             hedges.correction = TRUE, 
                             conf.level = 0.95, na.rm = TRUE)
      
      # Extract the mean and confidence intervals
      mean <- as.numeric(es$estimate)
      inf <- as.numeric(es$conf.int[1])
      sup <- as.numeric(es$conf.int[2])
      
      # Store the results in a list
      cohen_d_col_annoation <- rbind(cohen_d_col_annoation, 
                                     c(el, mean, inf, sup))
      
    }
    
    
    # Convert the list to a data frame
    cohen_d_col_annoation <- as.data.frame(cohen_d_col_annoation)
    colnames(cohen_d_col_annoation) <- c("col_annoation", "cohen_d", "CI_inf", "CI_sup")
    
    # Convert columns to appropriate data types
    cohen_d_col_annoation$cohen_d <- as.numeric(cohen_d_col_annoation$cohen_d)
    cohen_d_col_annoation$CI_inf <- as.numeric(cohen_d_col_annoation$CI_inf)
    cohen_d_col_annoation$CI_sup <- as.numeric(cohen_d_col_annoation$CI_sup)
    
    # Claculate 95% CI and SE:
    cohen_d_col_annoation$CI <- cohen_d_col_annoation$CI_sup - cohen_d_col_annoation$cohen_d
    cohen_d_col_annoation$SE <- cohen_d_col_annoation$CI/1.96
    
    annotation_2 <- annotations[annotations$feature_id %in% colnames(matrix_superclass),]
    
    # 🔎 CHECK ORDER BEFORE ASSIGNMENT
    all(annotation_2$feature_id == cohen_d_col_annoation$col_annoation)
    
    # ❗ STOP if mismatch (recommended while debugging)
    if(!all(annotation_2$feature_id == cohen_d_col_annoation$col_annoation)){
      stop("Feature order mismatch between annotations and Cohen's d table")
    }
    
    annotation_2$cohen_d <- cohen_d_col_annoation$cohen_d
    annotation_2$var_d <- ((cohen_d_col_annoation$CI_sup - cohen_d_col_annoation$CI_inf) / (2 * 1.96))^2
    
    pooled_ES <- c()
    for(annot in unique(annotation_2[,col_annoation])){
      df_annotation <- annotation_2[annotation_2[,col_annoation] == annot,]
      
      res <- metafor::rma(yi = df_annotation$cohen_d, 
                          vi = df_annotation$var_d, 
                          method = "REML") # random-effects
      parameters_annot <- c(res$b[1], res$ci.lb, res$ci.ub, annot, i_sp)
      pooled_ES <- rbind(pooled_ES, parameters_annot)
    }
    
    pooled_ES <- as.data.frame(pooled_ES)
    colnames(pooled_ES) <- c("ES", "CI_sup", "CI_inf", 
                             col_annoation, "species")
    pooled_ES$ES <- as.numeric(pooled_ES$ES)
    pooled_ES$CI_inf <- as.numeric(pooled_ES$CI_inf)
    pooled_ES$CI_sup <- as.numeric(pooled_ES$CI_sup)
    pooled_ES$CI <- (pooled_ES$CI_inf-pooled_ES$CI_sup)/2
    pooled_ES$Effects <- paste(round(pooled_ES$ES, 2), "±", 
                               round(pooled_ES$CI, 2))
    
    # Display the final data frame
    print_ES <- pooled_ES %>% dplyr::select(species, col_annoation, Effects) %>% dplyr::arrange(col_annoation)
    print(print_ES)
    
    expression <- c()
    significance <- c()
    for(i in 1:nrow(pooled_ES)){
      sign <- ifelse(pooled_ES$CI_sup[i] > 0 | pooled_ES$CI_inf[i] < 0, 
                     "Significant", "Not significant")
      expr <- ifelse(
        pooled_ES$ES[i] < 0,
        paste("Overexpressed in", label2),
        paste("Overexpressed in", label1)
      )
      significance <- c(significance, sign)
      expression <- c(expression, expr)
    }
    
    pooled_ES$expression <- expression
    pooled_ES$significance <- significance
    
    cohen_tot <- rbind(cohen_tot, pooled_ES)
    
  }
  
  cohen_tot <- as.data.frame(cohen_tot)
  #cohen_tot <- cohen_tot[cohen_tot[,col_annoation] != "Unannotated",] to have also the unannotated
  #cohen_tot[, col_annoation] <- factor(cohen_tot[, col_annoation], levels = rev(levels))
  colnames(cohen_tot)[4] <- c("moleculare_fam")
  cohen_tot$species_label <- paste0("italic('", cohen_tot$species, "')")
  # Plot effect sizes
  plot_final <- ggplot(cohen_tot, aes(x = ES, y = moleculare_fam, 
                                      fill = expression, 
                                      color = expression, 
                                      alpha = significance)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_bar(stat = "identity") +
    facet_wrap(
      vars(species_label),
      labeller = label_parsed
    ) +
    geom_errorbar(aes(xmin = CI_inf, xmax = CI_sup), width = 0, 
                  linewidth = 1.4, position = position_dodge(0.7)) +
    scale_color_manual(values = c(rep("black", 100))) +
    my_theme + theme(legend.position = "bottom") +
    labs(x = "Cohen's d", y = NULL) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    guides(color=guide_legend(nrow=2,byrow=TRUE)) +
    guides(alpha=guide_legend(nrow=2,byrow=TRUE))
  
  #print(plot_final)
  
  return(plot_final)
}


# run the function

annotations$feature_id <- as.character(annotations$feature_id)

# all plant species, control Vs H5
selected_plants <- c("L. repens", "L. grandiflora", 
                     "M. hippuroides", "P. palustris", "M. aquaticum")

metabo1 <- metabo %>%
  dplyr::filter(ATTTRIBUTE_species %in% selected_plants) %>%
  dplyr::filter(
    ATTTRIBUTE_species != "M. aquaticum" |
      Treatment %in% c("Control", "Herbivory 5 days")
  )

metabo1$Treatment <- factor(
  metabo1$Treatment,
  levels = c("Control", "Herbivory 5 days")
)
Figure4 <- df_annotation_treatment(
  metabo = metabo1,
  annotations = annotations,
  col_annoation = "pathways",
  col_feature_table = col_feature_table,
  col_metadata = col_metadata
);Figure4

metabo_2 <- metabo %>%
  filter(
    ATTTRIBUTE_species == "M. aquaticum",
    Treatment %in% c("Control", "Herbivory 1 day (plant material collected after 120h)")
  )
metabo_2$Treatment <- factor(
  metabo_2$Treatment,
  levels = c("Control", "Herbivory 1 day (plant material collected after 120h)")
)
metabo_PF_df  <- subset(metabo, ATTTRIBUTE_species == "M. aquaticum")
metabo_PF <- metabo_PF_df[, (col_metadata + 1):ncol(metabo_PF_df)]
annotations_PF <- annotations[
  annotations$feature_id %in% colnames(metabo_PF),
]
col_annoation="pathways"


Figure5 <- df_annotation_treatment(metabo_2, annotations_PF,
                                   col_annoation, col_feature_table, 
                                   col_metadata, levels); Figure5


#### Stats ####

get_cohen_stats <- function(metabo, annotations, col_feature_table, col_metadata) {
  
  cohen_tot <- list()
  
  for(i_sp in unique(metabo$ATTTRIBUTE_species)){
    
    metabo_i_sp <- subset(metabo, ATTTRIBUTE_species == i_sp)
    
    matrix_superclass <- metabo_i_sp[, col_feature_table:ncol(metabo_i_sp)]
    matrix_superclass <- matrix_superclass[, colSums(matrix_superclass) != 0]
    
    matrix_scaled <- as.data.frame(cbind(
      metabo_i_sp[,1:col_metadata],
      scale(matrix_superclass)
    ))
    
    df_scaled <- tidyr::gather(
      matrix_scaled,
      key = "feature_id",
      value = "Expression",
      all_of(colnames(matrix_superclass))
    )
    
    cohen_d_list <- list()
    
    for(el in unique(df_scaled$feature_id)){
      
      df_sub <- df_scaled %>% filter(feature_id == el)
      
      es <- effsize::cohen.d(
        df_sub$Expression,
        df_sub$Treatment,
        hedges.correction = TRUE,
        na.rm = TRUE
      )
      
      cohen_d_list[[el]] <- data.frame(
        feature_id = el,
        cohen_d = as.numeric(es$estimate),
        CI_inf = es$conf.int[1],
        CI_sup = es$conf.int[2]
      )
    }
    
    cohen_d_df <- bind_rows(cohen_d_list)
    
    annotation_2 <- annotations %>%
      filter(feature_id %in% colnames(matrix_superclass))
    
    annotation_2$cohen_d <- cohen_d_df$cohen_d
    annotation_2$var_d <- ((cohen_d_df$CI_sup - cohen_d_df$CI_inf)/(2*1.96))^2
    
    pooled_ES <- list()
    
    for(annot in unique(annotation_2$pathways)){
      
      df_annot <- annotation_2 %>% filter(pathways == annot)
      
      res <- metafor::rma(
        yi = df_annot$cohen_d,
        vi = df_annot$var_d,
        method = "REML"
      )
      
      pooled_ES[[annot]] <- data.frame(
        ES = as.numeric(res$b),
        CI_inf = res$ci.lb,
        CI_sup = res$ci.ub,
        p_value = res$pval,
        z_value = res$zval,
        species = i_sp,
        pathway = annot
      )
    }
    
    cohen_tot[[i_sp]] <- bind_rows(pooled_ES)
  }
  
  bind_rows(cohen_tot)
}


Table5 <- get_cohen_stats(
  metabo1,
  annotations,
  col_feature_table,
  col_metadata
)

Table5 <- Table5 %>%
  mutate(
    p_sig = ifelse(p_value < 0.05, "< 0.05", "≥ 0.05")
  )

Table5 <- Table5 %>%
  group_split(species) %>%
  purrr::map(~ {
    
    gt(.x) %>%
      
      cols_label(
        species = "Species",
        ES = "Effect size (ES)",
        CI_inf = "CI lower",
        CI_sup = "CI upper",
        p_sig = "p-value",
        z_value = "z-value",
        pathway = "Pathway"
      ) %>%
      
      tab_style(
        style = cell_text(style = "italic"),
        locations = cells_body(columns = species)
      ) %>%
      
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(
          columns = p_sig,
          rows = p_sig == "< 0.05"
        )
      ) %>%
      
      fmt_number(
        columns = c(ES, CI_inf, CI_sup, z_value),
        decimals = 3
      ) %>%
      
      cols_move(
        columns = p_sig,
        after = z_value
      ) %>%
      
      cols_hide(columns = p_value)
  })

Table5_LG = Table5[1]; Table5_LG
Table5_LR = Table5[2]; Table5_LR
Table5_PF = Table5[3]; Table5_PF
Table5_MH = Table5[4]; Table5_MH
Table5_PP = Table5[5]; Table5_PP
# Export 800 x 420

Table6 <- get_cohen_stats(
  metabo_2,
  annotations,
  col_feature_table,
  col_metadata
)

Table6 <- Table6 %>%
  mutate(
    p_sig = ifelse(p_value < 0.05, "< 0.05", "≥ 0.05")
  )

Table6 <- Table6 %>%
  gt() %>%
  
  cols_label(
    species = "Species",
    pathway = "Pathway",
    ES = "Effect size (ES)",
    CI_inf = "CI lower",
    CI_sup = "CI upper",
    z_value = "z-value",
    p_sig = "p-value"
  ) %>%
  
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = species)
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p_sig,
      rows = p_sig == "< 0.05"
    )
  ) %>%
  
  fmt_number(
    columns = c(ES, CI_inf, CI_sup, z_value),
    decimals = 3
  ) %>%
  
  cols_move(
    columns = p_sig,
    after = z_value
  ) %>%
  
  cols_hide(columns = p_value); Table6
# Export 800 x 420

# Print session information for reproducibility (R version + package versions)
sessionInfo()

#R version 4.4.2 (2024-10-31 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 10 x64 (build 19045)

#Matrix products: default


#locale:
#  [1] LC_COLLATE=French_Switzerland.utf8  LC_CTYPE=French_Switzerland.utf8    LC_MONETARY=French_Switzerland.utf8
#[4] LC_NUMERIC=C                        LC_TIME=French_Switzerland.utf8    

#time zone: Europe/Paris
#tzcode source: internal

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] gt_1.3.0            metafor_4.8-0       numDeriv_2016.8-1.1 metadat_1.4-0       effsize_0.8.1      
#[6] mixOmics_6.30.0     MASS_7.3-61         ggpubr_0.6.0        factoextra_1.0.7    vegan_2.7-2        
#[11] permute_0.9-8       lmerTest_3.1-3      lme4_1.1-35.5       Matrix_1.7-1        emmeans_1.10.5     
#[16] Rmisc_1.5.1         plyr_1.8.9          lattice_0.22-6      DHARMa_0.5.0        car_3.1-3          
#[21] carData_3.0-5       lubridate_1.9.3     forcats_1.0.0       stringr_1.6.0       dplyr_1.1.4        
#[26] purrr_1.2.2         readr_2.1.5         tidyr_1.3.2         tibble_3.2.1        ggplot2_4.0.2      
#[31] tidyverse_2.0.0    

#loaded via a namespace (and not attached):
#  [1] tidyselect_1.2.1     farver_2.1.2         S7_0.2.0             fastmap_1.2.0        TH.data_1.1-2       
#[6] mathjaxr_2.0-0       digest_0.6.39        timechange_0.3.0     estimability_1.5.1   lifecycle_1.0.5     
#[11] cluster_2.1.8.2      survival_3.7-0       magrittr_2.0.3       compiler_4.4.2       sass_0.4.9          
#[16] rlang_1.2.0          tools_4.4.2          igraph_2.1.4         ggsignif_0.6.4       labeling_0.4.3      
#[21] rARPACK_0.11-0       xml2_1.3.6           RColorBrewer_1.1-3   multcomp_1.4-26      abind_1.4-8         
#[26] BiocParallel_1.40.0  withr_3.0.2          grid_4.4.2           xtable_1.8-4         scales_1.4.0        
#[31] dichromat_2.0-0.1    cli_3.6.5            mvtnorm_1.3-2        ellipse_0.5.0        pairwiseAdonis_0.4.1
#[36] generics_0.1.4       rstudioapi_0.18.0    RSpectra_0.16-2      reshape2_1.4.5       tzdb_0.4.0          
#[41] minqa_1.2.8          splines_4.4.2        parallel_4.4.2       matrixStats_1.4.1    vctrs_0.6.5         
#[46] sandwich_3.1-1       boot_1.3-31          hms_1.1.3            rstatix_0.7.2        ggrepel_0.9.6       
#[51] Formula_1.2-5        glue_1.8.0           nloptr_2.1.1         codetools_0.2-20     stringi_1.8.7       
#[56] gtable_0.3.6         pillar_1.11.1        htmltools_0.5.8.1    R6_2.6.1             backports_1.5.1     
#[61] broom_1.0.7          corpcor_1.6.10       Rcpp_1.0.13-1        coda_0.19-4.1        gridExtra_2.3       
#[66] nlme_3.1-166         mgcv_1.9-4           fs_1.6.6             zoo_1.8-12           pkgconfig_2.0.3
