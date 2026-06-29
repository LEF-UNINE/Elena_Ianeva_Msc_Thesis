####### Code for Random Forest Analysis on Myriophyllum hippuroides Samples ####### 


# Code last update: 2026-06-23
# R version 4.4.2 (2024-10-31 ucrt)
# Details of the R session at the end of the script

# Author:
# Elena Ianeva

# Clean R memory
rm(list =ls()) 

#### Load Packages #### 

# Load required packages for data manipulation, statistics, and visualization
library(readr) # import csv and text files
library(tidyverse) # data manipulation
library(dplyr) # To modify the data
library(ggplot2) # for visualization 
library(randomForest) # Random Forest modelling
library(car) # stat tests and ANOVA functions
library(gt) # to create tables in HTML
library(hrbrthemes) # ggplot themes
library(viridis) # color palettes for plots
library(patchwork) # to combine multiple ggplot objects





#### Load data ####

# import annotations
setwd("C:/Users/Elena Ianeva/Desktop/Treatment")
formatted_annotation <- read.csv("formatted_annotation.csv", sep = ";")   
formatted_annotation <- as.data.frame(formatted_annotation)


# import features
setwd("C:/Users/Elena Ianeva/Desktop/Treatment")
features <- read.csv("formatted_peak_table.csv", sep=",")
features$Treatment <- as.character(features$Treatment)
features$Treatment[features$Treatment == "Induced"] <- "Herbivory 5 days"
# keep only samples from the Control and the 5-day herbivory treatments
features_filtered <- features %>%
  filter(Treatment %in% c("Control", "Herbivory 5 days"))


# import metadata
setwd("C:/Users/Elena Ianeva/Desktop/master_thesis/metabolo_26")

MetaData <- read.csv("MetaData_Elena_26_with_survival.csv", sep=";")
MetaData$X.1 <- NULL
MetaData$X <- NULL
MetaData$Treatment[MetaData$Treatment == "Induced"] <- "Herbivory 5 days"
MetaData <- MetaData %>%
  filter(Treatment %in% c("Control", "Herbivory 5 days"))
MetaData <- 
  MetaData[ !(MetaData$Plant_species %in% 
                c("M. verticillatum", "M. heterophyllum")), ]
MetaData$correct_name <- sub(".*sample_|\\.mzML$", "", MetaData$correct_name)
MetaData <- MetaData[MetaData$correct_name != "",]
MetaData$correct_name <- sub("\\.mzML$", "", MetaData$correct_name)


# merge both 
df <- merge(MetaData, features_filtered,
            by.x = "correct_name", by.y = "sample_id")


df$Treatment.x <- as.factor(df$Treatment.x)
df$Plant_species <- as.factor(df$Plant_species)
df$Survival <- as.factor(df$Survival)

species_name <- "M. hippuroides" # PUT HERE THE PLANT SPECIES I WANT 
df <- df %>%
  filter(Plant_species %in% species_name)

df_new <- as.data.frame(df[,-c(1:10)])
df_new <- as.data.frame(df_new[,-c(2:10)])
df_new$Survival= as.factor(df_new$Survival)


#### Random Forest ####

df_new[-1] <- lapply(df_new[-1], function(x) as.numeric(as.character(x)))

# Keep only features detected in at least one sample
df_new <- df_new[, c(TRUE, colSums(df_new[-1], na.rm = TRUE) != 0)]

set.seed(1)
rf_model_survival <- randomForest(
  Survival ~ ., 
  data = df_new,
  ntree = 1000,
  importance = TRUE
)

varImpPlot(
  rf_model_survival,
  type = NULL,
  n.var = 30,
  scale = TRUE,
  cex = 0.8,
  main = species_name
)

# only Mean Decrease Accuracy
varImpPlot(
  rf_model_survival,
  type = 1,   
  n.var = 30,
  scale = TRUE,
  cex = 0.8,
  main = species_name
)

# same plot but using ggplot
importance_df <- as.data.frame(importance(rf_model_survival))
importance_df$Feature <- rownames(importance_df)
importance_df <- importance_df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  slice(1:30)
importance_df$Feature_clean <- gsub("feature([0-9]+).*" , "  Feature \\1", 
                                    importance_df$Feature)
importance_df$Feature_clean <- trimws(importance_df$Feature_clean)


FigS7 = ggplot(importance_df, aes(x = reorder(Feature_clean, 
                                              MeanDecreaseAccuracy),
                              y = MeanDecreaseAccuracy)) +
  geom_col(fill = "lightsteelblue") +
  coord_flip() +  
  labs(
    # title = species_name,
    x = " ",
    y = "Mean Decrease Accuracy"
  ) +
  theme_minimal(base_size = 12)



#### Select the 30 most important features from the Random Forest model #### 
imp <- importance(rf_model_survival, type = 1)  # Mean decrease in accuracy
imp_df <- data.frame(
  feature = rownames(imp),
  importance = imp[,1]
)

top30 <- imp_df[order(imp_df$importance, decreasing = TRUE), ][1:30, "feature"]


#### create dataset containing survival, treatment and top features #### 
df_model <- df %>%
  dplyr::select(Survival, Treatment.x, dplyr::all_of(top30)) %>%
  tidyr::drop_na()

final_results <- data.frame()

#### test each feature individually #### 
for (f in top30) {
  
  if (!f %in% colnames(df_model)) next # skip if feature is missing
  
  df_tmp <- df_model[, c("Survival", "Treatment.x", f)]
  df_tmp <- df_tmp[complete.cases(df_tmp), ]
  
  # skip features with insufficient variation
  if (length(unique(df_tmp[[f]])) < 2) next
  if (length(unique(df_tmp$Treatment.x)) < 2) next
  if (length(unique(df_tmp$Survival)) < 2) next
  
  #### 1. TEST whether feature abundance differs between treatments
  
  sh <- shapiro.test(df_tmp[[f]])
  
  if (sh$p.value > 0.05) {
    test <- t.test(df_tmp[[f]] ~ df_tmp$Treatment.x)
    p_abund <- test$p.value
  } else {
    test <- wilcox.test(df_tmp[[f]] ~ df_tmp$Treatment.x)
    p_abund <- test$p.value
  }
  
  means <- tapply(df_tmp[[f]], df_tmp$Treatment.x, mean)
  
  if (any(is.na(means))) next
  
  higher_in_herbivory <- means["Herbivory 5 days"] > means["Control"]
  
  #### 2. Fit logistic regression model for survival 
  # test if Treatment and feature predict survival 
  model <- try(glm(Survival ~ ., data = df_tmp, family = binomial), 
               silent = TRUE) 
  if (class(model)[1] == "try-error") next
  
  summ <- summary(model)
  print(summ)
  if (nrow(summ$coefficients) < 2) next
  
  coef <- summ$coefficients[3,1] 
  # if the coef is < 0 --> the more the abundance is high 
  #the more the survival probability is low
  print(coef)
  
  #### 3. Test significance of metabolite effect on survival
  
  anova_res <- try(car::Anova(model, type = "II"), silent = TRUE)
  if (class(anova_res)[1] == "try-error") next
  print(anova_res)
  
  p_survival <- anova_res$`Pr(>Chisq)`[2]
  print(p_survival)
  
  #### 4. Keep metabolites meeting aall selection criteria
  
  if (p_abund < 0.05 & # features that have sign diff abundances
      higher_in_herbivory & # higher abundance under Herbivory
      p_survival < 0.05 & # sign effect on survival
      coef < 0) { # negative association with survival 
    
    final_results <- rbind(final_results, data.frame(
      feature = f,
      p_abundance = p_abund,
      p_survival = p_survival,
      coef = coef,
      mean_control = means["Control"],
      mean_herbivory = means["Herbivory 5 days"]
    ))
  }
}

#### Selected features #### 
final_results
final_features = final_results$feature

final_features <- c(
  "feature2351_2.772_659.2702",
  "feature1344_2.195_714.2861"
)

##### Fit indiv logistic regression models for each selected feature ##### 

run_glm_features <- function(df, final_features) {
  

  results <- data.frame()
  
  for (f in final_features) {
    
    # skip if feature not present 
    if (!f %in% colnames(df)) next
    
    df_tmp <- df[, c("Survival", f)]
    df_tmp <- df_tmp[complete.cases(df_tmp), ]
    
    # skip if not enough variation
    if (length(unique(df_tmp[[f]])) < 2) next
    if (length(unique(df_tmp$Survival)) < 2) next
    # fit logistic regression model 
    form <- as.formula(paste("Survival ~", f))
    
    model <- try(glm(form, data = df_tmp, family = binomial), silent = TRUE)
    if (inherits(model, "try-error")) next
    
    summ <- summary(model)
    
    coef_table <- summ$coefficients
    
    
    if (nrow(coef_table) < 2) next
    # extract feature effect estimates
    results <- rbind(results, data.frame(
      Feature = f,
      Estimate = coef_table[2, "Estimate"],
      Std.Error = coef_table[2, "Std. Error"],
      Z.value = coef_table[2, "z value"],
      p.value = coef_table[2, "Pr(>|z|)"]
    ))
  }
  
  return(results)
}
# run the function and extract feature effects
glm_results <- run_glm_features(df_model, final_features)

Table5 = gt(glm_results) %>%
  
  fmt_scientific(
    columns = c(Estimate, Std.Error, Z.value),
    decimals = 3
  ) %>%
  
  fmt_number(
    columns = p.value,
    decimals = 2
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p.value,
      rows = p.value < 0.05
    )
  ); Table5

test1 = glm(Survival ~ feature2351_2.772_659.2702 + Treatment.x, 
            data = df, family = binomial) # ok
test2 = glm(Survival ~ feature1344_2.195_714.2861 + Treatment.x, 
            data = df, family = binomial) # ok

Anova(test1)
Anova(test2)


# Check abundances of the two features in the two treatments
df_long <- df %>%
  dplyr::select(Treatment.x, dplyr::all_of(final_features)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(final_features),
    names_to = "feature",
    values_to = "abundance"
  )
df_long <- df_long %>%
  dplyr::rename(Treatment = Treatment.x)

df_long <- df_long %>%
  dplyr::mutate(
    feature = sub("feature([0-9]+)_.*", "Feature \\1", feature)
  )

df_summary <- df_long %>%
  dplyr::group_by(Treatment, feature) %>%
  dplyr::summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    se = sd(abundance, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

abundance = ggplot(df_long, aes(x = Treatment, y = abundance, 
                                fill = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ feature, scales = "free_y") +
  labs(x = NULL, y = "Abundance")+
  theme_bw()+
  scale_fill_brewer(palette = "Pastel1"); abundance


# check if these features are present in the other plant species

plot_feature <- function(feature_name) {
  
  feature_id <- sub("feature([0-9]+).*", "\\1", feature_name)
  name <- paste0("Feature ", feature_id)
  
  ggplot(features, aes_string(
    x = "ATTTRIBUTE_species",
    y = feature_name,
    fill = "ATTTRIBUTE_species"
  )) +
    geom_boxplot() +
    scale_fill_viridis_d(alpha = 0.6) +
    geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, 
                                 face = "italic")
    ) +
    ggtitle(name) +
    xlab("")+
    ylab("Abundance")
}


p1 = plot_feature("feature1344_2.195_714.2861") 
p5= plot_feature("feature2351_2.772_659.2702") 

Fig6 = (p1 | p5) ; Fig6 # 1129 x 400


# Print session information for reproducibility (R version + package versions)
#sessionInfo()
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
#  [1] patchwork_1.3.2      viridis_0.6.5        viridisLite_0.4.3    hrbrthemes_0.8.7     randomForest_4.7-1.2
#[6] gt_1.3.0             metafor_4.8-0        numDeriv_2016.8-1.1  metadat_1.4-0        effsize_0.8.1       
#[11] mixOmics_6.30.0      MASS_7.3-61          ggpubr_0.6.0         factoextra_1.0.7     vegan_2.7-2         
#[16] permute_0.9-8        lmerTest_3.1-3       lme4_1.1-35.5        Matrix_1.7-1         emmeans_1.10.5      
#[21] Rmisc_1.5.1          plyr_1.8.9           lattice_0.22-6       DHARMa_0.5.0         car_3.1-3           
#[26] carData_3.0-5        lubridate_1.9.3      forcats_1.0.0        stringr_1.6.0        dplyr_1.1.4         
#[31] purrr_1.2.2          readr_2.1.5          tidyr_1.3.2          tibble_3.2.1         ggplot2_4.0.2       
#[36] tidyverse_2.0.0     

#loaded via a namespace (and not attached):
#  [1] gridExtra_2.3           sandwich_3.1-1          rlang_1.2.0             magrittr_2.0.3         
#[5] multcomp_1.4-26         matrixStats_1.4.1       compiler_4.4.2          mgcv_1.9-4             
#[9] systemfonts_1.3.2       vctrs_0.6.5             reshape2_1.4.5          pkgconfig_2.0.3        
#[13] fastmap_1.2.0           backports_1.5.1         labeling_0.4.3          tzdb_0.4.0             
#[17] nloptr_2.1.1            BiocParallel_1.40.0     broom_1.0.7             parallel_4.4.2         
#[21] cluster_2.1.8.2         R6_2.6.1                stringi_1.8.7           RColorBrewer_1.1-3     
#[25] extrafontdb_1.0         boot_1.3-31             estimability_1.5.1      Rcpp_1.0.13-1          
#[29] zoo_1.8-12              extrafont_0.19          pairwiseAdonis_0.4.1    splines_4.4.2          
#[33] igraph_2.1.4            timechange_0.3.0        tidyselect_1.2.1        rstudioapi_0.18.0      
#[37] dichromat_2.0-0.1       abind_1.4-8             codetools_0.2-20        withr_3.0.2            
#[41] rARPACK_0.11-0          S7_0.2.0                coda_0.19-4.1           survival_3.7-0         
#[45] xml2_1.3.6              pillar_1.11.1           ellipse_0.5.0           generics_0.1.4         
#[49] mathjaxr_2.0-0          hms_1.1.3               scales_1.4.0            minqa_1.2.8            
#[53] xtable_1.8-4            glue_1.8.0              gdtools_0.5.0           tools_4.4.2            
#[57] RSpectra_0.16-2         ggsignif_0.6.4          fs_1.6.6                mvtnorm_1.3-2          
#[61] grid_4.4.2              Rttf2pt1_1.3.12         nlme_3.1-166            Formula_1.2-5          
#[65] cli_3.6.5               fontBitstreamVera_0.1.1 corpcor_1.6.10          gtable_0.3.6           
#[69] rstatix_0.7.2           fontquiver_0.2.1        sass_0.4.9              digest_0.6.39          
#[73] ggrepel_0.9.6           TH.data_1.1-2           farver_2.1.2            htmltools_0.5.8.1      
#[77] lifecycle_1.0.5         fontLiberation_0.1.0  

