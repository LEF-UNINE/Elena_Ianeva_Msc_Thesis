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

