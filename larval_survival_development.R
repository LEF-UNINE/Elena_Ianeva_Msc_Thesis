####### Code for Lysathia cilliersae Survival and Development Analysis ####### 


# Code last update: 2026-06-22
# R version 4.4.2 (2024-10-31 ucrt)
# Details of the R session at the end of the script

# Author:
# Elena Ianeva

# Clean R memory
rm(list =ls()) 

#### Load Packages #### 

# Load required packages for data manipulation, statistics, and visualization

library(readxl) # import Excel files
library(survival) # provides functions for survival analysis
# (Kaplan-Meier curves and Cox proportional hazards models)
library(survminer)# visualize survival analysis results
library(gtsummary) # for summary and regression tables
library(blme) # for Bayesian mixed models
library(arm) # for regression and hierarchical models
library(lmerTest) # p-values for mixed models
library(lattice) # graphics
library(effects) # to visualize variable effects
library(emmeans) # post hoc comparisons and estimated marginal means
library(car) # ANOVA
library(DHARMa) # residual diagnostics
library(ggpubr) # plots based on ggplot2
library(patchwork) # combine multiple ggplot2 figures into a single layout
library(tibble) # data frames
library(dplyr) # data manipulation
library(ggplot2) # data visualization
library(tidyr) # reshape data
library(gt) # tables for reports

#### Bayesian generalized linear model #### 

setwd("C:/Users/Elena Ianeva/Desktop/master_thesis/data")
data <- read.csv(
  "data_larvae.csv",
  sep = ";",
  stringsAsFactors = TRUE
)
data$Batch = as.factor(data$Batch)
levels(data$Treatment)[2] <- "Herbivory 5 days"
data$Plant_species <- relevel(data$Plant_species, ref = "M. aquaticum")
summary(data)

# explore the data
xtabs(~ Treatment + Plant_species, data = data)
xtabs(~ Treatment + Plant_species + Batch, data = data) 
xtabs(~ Treatment + Plant_species + Survival, data = data)  
# because of these zeros (complete separation),
# glmer does not work
plot(party::ctree(factor(Survival) ~ Treatment + Plant_species , data = data))  
# right branch is M.aquaticum
# left branch: L. grandiflora, L. repens, M. hippuroides, P. palustris



model1 = bayesglm(Survival ~ Batch + Treatment * Plant_species,
                  family = binomial,     data = data)
summary(model1)
Anova(model1)

##### Table 1 ##### 

table1 <- Anova(model1, type = 2) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  mutate(
    p.display = ifelse(`Pr(>Chisq)` < 0.05, "<0.05", ">0.05")
  )

table1 %>%
  select(term, `LR Chisq`, Df, p.display) %>%
  gt() %>%
  cols_label(
    term = "Term",
    `LR Chisq` = "LR Chi²",
    Df = "Df",
    p.display = "p-value"
  ) %>%
  fmt_number(columns = `LR Chisq`, decimals = 2) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p.display,
      rows = table1$`Pr(>Chisq)` < 0.05
    )
  ); table1


xtabs(~ factor(predict(model1, type="response")>= 0.5) + data$Survival) 
# 20% errors

# see effects : 
plot(effect("Treatment:Plant_species", model1), x.var="Treatment")


# takes the emmeans and compares the levels of treatment within each plant species (link scale by default)
pairs(emmeans(model1, "Treatment", by="Plant_species"))  # only in M. hippuroides we have Treat effect!!!!!!!
# on the response scale (odds ratio)
pairs(emmeans(model1, "Treatment", by = "Plant_species"), type = "response")

# estimated marginal meand (EMMs) of treatment within each plant species, on the response scale
emmeans(model1, "Treatment", by="Plant_species", type = "response")


emm_pairs <- pairs(
  emmeans(model1, ~ Treatment | Plant_species),
  type = "response"
)

##### Table 2 ##### 
df_pairs <- as.data.frame(emm_pairs)
df_pairs_clean <- df_pairs %>%
  mutate(
    p_label = ifelse(p.value < 0.05, "< 0.05", "> 0.05"),
    odds.ratio = round(odds.ratio, 3),
    SE = round(SE, 3),
    z.ratio = round(z.ratio, 3)
  )


table2 <- df_pairs_clean %>%
  gt(groupname_col = "Plant_species") %>%
  
  cols_label(
    contrast = "Comparison",
    odds.ratio = "Odds ratio",
    SE = "SE",
    p_label = "p-value",
    z.ratio = "z.ratio"
  ) %>%
  
  tab_header(
    title = "Pairwise comparisons of Treatment",
    subtitle = "emmeans contrasts (response scale)"
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p_label,
      rows = df_pairs_clean$p.value < 0.05
    )
  ) %>%
  cols_hide(columns = c(p.value)) %>%
  
  tab_options(
    row_group.font.weight = "bold"
  ); table2



# model assumptions check:
#devtools::install_github(repo = "florianhartig/DHARMa", subdir = "DHARMa", dependencies = T)

res <- simulateResiduals(model1)

##### Figure S3 ####
par(mfrow = c(3, 1))  

plot(res)
plotResiduals(res, form = model.frame(model1)$Treatment)
plotResiduals(res, form = model.frame(model1)$Batch)
plotResiduals(res, form = model.frame(model1)$Plant_species)

par(mfrow = c(1, 1))  


##### Figure 1 ##### 
# Plot (Obs survival for each plant species and each treatment)

control_col <- "#A8D5BA"
herb_col    <- "#E8A6A6"

batch_data <- data %>%
  filter(Plant_species %in% c("M. aquaticum", "M. hippuroides")) %>%
  group_by(Plant_species, Treatment, Batch) %>%
  summarise(
    survival_batch = mean(Survival),
    .groups = "drop"
  )

batch_summary <- batch_data %>%
  group_by(Plant_species, Treatment) %>%
  summarise(
    survival = mean(survival_batch) * 100,
    se = sd(survival_batch) / sqrt(n()) * 100,
    .groups = "drop"
  )

other_summary <- data %>%
  filter(!Plant_species %in% c("M. aquaticum", "M. hippuroides")) %>%
  group_by(Plant_species, Treatment) %>%
  summarise(
    survival = mean(Survival) * 100,
    se = NA,
    .groups = "drop"
  )

surv_data <- bind_rows(batch_summary, other_summary)

surv_data$Plant_species <- factor(surv_data$Plant_species)


annot <- surv_data %>%
  group_by(Plant_species) %>%
  summarise(
    y = max(survival + coalesce(se, 0), na.rm = TRUE) + 5,
    label = ifelse(Plant_species == "M. hippuroides",
                   "p < 0.05",
                   "n.s"),
    .groups = "drop"
  )

dodge <- position_dodge(width = 0.75)

species_levels <- levels(surv_data$Plant_species)
annot$x <- match(annot$Plant_species, species_levels)

annot_lines <- data.frame(
  x = annot$x - 0.2,
  xend = annot$x + 0.2,
  y = annot$y,
  yend = annot$y
)

Figure1= ggplot(surv_data,
       aes(x = Plant_species,
           y = survival,
           fill = Treatment)) +
  
  # bars
  geom_col(position = dodge,
           width = 0.7) +
  
  # error bars
  geom_errorbar(
    aes(ymin = survival - se,
        ymax = survival + se),
    position = dodge,
    width = 0.2,
    na.rm = TRUE
  ) +
  
  # colors
  scale_fill_manual(values = c(
    "Control" = control_col,
    "Herbivory 5 days" = herb_col
  )) +
  
  # HORIZONTAL BRACKET
  geom_segment(
    data = annot_lines,
    aes(x = x,
        xend = xend,
        y = y,
        yend = y),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  
  # LABEL
  geom_text(
    data = annot,
    aes(x = Plant_species,
        y = y + 3,
        label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.15))) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "italic", hjust = 0.5),
  ) +
  
  labs(
    x = NULL,
    y = "Survival (%)",
    fill = "Treatment"
  );Figure1


#### Kaplan-Meier curves + Cox models #### 
##### Batch 1 #####


###### M. aquaticum ###### 

# data cox model 
data_B1_PF <- 
  read_excel("C:/Users/Elena Ianeva/Desktop/master_thesis/data/datasheet.xlsx", 
                         sheet = "cox_model", range = "A1:G25")
data_B1_PF$Treatment= as.factor(data_B1_PF$Treatment)
levels(data_B1_PF$Treatment) <- c("Control", "Herbivory\n5 days")

# Cox proportional hazards model 
cox_model_1 = coxph(Surv(time, status) ~ Treatment, data = data_B1_PF)

# compute kaplan meier survival curves (only descriptive, not a model test)
fit <- survfit(Surv(time, status) ~ Treatment, data = data_B1_PF)
summary(fit)

plot_PF_B1 <- ggsurvplot(
  fit,
  data = data_B1_PF,
  palette = c("#40AA5F", "#F15854"),
  ggtheme = theme_minimal(),
  #risk.table = TRUE,      # ajoute number at risk
  #pval = TRUE, 
  #pval.method = TRUE,
  pval.coord = c(15, 0.2),
  xlim = c(0, 25),
  linetype = c("solid", "dashed"),
  break.time.by = 5,
  ylab = "Survival probability",
  #legend = "right",
  legend.title = "Treatment",
  legend.labs = c("Control", "Herbivory 5 days"),
  
  title = expression(italic("M. aquaticum"))
); plot_PF_B1


######  P. palustris  ######  

# Data cox
data_B1_PP = read_excel("C:/Users/Elena Ianeva/Desktop/master_thesis/data/datasheet.xlsx", 
                        sheet = "cox_model", range = "A30:G54")
data_B1_PP$Treatment= as.factor(data_B1_PP$Treatment)
levels(data_B1_PP$Treatment)  <- c("Control", "Herbivory\n5 days")

cox_model_2 = coxph(Surv(time, status) ~ Treatment, data = data_B1_PP)

fit <- survfit(Surv(time, status) ~ Treatment, data = data_B1_PP)
summary(fit)

plot_PP_B1 <- ggsurvplot(
  fit,
  data = data_B1_PP,
  palette = c("#40AA5F", "#F15854"),
  linetype = c("solid", "dashed"),
  #pval = TRUE, 
  pval.coord = c(15, 0.25),
  ggtheme = theme_minimal(),
  xlim = c(0, 25),
  break.time.by = 5,
  ylab = "Survival probability",
  #legend = "right",
  legend.title = "Treatment",
  legend.labs = c("Control", "Herbivory 5 days"),
  
  title = expression(italic("P. palustris"))
); plot_PP_B1


######  M. hippuroides  ######  

# Data cox
data_B1_MH = read_excel("C:/Users/Elena Ianeva/Desktop/master_thesis/data/datasheet.xlsx", 
                        sheet = "cox_model", range = "A59:G83")
data_B1_MH$Treatment= as.factor(data_B1_MH$Treatment)
levels(data_B1_MH$Treatment) <- c("Control", "Herbivory\n5 days")

library(survival)


cox_model_3 = coxph(Surv(time, status) ~ Treatment, data = data_B1_MH)


fit <- survfit(Surv(time, status) ~ Treatment, data = data_B1_MH)
summary(fit)

plot_MH_B1 <- ggsurvplot(
  fit,
  data = data_B1_MH,
  palette = c("#40AA5F", "#F15854"), 
  linetype = c("solid", "dashed"),
  #pval = TRUE, 
  pval.coord = c(15, 0.2),
  ggtheme = theme_minimal(),
  xlim = c(0, 25),
  break.time.by = 5,
  ylab = "Survival probability",
  #legend = "right",
  legend.title = "Treatment",
  legend.labs = c("Control", "Herbivory 5 days"),
  
  title = expression(italic("M. hippuroides"))
); plot_MH_B1


##### Batch 2 #####

######  M. aquaticum  ######  

# Data cox
data_B2_PF = read_excel("C:/Users/Elena Ianeva/Desktop/master_thesis/data/datasheet.xlsx", 
                        sheet = "cox_model", range = "A85:G109")
data_B2_PF$Treatment= as.factor(data_B2_PF$Treatment)
levels(data_B2_PF$Treatment) <- c("Control", "Herbivory\n5 days")

library(survival)


cox_model_4 = coxph(Surv(time, status) ~ Treatment, data = data_B2_PF)

fit <- survfit(Surv(time, status) ~ Treatment, data = data_B2_PF)
summary(fit)

plot_PF_B2 <- ggsurvplot(
  fit,
  data = data_B2_PF,
  palette = c("#40AA5F", "#F15854"),
  linetype = c("solid", "dashed"),
  ggtheme = theme_minimal(),
  xlim = c(0, 30),
  #pval = TRUE, 
  pval.coord = c(20, 0.2),
  break.time.by =5,
  ylab = "Survival probability",
  #legend = "right",
  legend.title = "Treatment",
  legend.labs = c("Control", "Herbivory 5 days"),
  
  title = expression(italic("M. aquaticum"))
); plot_PF_B2



######  L. repens  ######  

# Data cox
data_B2_LR = read_excel("C:/Users/Elena Ianeva/Desktop/master_thesis/data/datasheet.xlsx", 
                        sheet = "cox_model", range = "A113:G137")
data_B2_LR$Treatment= as.factor(data_B2_LR$Treatment)
levels(data_B2_LR$Treatment) <- c("Control", "Herbivory\n5 days")

library(survival)


cox_model_5 = coxph(Surv(time, status) ~ Treatment, data = data_B2_LR)

fit <- survfit(Surv(time, status) ~ Treatment, data = data_B2_LR)
summary(fit)

plot_LR_B2 <- ggsurvplot(
  fit,
  data = data_B2_LR,
  palette = c("#40AA5F", "#F15854"),  linetype = c("solid", "dashed"),
  
  ggtheme = theme_minimal(),
  xlim = c(0, 45),
  break.time.by = 5,
  #pval = TRUE, 
  pval.coord = c(30, 0.1),
  ylab = "Survival probability",
  #legend = "right",
  legend.title = "Treatment",
  legend.labs = c("Control", "Herbivory 5 days"),
  
  title = expression(italic("L. repens"))
); plot_LR_B2


######  L. grandiflora  ######  

# Data cox
data_B2_LG = read_excel("C:/Users/Elena Ianeva/Desktop/master_thesis/data/datasheet.xlsx", 
                        sheet = "cox_model", range = "A142:G166")
data_B2_LG$Treatment= as.factor(data_B2_LG$Treatment)
levels(data_B2_LG$Treatment) <- c("Control", "Herbivory\n5 days")

library(survival)


cox_model_6 = coxph(Surv(time, status) ~ Treatment, data = data_B2_LG)

fit <- survfit(Surv(time, status) ~ Treatment, data = data_B2_LG)
summary(fit)

plot_LG_B2 <- ggsurvplot(
  fit,
  data = data_B2_LG,
  palette = c("#40AA5F", "#F15854"),
  linetype = c("solid", "dashed"),
  ggtheme = theme_minimal(),
  xlim = c(0, 35),
  #pval = TRUE, 
  pval.coord = c(10, 0.1),
  break.time.by = 5,
  ylab = "Survival probability",
  #legend = "right",
  legend.title = "Treatment",
  legend.labs = c("Control", "Herbivory 5 days"),
  
  title = expression(italic("L. grandiflora"))
); plot_LG_B2


##### Batch 3 #####

######  M. aquaticum  ###### 

# Data cox
data_B3_PF = read_excel("C:/Users/Elena Ianeva/Desktop/master_thesis/data/datasheet.xlsx", 
                        sheet = "cox_model", range = "A194:G218")
data_B3_PF$Treatment= as.factor(data_B3_PF$Treatment)
levels(data_B3_PF$Treatment) <- c("Control", "Herbivory\n5 days")

library(survival)


cox_model_7 = coxph(Surv(time, status) ~ Treatment, data = data_B3_PF)

fit <- survfit(Surv(time, status) ~ Treatment, data = data_B3_PF)
summary(fit)

plot_PF_B3 <- ggsurvplot(
  fit,
  data = data_B3_PF,
  palette = c("#40AA5F", "#F15854"),
  linetype = c("solid", "dashed"),
  
  ggtheme = theme_minimal(),
  xlim = c(0, 30),
  break.time.by = 5,
  #pval = TRUE, 
  pval.coord = c(20, 0.2),
  ylab = "Survival probability",
  #legend = "right",
  legend.title = "Treatment",
  legend.labs = c("Control", "Herbivory 5 days"),
  
  title = expression(italic("M. aquaticum"))
); plot_PF_B3



######  M. hippuroides  ###### 


# Data cox
data_B3_MH = read_excel("C:/Users/Elena Ianeva/Desktop/master_thesis/data/datasheet.xlsx", 
                        sheet = "cox_model", range = "A168:G190")
data_B3_MH$Treatment= as.factor(data_B3_MH$Treatment)
levels(data_B3_MH$Treatment) <- c("Control", "Herbivory\n5 days")

library(survival)


cox_model_8 = coxph(Surv(time, status) ~ Treatment, data = data_B3_MH)

fit <- survfit(Surv(time, status) ~ Treatment, data = data_B3_MH)
summary(fit)

plot_MH_B3 <- ggsurvplot(
  fit,
  data = data_B3_MH,
  linetype = c("solid", "dashed"),
  
  palette = c("#40AA5F", "#F15854"),
  ggtheme = theme_minimal(),
  xlim = c(0, 25),
  break.time.by = 5,
  #pval = TRUE, 
  pval.coord = c(15, 0.2),
  ylab = "Survival probability",
  #legend = "right",
  legend.title = "Treatment",
  legend.labs = c("Control", "Herbivory 5 days"),
  
  title = expression(italic("M. hippuroides"))
); plot_MH_B3


##### Figure 2 ##### 

style_surv <- function(p) {
  p$plot <- p$plot +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      plot.margin = margin(5, 5, 5, 5)
    )
  p$plot$layers[[length(p$plot$layers)]]$aes_params$size <- 8/3
  return(p)
}
Fig_2 = ggarrange(
  style_surv(plot_PF_B1)$plot,
  style_surv(plot_PP_B1)$plot,
  style_surv(plot_MH_B1)$plot,
  style_surv(plot_PF_B2)$plot,
  style_surv(plot_LR_B2)$plot,
  style_surv(plot_LG_B2)$plot,
  style_surv(plot_PF_B3)$plot,
  style_surv(plot_MH_B3)$plot, common.legend = T, labels = LETTERS[1:8]
) ; Fig_2 # 700 x 550 

##### Table 3 ##### 

tbl1 <- tbl_regression(cox_model_1, exponentiate = TRUE)%>% bold_p()
tbl2 <- tbl_regression(cox_model_2, exponentiate = TRUE)%>% bold_p()
tbl3 <- tbl_regression(cox_model_3, exponentiate = TRUE)%>% bold_p()
tbl4 <- tbl_regression(cox_model_4, exponentiate = TRUE)%>% bold_p()
tbl5 <- tbl_regression(cox_model_5, exponentiate = TRUE)%>% bold_p()
tbl6 <- tbl_regression(cox_model_6, exponentiate = TRUE)%>% bold_p()
tbl7 <- tbl_regression(cox_model_7, exponentiate = TRUE)%>% bold_p()
tbl8 <- tbl_regression(cox_model_8, exponentiate = TRUE)%>% bold_p()

Table3 <- tbl_merge(
  tbls = list(tbl1, tbl2, tbl3, tbl4, tbl5, tbl6, tbl7, tbl8),
  tab_spanner = c(
    "M. aquaticum",
    "P. palustris",
    "M. hippuroides",
    "M. aquaticum",
    "L. repens",
    "L. grandiflora",
    "M. aquaticum",
    "M. hippuroides"
  )
); Table3



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
#  [1] tidyr_1.3.2         tibble_3.2.1        gt_1.3.0            gtsummary_2.5.0    
#[5] dplyr_1.1.4         multcompView_0.1-10 patchwork_1.3.2     DHARMa_0.5.0       
#[9] car_3.1-3           emmeans_1.10.5      effects_4.2-2       carData_3.0-5      
#[13] lattice_0.22-6      lmerTest_3.1-3      arm_1.14-4          MASS_7.3-61        
#[17] blme_1.0-7          lme4_1.1-35.5       Matrix_1.7-1        survminer_0.5.1    
#[21] ggpubr_0.6.0        ggplot2_4.0.2       survival_3.7-0      readxl_1.4.3       

#loaded via a namespace (and not attached):
#  [1] Rdpack_2.6.2         DBI_1.2.3            gridExtra_2.3        sandwich_3.1-1      
#[5] rlang_1.2.0          magrittr_2.0.3       multcomp_1.4-26      matrixStats_1.4.1   
#[9] compiler_4.4.2       mgcv_1.9-4           vctrs_0.6.5          stringr_1.6.0       
#[13] pkgconfig_2.0.3      fastmap_1.2.0        backports_1.5.1      labeling_0.4.3      
#[17] KMsurv_0.1-6         promises_1.3.2       markdown_1.13        haven_2.5.4         
#[21] nloptr_2.1.1         purrr_1.2.2          xfun_0.49            modeltools_0.2-23   
#[25] labelled_2.16.0      later_1.4.1          broom_1.0.7          parallel_4.4.2      
#[29] gap.datasets_0.0.6   R6_2.6.1             qgam_1.3.4           coin_1.4-3          
#[33] stringi_1.8.7        RColorBrewer_1.1-3   boot_1.3-31          cellranger_1.1.0    
#[37] numDeriv_2016.8-1.1  estimability_1.5.1   iterators_1.0.14     Rcpp_1.0.13-1       
#[41] knitr_1.49           zoo_1.8-12           httpuv_1.6.17        splines_4.4.2       
#[45] nnet_7.3-19          tidyselect_1.2.1     rstudioapi_0.18.0    dichromat_2.0-0.1   
#[49] abind_1.4-8          doParallel_1.0.17    codetools_0.2-20     plyr_1.8.9          
#[53] shiny_1.9.1          withr_3.0.2          S7_0.2.0             coda_0.19-4.1       
#[57] evaluate_1.0.1       survey_4.4-2         xml2_1.3.6           survMisc_0.5.6      
#[61] pillar_1.11.1        party_1.3-17         gap_1.6              foreach_1.5.2       
#[65] stats4_4.4.2         insight_1.5.1        generics_0.1.4       hms_1.1.3           
#[69] commonmark_1.9.2     scales_1.4.0         minqa_1.2.8          xtable_1.8-4        
#[73] glue_1.8.0           tools_4.4.2          data.table_1.16.2    ggsignif_0.6.4      
#[77] forcats_1.0.0        fs_1.6.6             mvtnorm_1.3-2        cowplot_1.1.3       
#[81] grid_4.4.2           mitools_2.4          rbibutils_2.3        libcoin_1.0-10      
#[85] cards_0.7.1          colorspace_2.1-1     nlme_3.1-166         Formula_1.2-5       
#[89] cli_3.6.5            km.ci_0.5-6          broom.helpers_1.22.0 rematch_2.0.0       
#[93] strucchange_1.5-4    gtable_0.3.6         rstatix_0.7.2        sass_0.4.9          
#[97] digest_0.6.39        TH.data_1.1-2        farver_2.1.2         htmltools_0.5.8.1   
#[101] lifecycle_1.0.5      mime_0.12 

