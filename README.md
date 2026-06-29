# Elena_Ianeva_Msc_Thesis_2026
This is a repository for the code and the data used for my Master's Thesis

## Datasets
There are ... datasets

    - formatted_peak_table.csv
    - formatted_annotation.csv
    - MetaData_Elena_26.csv
    - MetaData_Elena_26_with_survival.csv
    - data_larvae.csv
    - cox_model.xlsx
    - feature.csv
    - MetaData - Copie.csv
    
    
## Description of the different datasets  

### 1. formatted_peak_table.csv
   
This table contains metabolite features and their corresonding peak intensities 

### 2. formatted_annotation.csv

This table contains metabolite identification information.

### 3. MetaData_Elena_26.csv
   
This dataset includes the following variables :
 - the filename of the sample
 - type of the sample (qc/blank/sample)
 - run date (date of the chemical analysis)
 - correct name of the sample (as there was an issue with the labelling)
 - start date of the experiment (29.08.2025 / 03.10.2025 / 07.11.2025/ 15.01.2026) 
 - plant species (*M. aquaticum* / *M. verticillatum* / *M. heterophyllum* / *M. hippuroides* / *L. repens* / *L. grandiflora* / *P. palustris*)
 - exposure time (72h / 120h)
 - treatment (Control / Herbivory 1 day / Herbivory 5 days / Mechanical Damage)
 - wet weight of the sample (in mg)
 - quantity used for the extraction (in mg)

### 4. MetaData_Elena_26_with_survival.csv

This dataset is identical to MetaData_Elena_26.csv, with the addition of a binary Survival variable for *M. hippuroides* samples: 0 = larva did not reach the pupal stage, 1 = larva reached the pupal stage

### 5. data_larvae.csv

This dataset includes the following variables :
 - name of the sample
 - plant species (*M. aquaticum* / *M. hippuroides* / *L. repens* / *L. grandiflora* / *P. palustris*)
 - survival (binary : 0 (dead) OR 1 (alive))
 - treatment (Control / Induced). NB: Induced is equivalent to Herbivory 5 days and is recoded as such in the R script.
 - Batch (1/2/3)

### 6. cox_model.xlsx

This dataset includes the following variables : 
 - sample id
 - time : duration (in days) for which the larva was observed alive. If the larva reached the pupal stage, "time" corresponds to the time required to reach that stage
 - status : event indicator (binary) --> 0 (the larva reached the pupal stage = censored observation), 1 (the larva died before reaching the pupal stage)
 - plant species (*M. aquaticum* / *M. hippuroides* / *L. repens* / *L. grandiflora* / *P. palustris*)
 - treatment (Control / Hebivory 5 days)
 - Survival to pupa (binary : 0 (yes) 1 (no))


## R scripts    
There are 5 scripts

    - chem_richness_shannon_diversity_PLSDA.R 
    - larval_survival_development.R
    - random_forest_analysis.R
    - pathways.R
    - supp_mat.R

## Description of the R scripts 
### 1. chem_richness_shannon_diversity_PLSDA.R

Aim: Analyze and visualize chemical richness, Shannon diversity, and partial least squares discriminant analysis (PLS-DA) to compare control and 5-day herbivory treatments within each plant species

Outputs generated: 
 - Figure 3 : Chemical Richness, Shannon diversity and PLS-DA comparisons between control and 5-day herbivory treatments within the following plant species : *M. aquaticum*, *M. hippuroides*, *L. repens*, *L. grandiflora* and *P. palustris*
 - Table S3: Statistical results for chemical richness
 - Table S4: Statistical results for Shannon diversity

Files used in this script: 

    - formatted_peak_table.csv
    - MetaData_Elena_26.csv

### 2. larval_survival_development.R

Aim: Analyse and visualize Lysathia cilliersae survival and development 

Outputs generated: 
 - Figure 1 : Observed survival (%) on each plant species (*M. aquaticum*, *M. hippuroides*, *L. repens*, *L. grandiflora* and *P. palustris*) under Control and 5-day herbivory treatments 
 - Table 1: Analysis of fixed effects from the Bayesian generalized model used to analyze the observed survival on each plant species under the two conditions
 - Table 2: Pairwise comparisons of Treatment effects within each plant species based on estimated marginal means (emmeans) from the Bayesian generalized linear model
 - Figure S3: Diagnostic plots of scaled residuals from the Bayesian generalized linear model 
 - Figure 2: Kaplan-Meier curves (for each plant species and each batch)
 - Table 3: Statistical results of the cox model
   
Files used in this script: 

    - data_larvae.csv
    - cox_model.xlsx

### 3. random_forest_analysis.R

Aim: Random forest analysis to identify features associated with Lysathia cilliersae survival on *Myriophyllum hippuroides*

Outputs generated: 
- Figure 6 : Abundance of the two features retained after the filtering across the five tested plant species
- Table 5: Coefficients from generalized linear models testing the effect of each selected feature on survival
- Figure S7: The 30 most important features identified by the random forest model, ranked according to their contribution to classification performance 

Files used in this script: 

    - formatted_annotation.csv
    - formatted_peak_table.csv
    - MetaData_Elena_26_with_survival.csv

### 4. pathways.R

Aim: Overview of pathway regulation across plant species and treatments 

Outputs generated: 
- Figure 4 : Species-specific effect sizes of herbivory-induced changes in metabolomic pathways after 5 days of exposure to *L. cilliersae*
- Figure 5: Effect sizes of herbivory-induced changes in metabolomic pathways after 1 day of exposure to *L. cilliersae* in *M. aquaticum*
- Table S5: Pathway-level pooled effect sizes (Hedges’g) for herbivory-treated (5 days) versus control samples across plant species
- Table S6: Pathway-level pooled effect sizes (Hedges’g) for herbivory-treated (1 day) versus control samples in *M. aquaticum*

Files used in this script: 

    - formatted_annotation.csv
    - formatted_peak_table.csv

### 5. supp_mat.R

Aim: Generate supplementary Figures and Tables

Outputs generated: 
- Figure S4 : Partial Least Squares Discriminant Analysis (PLS-DA) score plot showing discrimination between control samples (0 h) and samples subjected to different exposure durations to *L. cilliersae*
- Figure S5: Partial Least Squares Discriminant Analysis (PLS-DA) score plot showing discrimination between the 5-day herbivory treatment and the mechanical damage treatments
- Figure S6: Non-metric Multidimensional Scaling (NMDS) score plot showing discrimination of metabolomic profiles according to plant species and treatment conditions
- Table S2: Pairwise PERMANOVA comparisons between Exposure Time groups

Files used in this script:

    - formatted_annotation.csv
    - formatted_peak_table.csv
    - feature.csv
    - MetaData.csv
    






    
  


   
