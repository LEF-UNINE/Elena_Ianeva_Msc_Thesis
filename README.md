# Elena_Ianeva_Msc_Thesis_2026
This is a repository for the code and the data used for my Master's Thesis

## Datasets
There are ... datasets

    - formatted_peak_table.csv
    - MetaData_Elena_26.csv
    - data_larvae.csv
    - cox_model.csv
    
## Description of the different datasets    
### 1. formatted_peak_table.csv
   
This table contains metabolite features and their corresonding peak intensities 

### 2. MetaData_Elena_26.csv
   
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

### 3. data_larvae.csv

This dataset includes the following variables :
 - name of the sample
 - plant species (*M. aquaticum* / *M. hippuroides* / *L. repens* / *L. grandiflora* / *P. palustris*)
 - survival (binary : 0 (dead) OR 1 (alive))
 - treatment (Control / Induced). NB: Induced is equivalent to Herbivory 5 days and is recoded as such in the R script.
 - Batch (1/2/3)

### 4. cox_model.csv

This dataset includes the following variables : 
 - sample id
 - time : duration (in days) for which the larva was observed alive. If the larva reached the pupal stage, "time" corresponds to the time required to reach that stage
 - status : event indicator (binary) --> 0 (the larva reached the pupal stage = censored observation), 1 (the larva died before reaching the pupal stage)
 - plant species (*M. aquaticum* / *M. hippuroides* / *L. repens* / *L. grandiflora* / *P. palustris*)
 - treatment (Control / Hebivory 5 days)
 - Survival to pupa (binary : 0 (yes) 1 (no))


## R scripts    
There are ... scripts

    - chem_richness_shannon_diversity_PLSDA.R 
    - larval_survival_development.R

## Description of the R scripts 
### 1. chem_richness_shannon_diversity_PLSDA.R

Aim: Analyze and visualize chemical richness, Shannon diversity, and partial least squares discriminant analysis (PLS-DA) to compare control and 5-day herbivory treatments within each plant species

Outputs generated: 
 - Fig 3 : Chemical Richness, Shannon diversity and PLS-DA comparisons between control and 5-day herbivory treatments within the following plant species : *M. aquaticum*, *M. hippuroides*, *L. repens*, *L. grandiflora* and *P. palustris*
 - Table S3: Statistical results for chemical richness
 - Table S4: Statistical results for Shannon diversity

Files used in this script: 

    - formatted_peak_table.csv
    - MetaData_Elena_26.csv

### 2. larval_survival_development.R

Aim: Analyse and visualize Lysathia cilliersae survival and development 

Outputs generated: 
 - Fig 1 : Observed survival (%) on each plant species (*M. aquaticum*, *M. hippuroides*, *L. repens*, *L. grandiflora* and *P. palustris*) under Control and 5-day herbivory treatments 
 - Table 1: Analysis of fixed effects from the Bayesian generalized model used to analyze the observed survival on each plant species under the two conditions
 - Table 2: Pairwise comparisons of Treatment effects within each plant species based on estimated marginal means (emmeans) from the Bayesian generalized linear model
 - Figure S3: Diagnostic plots of scaled residuals from the Bayesian generalized linear model 
 - Figure 2: Kaplan-Meier curves (for each plant species and each batch)
 - Table 3: Statistical results of the cox model
   
Files used in this script: 

    - data_larvae.csv
    - cox_model.csv





    
  


   
