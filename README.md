# Multinomial Logistic Regression Analysis of Dermatological Diseases

## Overview

This project applies generalized linear modeling techniques to investigate the relationship between clinical and dermatological characteristics and the classification of dermatological diseases.

The analysis uses the Dermatology dataset from the UCI Machine Learning Repository, which contains clinical and histopathological attributes of patients diagnosed with different dermatological conditions.

The primary objective is to examine how selected clinical characteristics contribute to distinguishing among different dermatological disease categories using multinomial logistic regression.

## Research Objectives

- Examine the associations between selected clinical characteristics and dermatological disease classification.
- Evaluate the influence of patient age, erythema, scaling, definite borders, itching, and the Koebner phenomenon on disease classification.
- Investigate whether the association between itching and disease classification varies with age through an interaction term.
- Assess the performance and interpretability of the multinomial logistic regression model.

## Statistical Methods

- Exploratory data analysis (EDA)
- Multinomial logistic regression
- Generalized linear models (GLMs)
- Categorical predictor analysis
- Interaction modeling (age-by-itching)
- Model assessment and diagnostic evaluation
- Interpretation of regression coefficients and relative odds

## Tools and Technologies

- R
- R Markdown
- tidyverse
- nnet
- ggplot2
- dplyr

## Dataset

The analysis uses the Dermatology dataset from the UCI Machine Learning Repository.

The response variable consists of six dermatological disease categories:

- Psoriasis
- Seborrheic dermatitis
- Lichen planus
- Pityriasis rosea
- Chronic dermatitis
- Pityriasis rubra pilaris
  
#### The analysis focuses on selected predictors, including:

- Age
- Erythema
- Scaling
- Definite borders
- Itching
- Koebner phenomenon
- Observations with missing values in the age variable are excluded from the modeling analysis.

## Repository Structure

- Dermatology_Analysis.Rmd: Source file containing the statistical analysis and R code.
- Dermatology_Analysis.pdf: Knitted PDF report presenting the analysis, results, and interpretations.
- README.md: Project overview and documentation.

## Reproducibility

The R Markdown source file documents the analytical workflow, including data preparation, exploratory analysis, model specification, estimation, and interpretation.

Required R packages should be installed before running the analysis. The original dataset can be obtained from the UCI Machine Learning Repository.

## Data Source

- UCI Machine Learning Repository: Dermatology Dataset.

- https://archive.ics.uci.edu/dataset/33/dermatology

## Author

- Bernard Sefah
- MS Statistics
- Department of Mathematical Sciences
- Montana State University

