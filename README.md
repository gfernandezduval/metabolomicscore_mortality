# A multi-metabolite signature robustly predicts long-term mortality in the PREDIMED trial and several US cohorts

The R-code Scripts in this repository cover the main analysis:

-  *0_predimed_data.R* : Process baseline PREDIMED metabolomic data by handling missing values and applying transformations
-  *1_metabolome_wide_association.R* : Association of each individual metabolite with mortality
-  *2_metabolite_selection.R* : Feature selection of metabolites associated with mortality using Elastic Net Cox Regression
-  *3_baseline_predimed_score.R* : Calculate baseline PREDIMED Score using 10-fold LOFO (Leave-one-fold out) with Elastic Net Cox Regression
-  *4_score_mortality.R* : Association of baseline PREDIMED score with all-cause mortality
-  *5_validation_predimed1year.R* : Validation of 38-multimetabolite score in PREDIMED at 1 year of follow-up

If you have any questions, please contact Gonzalo Fernández-Duval ( ghfernandezd@unav.es )
