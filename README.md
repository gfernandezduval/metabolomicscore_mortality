A multi-metabolite signature robustly predicts long-term mortality in the PREDIMED trial and several US cohorts

The R-code Scripts in this repository cover the main analysis:

0) Baseline PREDIMED data: Process metabolomic data by handling missing values and applying transformations
1) Metabolome Wide Association: Association of each individual metabolite with mortality
2) Feature Selection of metabolites: Feature selection of metabolites associated with mortality using Elastic Net Cox Regression
3) Baseline PREDIMED Score: Calculate baseline PREDIMED Score using 10-fold LOFO (Leave-one-fold out) with Elastic Net Cox Regression
4) Score association with all-cause mortality: Association of baseline PREDIMED score and mortality
5) Validation example (PREDIMED 1-year of follow-up): Validation of 38-multimetabolite score in PREDIMED at 1 year of follow-up

For any questions, please contact Gonzalo Fernández-Duval ( ghfernandezd@unav.es )
