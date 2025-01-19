# Author: Gonzalo Fernandez-Duval
# Title: Metabolite selection
# Description: Feature selection of metabolites associated with mortality using Elastic Net Cox Regression

# Libraries
library(dplyr)
library(survival)
library(glmnet)

# Data
#source("0_predimed_data.R")
load("./output_data/df_mortality.RData")

# Define metabolite column range
first_m <- 2
last_m <- 338

# Only metabolites for modeling
x_metabolites <- df_mortality %>%
  select(first_m:last_m) %>%              # Select metabolite columns
  as.matrix()

#==========================================================
#============= Metabolites Feature Selection  =============
#==========================================================

# Set seed
set.seed(123)

# Cross-validated Elastic Net Cox Regression
cv_model <- cv.glmnet(x = x_metabolites, 
                      y = Surv(df_mortality$follow_up19, df_mortality$death19), 
                      alpha = 0.5, 
                      nfolds = 10,
                      family = "cox")

# Extract coefficients using 1 standard error lambda
model_coef <- coef(cv_model, s= cv_model$lambda.1se)

# Create data frame of selected metabolites and their coefficients
selected_metabolites <- data.frame(
  metabolites = colnames(x_metabolites)[model_coef@i + 1],
  coeff_original = round(model_coef@x, digits = 4))

# Save selected metabolites
save(selected_metabolites, file = "./output_data/selected_metabolites.RData")
