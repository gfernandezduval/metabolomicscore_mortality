# Author: Gonzalo Fernandez-Duval
# Title: Validation PREDIMED +1 year
# Description: Validation of 38-multimetabolite score in PREDIMED at 1 year of follow-up

# Libraries
library(dplyr)
library(stringr)
library(survival)

#=====================================================
#================== Preparing Data  ==================
#=====================================================

# Load required data
load("./input_data/predimed_1year.RData")  # 1-year follow-up data
load("./output_data/selected_metabolites.RData")  # Selected metabolites and coefficients

# Create 1-year metabolite names by adding "1"
metabolites_1year <- selected_metabolites %>%
  mutate(metabolites1 = str_c(metabolites, "1")) 

# Extract metabolite names and coefficients
list_metabolites1year <- metabolites_1year$metabolites1
coeff_metabolites <- as.numeric(metabolites_1year$coeff_original)

#=================================================
#=============== Score Calculation ===============
#=================================================

# Calculate metabolite score and final data
df_predimed_1year <- predimed1 %>%
  select(id, 
         all_of(list_metabolites1year),  # same metabolites, same order
         age0:death19) %>%
  mutate(score1 = as.numeric(as.matrix(select(., all_of(list_metabolites1year))) %*% coeff_metabolites),  # multiply coefficients with metabolite value
         score1_q = factor(ntile(score1, 5), levels = c("1", "2", "3", "4","5")),
         score1_stand = as.numeric(scale(score1)))

#=====================================================
#================== Cox Models ======================
#===================================================== 

# Define covariate sets
basic_covariates <- c("age10_1", "sex")

full_covariates <- c("age10_1", "sex", "diabetes0", "wth_cat_1", 
                     "smoking1", "alcohol_cat1", "hyperten0", 
                     "educ2", "fam_history", "dyslip0", "bmi_cat_1",
                     "energyt1_stand", "getotal1_stand", "p14_1", "group_int")

# Function to create formula string
create_formula <- function(score_var, covariates = NULL) {
  if (is.null(covariates)) {
    paste("Surv(follow_up19_1, death19) ~", score_var)
  } 
  else {
    paste("Surv(follow_up19_1, death19) ~", score_var, "+", 
          paste(covariates, collapse = " + "))
  }
}

# Function to fit and analyze Cox model
analyze_cox_model <- function(formula, data, model_name = "") {
  # Fit model
  model <- coxph(as.formula(formula), data = data)
  
  # Get summary
  model_summary <- summary(model)
  
  # Check proportional hazards assumption
  ph_test <- cox.zph(model)
  
  return(list(model = model,
              summary = model_summary,
              ph_test = ph_test))
}

# Define models with formulas
models_1year <- list(score1_quintiles = create_formula("score1_q"),
                     score1_standardized = create_formula("score1_stand"),
                     score1_quintiles_agesex = create_formula("score1_q", basic_covariates),
                     score1_stand_agesex = create_formula("score1_stand", basic_covariates),
                     score1_quintiles_full = create_formula("score1_q", full_covariates),
                     score1_stand_full = create_formula("score1_stand", full_covariates))

# Analyze all models
results_1year <- list()
for (model_name in names(models_1year)) {
  results_1year[[model_name]] <- analyze_cox_model(formula = models_1year[[model_name]],
                                                   data = df_predimed_1year,
                                                   model_name = model_name)
}

# Save data
save(df_predimed_1year, file = "./output_data/df_predimed_1year.RData")
