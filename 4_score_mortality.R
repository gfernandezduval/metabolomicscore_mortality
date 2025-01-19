# Author: Gonzalo Fernandez-Duval
# Title: Association between score and mortality
# Description: Association of baseline PREDIMED score and mortality, adjusted by covariables

# Libraries
library(dplyr)
library(survival)

#=====================================================
#================== Preparing Data  ==================
#=====================================================

# Data
#source("0_predimed_data.R")
load("./output_data/df_mortality.RData")

# Define metabolite column range
first_m <- 2
last_m <- 338

# Covariables names
covariables <- names(df_mortality[, (last_m + 1):ncol(df_mortality)])

#Loading the calculated score
load(file = ".\\output_data\\df_score.Rdata")
data_score <- select(df_score, id, score, score_q, score_stand)

#Loading selected metabolites
load(file = ".\\output_data\\selected_metabolites.Rdata")
selected_metabolites <- c(selected_metabolites$metabolites)

# Creates final baseline data for analysis
df_mortality_selection <- df_mortality %>%
  select(id,
         all_of(selected_metabolites),
         all_of(covariables)) %>%
  left_join(data_score, by = "id")

#=====================================================
#================== Cox Models ======================
#===================================================== 

# Define covariate sets
basic_covariates <- c("age10", "sex")

full_covariates <- c("age10", "sex", "diabetes0", "wth_cat", 
                     "smoking0", "alcohol_cat", "hyperten0", 
                     "educ2", "fam_history", "dyslip0", "bmi_cat",
                     "energyt_stand", "getotal_stand", "p14", "group_int")

# Function to create formula string
create_formula <- function(score_var, covariates = NULL) {
  if (is.null(covariates)) {
    paste("Surv(follow_up19, death19) ~", score_var)
  } 
  else {
    paste("Surv(follow_up19, death19) ~", score_var, "+", 
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
models <- list(score_quintiles = create_formula("score_q"),
               score_standardized = create_formula("score_stand"),
               score_quintiles_agesex = create_formula("score_q", basic_covariates),
               score_stand_agesex = create_formula("score_stand", basic_covariates),
               score_quintiles_full = create_formula("score_q", full_covariates),
               score_stand_full = create_formula("score_stand", full_covariates))

# Analyze all models
results <- list()
for (model_name in names(models)) {
  results[[model_name]] <- analyze_cox_model(
    formula = models[[model_name]],
    data = df_mortality_selection,
    model_name = model_name)
}

# Save data
save(df_mortality_selection, file = "./output_data/df_mortality_selection.RData")
