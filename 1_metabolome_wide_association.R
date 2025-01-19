# Author: Gonzalo Fernandez-Duval
# Title: Metabolome-wide association
# Description: Association of each individual metabolite with mortality

# Libraries
library(dplyr)
library(survival)
library(readxl)

#=====================================================
#================== Preparing Data  ==================
#=====================================================

# Data
#source("0_predimed_data.R")
load("./output_data/df_mortality.RData")

# Define metabolite column range
first_m <- 2
last_m <- 338

# Extract metabolite names and create base dataframe
metabolites_names <- names(df_mortality[, first_m:last_m])

namesformal <- read_excel("input_data/metabolitenames.xlsx") 

metabonames <- data.frame(metabolites = metabolites_names) %>%
  left_join(namesformal, by = "metabolites") 

# Function to initialize result dataframes
initialize_df <- function(base_df, is_sex_specific = FALSE) {
  if (!is_sex_specific) {
    base_df %>% mutate(hr = NA,
                        lower_ci = NA,
                        upper_ci = NA,
                        pvalue = NA,
                        pvalue_adjusted = NA)
  } 
  else {
    base_df %>% mutate(hr_male = NA,
                        lower_ci_male = NA,
                        upper_ci_male = NA,
                        pvalue_male = NA,
                        pvalue_adjusted_male = NA,
                        hr_female = NA,
                        lower_ci_female = NA,
                        upper_ci_female = NA,
                        pvalue_female = NA,
                        pvalue_adjusted_female = NA)
  }
}

# Initial data without information
data_allcause <- initialize_df(metabonames)
data_cancer <- initialize_df(metabonames)
data_cvd <- initialize_df(metabonames)
data_other <- initialize_df(metabonames)
data_sex <- initialize_df(metabonames, is_sex_specific = TRUE)

#========================================================================
#================== Cox Regression to each metabolite  ==================
#========================================================================

# Function to fit each outcome Cox model
fit_cox_model <- function(metabolite, outcome, data) {
  coxph(formula = reformulate(termlabels = metabolite,
                              response = paste0('Surv(follow_up19, ', outcome, ')')),
        data = data)
}

# Function to extract Cox results
extract_cox_results <- function(model) {
  coeff <- summary(model)$coefficients
  conf_int <- exp(confint(model))
  c(hr = coeff[2],
    lower_ci = conf_int[1],
    upper_ci = conf_int[2],
    pvalue = coeff[5])
}

# Iteration for each metabolite
for (i in 1:nrow(data_allcause)) {
  
  # Console output iteration
  cat("Metabolite ", data_allcause[i,1], '\n')
  
  #Specific metabolite
  metabolite <- data_allcause[i, 1]
  
  # All-cause mortality
  cox_all <- fit_cox_model(metabolite, "death19", df_mortality)
  data_allcause[i, 3:6] <- extract_cox_results(cox_all)
  
  # Cancer mortality
  cox_cancer <- fit_cox_model(metabolite, "death_canc19", df_mortality)
  data_cancer[i, 3:6] <- extract_cox_results(cox_cancer)
  
  # CVD mortality
  cox_cvd <- fit_cox_model(metabolite, "death_cvd19", df_mortality)
  data_cvd[i, 3:6] <- extract_cox_results(cox_cvd)
  
  # Other causes mortality
  cox_other <- fit_cox_model(metabolite, "death_oth19", df_mortality)
  data_other[i, 3:6] <- extract_cox_results(cox_other)
  
  # Sex-stratified all-cause mortality
  cox_sex <- lapply(split(df_mortality, df_mortality$sex),
                    function(x) fit_cox_model(metabolite, "death19", x))
  
  # Extract results for males (1) and females (2)
  male_results <- extract_cox_results(cox_sex[[1]])
  female_results <- extract_cox_results(cox_sex[[2]])
  
  data_sex[i, 3:6] <- male_results
  data_sex[i, 8:11] <- female_results
  
}

# P value adjusted
data_allcause$pvalue_adjusted <- p.adjust(data_allcause$pvalue, method = "BH", n = nrow(data_allcause))
data_cancer$pvalue_adjusted <- p.adjust(data_cancer$pvalue, method = "BH", n = nrow(data_cancer))
data_cvd$pvalue_adjusted <- p.adjust(data_cvd$pvalue, method = "BH", n = nrow(data_cvd))
data_other$pvalue_adjusted <- p.adjust(data_other$pvalue, method = "BH", n = nrow(data_other))
data_sex$pvalue_adjusted_male <- p.adjust(data_sex$pvalue_male, method = "BH", n = nrow(data_sex))
data_sex$pvalue_adjusted_female <- p.adjust(data_sex$pvalue_female, method = "BH", n = nrow(data_sex))

