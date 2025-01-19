# Author: Gonzalo Fernandez-Duval
# Title: ETL Metabolomic Data
# Description: Process metabolomic data by handling missing values and applying transformations

# Libraries
library(dplyr)
library(rcompanion)

# Load metabolomic data
load("./input_data/predimed_metabolites.RData")

# Load survival variables and covariables
load("./input_data/predimed_variables.RData")

#========================================
#============== Remove NA's =============
#========================================

# 1. Remove metabolites with more than 20% NA's (columns)
na_count <- data.frame(name = names(predimed_metabolites)[-1],
                       n_null = sapply(predimed_metabolites[-1], function(x) sum(is.na(x)))) %>%
  mutate(perc = (n_null/nrow(predimed_metabolites)) * 100,
         if_twenty = ifelse(perc >= 20, 1, 0)) %>%
  arrange(desc(perc))

# Get list of metabolites to remove
list_colnull <- na_count %>%
  filter(if_twenty == 1) %>%
  pull(name)

# 2. Remove participants with more than 20% NA's (rows) 
metabolites_idnull <- predimed_metabolites %>%
  mutate(count_na = rowSums(is.na(.)),
         perc_na = (count_na/(ncol(predimed_metabolites) - 1)) * 100,
         if_twenty_na = ifelse(perc_na >= 20, 1, 0)) %>%
  filter(if_twenty_na == 1)

# Get list of participants to remove
list_rownull <- metabolites_idnull$id

# Remove both high-NA metabolites and participants
metabolites_clean <- predimed_metabolites %>%
  select(-all_of(list_colnull)) %>%         # columns
  filter(!(id %in% list_rownull))           # rows

#=========================================
#============== Impute NA's ==============
#=========================================

# Impute remaining NA's with half of minimum value for each metabolite
metabolites_imputed <- metabolites_clean

for(col in names(metabolites_imputed)[-1]) {  # Skip id column
  if(any(is.na(metabolites_imputed[[col]]))) {
    min_value <- min(metabolites_imputed[[col]], na.rm = TRUE) / 2
    metabolites_imputed[[col]][is.na(metabolites_imputed[[col]])] <- min_value
  }
}

# Apply Blom's transformation to each metabolite
df_metabolites <- metabolites_imputed %>%
  mutate(across(-id, ~blom(., method = "blom"))) %>%
  as.data.frame()

#===================================================
#============== Merge all needed data ==============
#===================================================

# Merge all data: id, metabolites, survival variables and covariables
df_mortality <- df_metabolites %>%
  inner_join(predimed_variables, by = "id") 

save(df_mortality, file = "./output_data/df_mortality.RData")
