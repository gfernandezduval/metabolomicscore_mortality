# Author: Gonzalo Fernandez-Duval
# Title: Baseline PREDIMED Score 
# Description: Calculate baseline PREDIMED Score using 10-fold LOFO (Leave-one-fold out) with Elastic Net Cox Regression

# Libraries
library(dplyr)
library(survival)
library(glmnet)
library(cvTools)

#=====================================================
#================== Preparing Data  ==================
#=====================================================

# Data
#source("0_predimed_data.R")
load("./output_data/df_mortality.RData")

# Define metabolite column range
first_m <- 2
last_m <- 338

# Prepare modeling data
df_metabolomic <- df_mortality %>%
  select(id,
         first_m:last_m,         # metabolites
         follow_up19, 
         death19)

# Create matrix for modeling
x_metabolites <- as.matrix(df_metabolomic[,first_m:last_m])
metabolites_names <- colnames(x_metabolites)

# Data for coefficients to each participant
df_coefficients <- select(df_metabolomic, -follow_up19, -death19)
df_coefficients[,first_m:last_m] <- NA      # To fill after each iteration

#Data to add the metabolite coefficients in each fold
fold_metabolites <- data.frame(metabolite = metabolites_names,
                               matrix(NA, 
                                      nrow = length(metabolites_names), 
                                      ncol = 10, 
                                      dimnames = list(NULL, paste0("fold", 1:10))))

#======================================================
#================= 10-iteration LOFO  =================
#======================================================

# Set seed
set.seed(321)

# Create 10 folds
k = 10
n = nrow(df_mortality)
folds = cvFolds(n, K = k)

# Seed for each fold
seeds <- c(58, 8, 3259, 94, 321, 5, 89, 6829, 99, 231)

# Iterate through folds
for (i in 1:k) {
  
  cat("Iteration Fold ", i, '\n')
  
  #Seed
  set.seed(seeds[i])
  
  # Train and test indices
  train_idx <- folds$subsets[folds$which != i]
  test_idx <- folds$subsets[folds$which == i]
  
  # Finding the best lambda
  cv_fit <- cv.glmnet(x = x_metabolites[train_idx, ], 
                      y = Surv(df_metabolomic$follow_up19[train_idx], df_metabolomic$death19[train_idx]), 
                      alpha = 0.5, 
                      nfolds = 10,
                      family = "cox") 
  
  # Get coefficients
  model_coef <- coef(cv_fit, s= cv_fit$lambda.1se)
  selected_idx <- model_coef@i + 1
  coef_values <- model_coef@x      #round(, digits = 4)
  
  # Store selected metabolites and their coefficients for this fold
  fold_metabolites[selected_idx, paste0("fold", i)] <- coef_values
  
  # Save coefficients for test set participants
  test_ids <- df_metabolomic$id[test_idx]
  for(test_id in test_ids) {
    id_row <- which(df_coefficients$id == test_id)
    df_coefficients[id_row, selected_idx + 1] <- coef_values
  }
  
}

#============================================================
#================= Baseline PREDIMED Score  =================
#============================================================

# Data frame with individual multiplications for each metabolite
df_multiplied <- df_coefficients
df_multiplied[, first_m:last_m] <- df_coefficients[, first_m:last_m] * x_metabolites

# Final baseline Score
df_score <- df_multiplied %>%
  mutate(score = rowSums(select(., all_of(first_m):all_of(last_m)), na.rm = T),    # Apply coeff to each metabolite
         score_q = ntile(score, 5),
         score_q = factor(score_q, levels = c("1", "2", "3", "4", "5")),
         score_stand = as.numeric(scale(score)),
         score = round(score, digits = 4))

# Final table of each fold 
df_fold_metabolites <- fold_metabolites %>%
  mutate(count = rowSums(!is.na(across(fold1:fold10)), na.rm = T),
         coeff_mean = round(rowMeans(across(fold1:fold10), na.rm = T), digits = 4)) %>%
  arrange(desc(count))

#Save results
save(df_score, file = "./output_data/df_score.RData") 
save(df_fold_metabolites, file = "./output_data/df_fold_metabolites.RData")  
