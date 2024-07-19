# title: Full Data Simulation 
# author: Jinyu Luo
# version: 2024-07-18
rm(list = ls())

# Required R packages ----------------------------------------------------------
library(tidyverse)
library(geepack)
library(parallel)
library(survival)
library(simsurv)
library(doParallel) 
source("qls_functions.R")
source("helperFunctions.R")

# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less core than available
cl <- makeCluster(num_cores)
registerDoParallel(cl)

load("Outputs/matched_data.RData")
# Formulas for model fit 
formula_red <- y ~ bav*visit 
formula_full <- y ~ bav*visit + age + male + bsa

# Model Specification 
model_specs <- list(
  ind_mdl_full = list(formula = formula_full, corstr = "independence", adjusted = TRUE),
  ind_mdl_red = list(formula = formula_red, corstr = "independence", adjusted = FALSE),
  ar1_mdl_full = list(formula = formula_full, corstr = "ar1", adjusted = TRUE),
  ar1_mdl_red = list(formula = formula_red, corstr = "ar1", adjusted = FALSE),
  exch_mdl_full = list(formula = formula_full, corstr = "exchangeable", adjusted = TRUE),
  exch_mdl_red = list(formula = formula_red, corstr = "exchangeable", adjusted = FALSE)
)

gee_fit_full = fit_gee(matched_df$matched_full)
gee_fit_dropout = fit_gee(matched_df$matched_ddat)

print("GEE Model fitting completed.")

# QLS --------------------------------------------------------------------------
# Function to transform matched data to QLS data 
to_QLS_data <- function(df) {
  df %>%
    arrange(matched) %>%
    mutate(clusterID = as.integer(factor(matched))) %>%
    select(-matched) %>%
    group_by(clusterID) %>%
    mutate(cluster.var = ifelse(bav == 0, 1, 2),
           order = row_number()) %>%
    relocate(c(clusterID, cluster.var, order), .after = id)
}

get_qls_results <- function(df, formula, corstr, adjusted) {
  
  if(adjusted){
    print(paste("Fitting with adjusted", corstr, "correlation structure"))
  }
  
  print(paste("Fitting with unadjusted", corstr, "correlation structure"))
  
  model <- tryCatch({
    qls(formula, data = df, time.var = df$visit, maxT = 6, corstr = corstr)
  }, error = function(e) {
    print(paste("Error fitting model:", e$message))
    NULL
  })
  
  if (is.null(model)) {
    print("Model fitting failed.")
    return(data.frame(term = NA, estimate = NA, std_error = NA, 
                      convergence = FALSE, rho = NA, tau = NA))
  }
  
  print("Model fitting succeeded.")
  
  lower <- model$coefficients - qnorm(1-(1-0.95)/2) * model$se
  upper <- model$coefficients + qnorm(1-(1-0.95)/2) * model$se
  N <- nrow(df)
  p <- length(model$coefficients) - 1
  adj_se <- sqrt(diag((N/(N-p))*model$vcov))
  adj_lower <- model$coefficients - qnorm(1-(1-0.95)/2) * adj_se
  adj_upper <- model$coefficients + qnorm(1-(1-0.95)/2) * adj_se
  
  result <- data.frame(term = names(model$se), 
                       estimate = model$coefficients, 
                       std_error = model$se, 
                       adj_std_error = adj_se, 
                       lower = lower, upper = upper, 
                       adj_lower = adj_lower, adj_upper = adj_upper,
                       convergence = TRUE, 
                       rho = model$alpha,
                       tau = model$tau) %>% cbind(model$vcov)
  
  return(result)
}

# Combine the transformation and fitting in a single foreach loop
combined_results <- foreach(i = seq_len(nrow(matched_df)), 
                            .packages = c('tidyverse', 'geepack', 'parallel', 'survival', 'simsurv'),
                            .export = c('to_QLS_data', 'get_qls_results', 'model_specs')) %dopar% {
                              
                              # Transformation
                              qls_full_data <- to_QLS_data(matched_df$matched_full[[i]])
                              qls_drop_data <- to_QLS_data(matched_df$matched_ddat[[i]])
                              
                              n_pairs_full <- n_distinct(qls_full_data$clusterID)
                              n_pairs_drop <- n_distinct(qls_drop_data$clusterID)
                              
                              # Fitting models for qls_full_data
                              qls_full_model_results <- list()
                              for (mdl in names(model_specs)) {
                                spec <- model_specs[[mdl]]
                                qls_full_model_results[[mdl]] <- get_qls_results(qls_full_data, spec$formula, spec$corstr, spec$adjusted)
                              }
                              
                              # Fitting models for qls_drop_data
                              qls_drop_model_results <- list()
                              for (mdl in names(model_specs)) {
                                spec <- model_specs[[mdl]]
                                qls_drop_model_results[[mdl]] <- get_qls_results(qls_drop_data, spec$formula, spec$corstr, spec$adjusted)
                              }
                              
                              list(
                                sim_id = i,
                                n_pairs_full = n_pairs_full,
                                n_pairs_drop = n_pairs_drop,
                                qls_full_model_results = qls_full_model_results,
                                qls_drop_model_results = qls_drop_model_results
                              )
                            }

# Combine results into tibbles
temp_qls <- tibble(
  sim_id = map(combined_results, "sim_id"),
  n_pairs_full = map_int(combined_results, "n_pairs_full"),
  n_pairs_drop = map_int(combined_results, "n_pairs_drop")
)

qls_full_out <- tibble(
  sim_id = map(combined_results, "sim_id"),
  ind_mdl_full = map(combined_results, ~ .x$qls_full_model_results$ind_mdl_full), 
  ind_mdl_red = map(combined_results, ~ .x$qls_full_model_results$ind_mdl_red), 
  ar1_mdl_full = map(combined_results, ~ .x$qls_full_model_results$ar1_mdl_full), 
  ar1_mdl_red = map(combined_results, ~ .x$qls_full_model_results$ar1_mdl_red), 
  exch_mdl_full = map(combined_results, ~ .x$qls_full_model_results$exch_mdl_full), 
  exch_mdl_red = map(combined_results, ~ .x$qls_full_model_results$exch_mdl_red)
)

qls_drop_out <- tibble(
  sim_id = map(combined_results, "sim_id"),
  ind_mdl_full = map(combined_results, ~ .x$qls_drop_model_results$ind_mdl_full), 
  ind_mdl_red = map(combined_results, ~ .x$qls_drop_model_results$ind_mdl_red), 
  ar1_mdl_full = map(combined_results, ~ .x$qls_drop_model_results$ar1_mdl_full), 
  ar1_mdl_red = map(combined_results, ~ .x$qls_drop_model_results$ar1_mdl_red), 
  exch_mdl_full = map(combined_results, ~ .x$qls_drop_model_results$exch_mdl_full), 
  exch_mdl_red = map(combined_results, ~ .x$qls_drop_model_results$exch_mdl_red)
)

# Stop the parallel backend
stopCluster(cl)

print("QLS Model fitting completed.")

save(temp_qls, 
     gee_fit_full, gee_fit_dropout, 
     qls_full_out, qls_drop_out,
     file = "Outputs/simfit.RData")

