# title: Fit QLS on Simulated Data 
# author: Jinyu Luo
# version: 2024-06-30

rm(list = ls())
# Required R packages ----------------------------------------------------------
library(tidyverse)
library(geepack)
library(parallel)
library(survival)
library(simsurv)
library(doParallel) 
load("data/realcoefs.RData")
load("Outputs/matched_data.RData")
source("qls_functions.R")

# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less core than available
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Global Variables -------------------------------------------------------------
maxT <- 6
alpha_ci <- 0.05
n_sim <- 1500

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

# GEE --------------------------------------------------------------------------
# Get the matched data
# matched_data <- matched_df %>% pull(matched_full)
matched_data <- matched_df %>% pull(matched_ddat)

# Function to fit GEE with different correlation structures on simulated data
get_gee_results <- function(df, formula, corstr, adjusted) {
  
  if(adjusted){
    print(paste("Fitting with adjusted", corstr, "correlation structure"))
  }
  
  print(paste("Fitting with unadjusted", corstr, "correlation structure"))
  
  model <- tryCatch({
      geeglm(formula, family = binomial('logit'), wave = factor(visit),
             corstr = corstr, id = id, data = df)
  }, error = function(e) {
    print(paste("Error fitting model:", e$message))
    NULL
  })
  
  if (is.null(model)) {
    print("Model fitting failed.")
    return(data.frame(term = NA, estimate = NA, std_error = NA, lower = NA, upper = NA, 
                      adj_lower = NA, adj_upper = NA, convergence = FALSE))
  }
  
  print("Model fitting succeeded.")
  
  fit <- summary(model)
  est <- fit$coefficients
  z <- qnorm(1 - alpha_ci / 2)
  lower <- est[, "Estimate"] - z * est[, "Std.err"]
  upper <- est[, "Estimate"] + z * est[, "Std.err"]
  
  p <- nrow(est)-1
  N <- nrow(df)
  v_cov <- vcov(model)
  V_df <- (N / (N - p)) * v_cov
  adj_se <- sqrt(diag(V_df))
  adj_lower <- est[, "Estimate"] - z * adj_se
  adj_upper <- est[, "Estimate"] + z * adj_se
  
  rho <- ifelse(corstr == "independence", 0, fit$corr[1,1])
  rho_se <- ifelse(corstr == "independence", 0, fit$corr[1,2])
  
  result <- data.frame(term = rownames(est), 
                       estimate = est[, "Estimate"], 
                       std_error = est[, "Std.err"], 
                       adj_std_error = adj_se,
                       lower = lower, upper = upper, 
                       adj_lower = adj_lower, adj_upper = adj_upper,
                       convergence =TRUE, 
                       rho = rho, 
                       rho_se = rho_se) 
  return(result)
}

# Initialize a list to store the results
gee_fits <- list()

# Loop through each dataset and fit all models
for (i in 1:length(matched_data)) {
  df <- matched_data[[i]]
  sim_id <- i
  model_results <- list()
  
  for (mdl in names(model_specs)) {
    m <- model_specs[[mdl]]
    model_results[[mdl]] <- get_gee_results(df = df, formula = m$formula, corstr = m$corstr, adjusted = m$adjusted)
  }
  
  gee_fits[[i]] <- c(list(sim_id = sim_id), model_results)
  
  print(paste("Completed simulation", i, "out of", n_sim))
}

# Convert the results list to a dataframe
gee_fits_df <- tibble(
  sim_id = map(gee_fits, "sim_id"),
  ind_mdl_full = map(gee_fits, "ind_mdl_full"),
  ind_mdl_red = map(gee_fits, "ind_mdl_red"), 
  ar1_mdl_full = map(gee_fits, "ar1_mdl_full"), 
  ar1_mdl_red = map(gee_fits, "ar1_mdl_red"), 
  exch_mdl_full = map(gee_fits, "exch_mdl_full"), 
  exch_mdl_red = map(gee_fits, "exch_mdl_red")
)

print("GEE Model fitting completed.")

# QLS --------------------------------------------------------------------------
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

qls_df <- matched_df %>% select(matched_ddat) %>%  
  mutate(sim_id = 1:n_sim, 
         matched_data = map(matched_ddat, ~ .x %>% 
                              arrange(matched) %>% 
                              mutate(clusterID = as.integer(factor(matched))) %>% 
                              select(-matched) %>% 
                              group_by(clusterID) %>% 
                              mutate(cluster.var = ifelse(bav == 0, 1, 2), 
                                     order = row_number()) %>% 
                              relocate(c(clusterID, cluster.var, order), .after = id)), 
         n_pairs = map_int(matched_data, ~ n_distinct(.x$clusterID)))


# Parallel processing using foreach
qls_results <- foreach(i = seq_len(nrow(qls_df)), .packages = c('tidyverse', 'geepack', 'parallel', 'survival', 'simsurv')) %dopar% {
  df <- qls_df$matched_data[[i]]
  sim_id <- qls_df$sim_id[i]
  
  model_results <- list()
  for (mdl in names(model_specs)) {
    spec <- model_specs[[mdl]]
    model_results[[mdl]] <- get_qls_results(df, spec$formula, spec$corstr, spec$adjusted)
  }
  
  c(list(sim_id = sim_id), model_results)
}

# Stop the parallel backend
stopCluster(cl)

# Convert the results to a tibble
qls_fits_df <- tibble(
  sim_id = map(qls_results, "sim_id"),
  ind_mdl_full = map(qls_results, "ind_mdl_full"), 
  ind_mdl_red = map(qls_results, "ind_mdl_red"), 
  ar1_mdl_full = map(qls_results, "ar1_mdl_full"), 
  ar1_mdl_red = map(qls_results, "ar1_mdl_red"), 
  exch_mdl_full = map(qls_results, "exch_mdl_full"), 
  exch_mdl_red = map(qls_results, "exch_mdl_red")
)

print("QLS Model fitting completed.")

save(gee_fits_df, qls_fits_df, file = "Outputs/simfit.RData")


