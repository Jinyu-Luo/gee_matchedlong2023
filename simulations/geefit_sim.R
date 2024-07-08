# title: GEE Fit 
# author: Jinyu Luo
# version: 2024-06-30
rm(list = ls())
# Required R packages ----------------------------------------------------------
library(tidyverse)
library(geepack)
library(lme4)
library(parallel)
library(survival)
library(simsurv)
library(R.utils)
load("data/realcoefs.RData")
load("Outputs/sim_data.RData")
# Global Variables -------------------------------------------------------------
maxT <- 6
rho <- 0.3
alpha_ci <- 0.05

true_coefs <- full_ar1$coefficients$Estimate
names(true_coefs) <- c("g0", "bav", "visit", "age", "male", "bsa", "bav_visit")
corr_alpha <- full_ar1$corr$Estimate
surv_coefs <- surv_coefs[,-4]
rownames(surv_coefs)[4] <- "male"

# Get the matched data
matched_df <- sim_df %>% pull(matched_data)

# Helper Functions -------------------------------------------------------------
# 1. GEE fit Function 
get_gee_results <- function(df, formula, corstr, adjusted) {
  
  if(adjusted){
    print(paste("Fitting with adjusted", corstr, "correlation structure"))
  }
  
  print(paste("Fitting with unadjusted", corstr, "correlation structure"))
  
  model <- tryCatch({
    withTimeout({
      geeglm(formula, family = binomial('logit'), wave = factor(visit),
             corstr = corstr, id = id, data = df)
    }, timeout = 60, onTimeout = "warning")
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
  
  p <- nrow(est)
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
                       lower = lower, upper = upper, 
                       adj_lower = adj_lower, adj_upper = adj_upper,
                       convergence =TRUE, 
                       rho = rho, 
                       rho_se = rho_se) 
  return(result)
}

# Define the model specifications
model_specs <- list(
  ind_mdl_full = list(formula = y ~ visit * bav + age + male + bsa, 
                      corstr = "independence", adjusted = TRUE),
  ind_mdl_red = list(formula = y ~ visit * bav, 
                     corstr = "independence", adjusted = FALSE),
  ar1_mdl_full = list(formula = y ~ visit * bav + age + male + bsa, 
                      corstr = "ar1", adjusted = TRUE),
  ar1_mdl_red = list(formula = y ~ visit * bav, 
                     corstr = "ar1", adjusted = FALSE),
  exch_mdl_full = list(formula = y ~ visit * bav + age + male + bsa, 
                       corstr = "exchangeable", adjusted = TRUE),
  exch_mdl_red = list(formula = y ~ visit * bav, 
                      corstr = "exchangeable", adjusted = FALSE)
)


# GEE Fit ----------------------------------------------------------------------
# Initialize a list to store the results
gee_fits <- list()

# Loop through each dataset and fit all models
for (i in 1:length(matched_df)) {
  df <- matched_df[[i]]
  sim_id <- i
  model_results <- list()
  
  for (mdl in names(model_specs)) {
    m <- model_specs[[mdl]]
    model_results[[mdl]] <- get_gee_results(df, m$formula, m$corstr, m$adjusted)
  }
  
  gee_fits[[i]] <- c(list(sim_id = sim_id), model_results)
  
  print(paste("Completed simulation", i, "out of", nrow(sim_df)))
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

# gee_df <- sim_df %>% select(sim_id, matched_data) %>% 
#   mutate(
#     ind_mdl_full = map(matched_data, function(df) {
#       get_gee_results(df, y ~ visit * bav + age + male + bsa, "independence")
#     }),
#     ind_mdl_red = map(matched_data, function(df) {
#       get_gee_results(df, y ~ visit * bav, "independence")
#     }),
#     ar1_mdl_full = map(matched_data, function(df) {
#       get_gee_results(df, y ~ visit * bav + age + male + bsa, "ar1")
#     }),
#     ar1_mdl_red = map(matched_data, function(df) {
#       get_gee_results(df, y ~ visit * bav, "ar1")
#     }),
#     exch_mdl_full = map(matched_data, function(df) {
#       get_gee_results(df, y ~ visit * bav + age + male + bsa, "exchangeable")
#     }),
#     exch_mdl_red = map(matched_data, function(df) {
#       get_gee_results(df, y ~ visit * bav, "exchangeable")
#     })
#   )

save(gee_fits_df, file = "Outputs/GEE_results.RData")
