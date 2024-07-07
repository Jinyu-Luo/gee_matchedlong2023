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

# 2. Function to calculate coverage for each simulation 
calculate_coverage <- function(summary_list, true_set) {
  summary_list %>%
    map_dfr(~.x %>%
              left_join(true_set, by = "term") %>%
              mutate(unadjusted = coefficient >= lower & coefficient <= upper,
                     adjusted = coefficient >= adj_lower & coefficient <= adj_upper)) %>%
    group_by(term) %>%
    summarize(unadjusted = mean(unadjusted, na.rm = TRUE),
              adjusted = mean(adjusted, na.rm = TRUE))
}

# 3. Function to calculate mean estimates for each term
calculate_mean_estimates <- function(results_list) {
  results_list %>%
    map_dfr(~.x %>% 
              select(term, estimate)) %>%
    group_by(term) %>%
    summarize(mean_estimate = mean(estimate, na.rm = TRUE))
}

# 4. Function to calculate mean rho and rho_se for AR1 and Exchangeable 
calculate_mean_rho <- function(results_list) {
  results_list %>%
    map_dfr(~ .x %>% select(rho, rho_se) %>% distinct()) %>%
    summarise(mean_rho = mean(rho), mean_rho_se = mean(rho_se))
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


## Check Convergence ------------------------------------------------------------
convergence <- gee_fits_df %>%
  mutate(ind_full = map_lgl(ind_mdl_full, ~any(.x$convergence == TRUE)),
         ar1_full = map_lgl(ar1_mdl_full, ~any(.x$convergence == TRUE)),
         exch_full = map_lgl(exch_mdl_full, ~any(.x$convergence == TRUE)),
         ind_red = map_lgl(ind_mdl_red, ~any(.x$convergence == TRUE)), 
         ar1_red = map_lgl(ar1_mdl_red, ~any(.x$convergence == TRUE)), 
         exch_red = map_lgl(exch_mdl_red, ~any(.x$convergence == TRUE))) %>% 
  select(sim_id, ind_full:exch_red) 

convergence_result <- convergence %>% select(-sim_id) %>% colMeans()

# Filter out non-converged simulations
non_converged <- convergence %>% 
  filter(if_all(ind_full:exch_red, ~ . == FALSE)) %>% 
  select(sim_id, ind_full:exch_red) 
nrow(non_converged)

sim_results <- gee_fits_df %>% filter(!(sim_id %in% non_converged$sim_id))
nrow(sim_results)

# Remove Simulations with SE > 5 
se_max <- 5 
extreme_sim <- sim_results %>% 
  mutate(ind_full_se = map_lgl(ind_mdl_full, ~any(.x$std_error > se_max)),
         ar1_full_se = map_lgl(ar1_mdl_full, ~any(.x$std_error > se_max)),
         exch_full_se = map_lgl(exch_mdl_full, ~any(.x$std_error > se_max)),
         ind_red_se = map_lgl(ind_mdl_red, ~any(.x$std_error > se_max)), 
         ar1_red_se = map_lgl(ar1_mdl_red, ~any(.x$std_error > se_max)), 
         exch_red_se = map_lgl(exch_mdl_red, ~any(.x$std_error > se_max))) %>% 
  select(sim_id, ind_full_se:exch_red_se) %>% 
  mutate(across(everything(), ~replace_na(., TRUE))) %>% 
  filter(if_any(-sim_id, ~ . == TRUE | is.na(.))) 

extreme_sim %>% select(-sim_id) %>% colMeans(na.rm = TRUE)

sim_results <- sim_results %>% filter(!(sim_id %in% extreme_sim$sim_id))
nrow(sim_results)
# [1] 958


## True Parameter Coverage Probability -----------------------------------------
true_set <- data.frame(term = c("(Intercept)", "visit", "bav", "age", "male", "bsa", "visit:bav"), 
                       coefficient = c(true_coefs[1], true_coefs[3], true_coefs[2], 
                                       true_coefs[4:7]), 
                       row.names = NULL)
# (1) Coverage probabilities for models without adjustment by age, male, and BSA
unadjusted_coverage <- calculate_coverage(sim_results$ind_mdl_red, true_set) %>% 
  rbind(calculate_coverage(sim_results$ar1_mdl_red, true_set)) %>% 
  rbind(calculate_coverage(sim_results$exch_mdl_red, true_set)) %>% 
  na.omit()

# (2) Coverage probabilities for models adjusted by age, male, and BSA
adjusted_coverage <- calculate_coverage(sim_results$ind_mdl_full, true_set) %>% 
  rbind(calculate_coverage(sim_results$ar1_mdl_full, true_set)) %>% 
  rbind(calculate_coverage(sim_results$exch_mdl_full, true_set)) %>% 
  na.omit()


## Mean GEE Fit ----------------------------------------------------------------
# (1) Unadjusted Models 
mean_est_unadj <- true_set %>% rename(True = coefficient) %>% 
  filter(!(term %in% c("age", "male", "bsa"))) %>% 
  left_join(calculate_mean_estimates(sim_results$ind_mdl_red) %>% rename(`Independence` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(sim_results$ar1_mdl_red) %>% rename(`AR1` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(sim_results$exch_mdl_red) %>% rename(`Exchangeable` = mean_estimate), by = "term")

# (2) Adjusted Models 
mean_est_adj <- true_set %>% rename(True = coefficient) %>% 
  left_join(calculate_mean_estimates(sim_results$ind_mdl_full) %>% rename(`Independence` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(sim_results$ar1_mdl_full) %>% rename(`AR1` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(sim_results$exch_mdl_full) %>% rename(`Exchangeable` = mean_estimate), by = "term")


# Calculate No.matched pairs in each simulation 
matched_info <- data.frame(sim_id = 1:n_sim, 
                           n_patients = map_int(matched_df, ~ n_distinct(.x$id)), 
                           pairs = map_int(matched_df, ~ n_distinct(.x$matched)))

matched_info %>% 
  ggplot(aes(x = n_patients)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of the number of matched patients among 1000 simulated datasets", 
       x = "Number of Distinct Patients", y = "Frequency")+
  theme_minimal()

matched_info %>% 
  ggplot(aes(x = pairs)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of the number of matched patients among 1000 simulated datasets", 
       x = "Number of Matched Pairs", y = "Frequency")+
  theme_minimal()


## Check the correlation parameters --------------------------------------------
mean_rho_values <- sim_results %>%
  select(ar1_mdl_full, ar1_mdl_red, exch_mdl_full, exch_mdl_red) %>% 
  summarise(across(everything(), calculate_mean_rho, .names = "mean_{col}")) %>%
  pivot_longer(everything(), names_to = "model_type", values_to = "mean_values") %>%
  unnest(mean_values) %>% 
  mutate(corstr = c("AR1", "AR1", "Exchangeable", "Exchangeable"), 
         covariates = c("Adjusted", "Unadjusted", "Adjusted", "Unadjusted"),
         true_rho = rep(rho, 4)) %>% 
  relocate(c(corstr, covariates), .after = model_type) %>% 
  select(-model_type)


save(convergence_result, true_set, adjusted_coverage, unadjusted_coverage, 
     mean_est_adj, mean_est_unadj, matched_info, mean_rho_values, 
     file = "Outputs/GEE_results.RData")
