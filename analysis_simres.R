# title: Analysis of GEE and QLS fit results on simulated data
# author: Jinyu Luo
# version: 2024-06-30

library(janitor)
library(tidyverse)

# Load simulation results ------------------------------------------------------
load("data/realcoefs.RData")
load("Outputs/GEE_results.RData")
load("Outputs/qls_results.RData")
true_coefs <- full_ar1$coefficients$Estimate
names(true_coefs) <- c("g0", "bav", "visit", "age", "male", "bsa", "bav_visit")
corr_alpha <- full_ar1$corr$Estimate
rho <- 0.3

# Calculate No.matched pairs in each simulation 
# matched_info <- data.frame(sim_id = 1:n_sim, 
#                            n_patients = map_int(matched_df, ~ n_distinct(.x$id)), 
#                            pairs = map_int(matched_df, ~ n_distinct(.x$matched)))
# 
# matched_info %>% 
#   ggplot(aes(x = n_patients)) +
#   geom_histogram(binwidth = 5, fill = "steelblue", color = "black", alpha = 0.7) +
#   labs(title = "Distribution of the number of matched patients among 1000 simulated datasets", 
#        x = "Number of Distinct Patients", y = "Frequency")+
#   theme_minimal()
# 
# matched_info %>% 
#   ggplot(aes(x = pairs)) +
#   geom_histogram(binwidth = 5, fill = "steelblue", color = "black", alpha = 0.7) +
#   labs(title = "Distribution of the number of matched patients among 1000 simulated datasets", 
#        x = "Number of Matched Pairs", y = "Frequency")+
#   theme_minimal()

# Helper functions -------------------------------------------------------------
# 1. Function to calculate coverage for each simulation 
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

# 2. Function to calculate mean estimates for each term
calculate_mean_estimates <- function(results_list) {
  results_list %>%
    map_dfr(~.x %>% 
              select(term, estimate)) %>%
    group_by(term) %>%
    summarize(mean_estimate = mean(estimate, na.rm = TRUE), 
              std = sd(estimate, na.rm=TRUE)) 
}


# 3. Function to calculate mean rho and rho_se for AR1 and Exchangeable 
calculate_mean_rho <- function(results_list) {
  results_list %>%
    map_dfr(~ .x %>% select(rho) %>% distinct()) %>%
    summarise(mean_rho = mean(rho))
}


## Analysis --------------------------------------------------------------------
se_max <- 5

### GEE ------------------------------------------------------------------------
# 1. Convergence 
gee_convergence_df <- gee_fits_df %>%
  mutate(ind_full = map_lgl(ind_mdl_full, ~any(.x$convergence == TRUE)),
         ar1_full = map_lgl(ar1_mdl_full, ~any(.x$convergence == TRUE)),
         exch_full = map_lgl(exch_mdl_full, ~any(.x$convergence == TRUE)),
         ind_red = map_lgl(ind_mdl_red, ~any(.x$convergence == TRUE)), 
         ar1_red = map_lgl(ar1_mdl_red, ~any(.x$convergence == TRUE)), 
         exch_red = map_lgl(exch_mdl_red, ~any(.x$convergence == TRUE))) %>% 
  select(sim_id, ind_full:exch_red) 

gee_convergence <- gee_convergence_df %>% select(-sim_id) %>% colMeans()

# Filter out non-converged simulations
gee_non_converged <- gee_convergence_df %>% 
  filter(if_all(ind_full:exch_red, ~ . == FALSE)) %>% 
  select(sim_id, ind_full:exch_red) 
print(paste("There are", nrow(gee_non_converged), "simulations did not converge for GEE fit."))

gee_sim_results <- gee_fits_df %>% filter(!(sim_id %in% gee_non_converged$sim_id))


# 2. Check Extreme standard errors 
# Remove Simulations with SE > 5 
se_max <- 5 
gee_extreme_se <- gee_sim_results %>% 
  mutate(ind_full_se = map_lgl(ind_mdl_full, ~any(.x$std_error > se_max)),
         ar1_full_se = map_lgl(ar1_mdl_full, ~any(.x$std_error > se_max)),
         exch_full_se = map_lgl(exch_mdl_full, ~any(.x$std_error > se_max)),
         ind_red_se = map_lgl(ind_mdl_red, ~any(.x$std_error > se_max)), 
         ar1_red_se = map_lgl(ar1_mdl_red, ~any(.x$std_error > se_max)), 
         exch_red_se = map_lgl(exch_mdl_red, ~any(.x$std_error > se_max))) %>% 
  select(sim_id, ind_full_se:exch_red_se) %>% 
  mutate(across(everything(), ~replace_na(., TRUE))) %>% 
  filter(if_any(-sim_id, ~ . == TRUE | is.na(.))) 

gee_extreme_se %>% select(-sim_id) %>% colMeans(na.rm = TRUE)

gee_sim_results <- gee_sim_results %>% filter(!(sim_id %in% gee_extreme_se$sim_id))
print(paste(1000-nrow(gee_sim_results), "simulations with extreme S.E. are removed for GEE fit."))


## 3. True Parameter Coverage Probability -----------------------------------------
true_set <- data.frame(term = c("(Intercept)", "visit", "bav", "age", "male", "bsa", "visit:bav"), 
                       coefficient = c(true_coefs[1], true_coefs[3], true_coefs[2], 
                                       true_coefs[4:7]), 
                       row.names = NULL)
# (1) Coverage probabilities for models without adjustment by age, male, and BSA
coverage_gee_unadj <- calculate_coverage(gee_sim_results$ind_mdl_red, true_set) %>% 
  rbind(calculate_coverage(gee_sim_results$ar1_mdl_red, true_set)) %>% 
  rbind(calculate_coverage(gee_sim_results$exch_mdl_red, true_set)) %>% 
  na.omit()

# (2) Coverage probabilities for models adjusted by age, male, and BSA
coverage_gee_adj <- calculate_coverage(gee_sim_results$ind_mdl_full, true_set) %>% 
  rbind(calculate_coverage(gee_sim_results$ar1_mdl_full, true_set)) %>% 
  rbind(calculate_coverage(gee_sim_results$exch_mdl_full, true_set)) %>% 
  na.omit()


# 4. Mean GEE Fit 
# (1) Unadjusted Models 

true_set %>% 
  right_join(calculate_mean_estimates(gee_sim_results$ind_mdl_red), by = "term")

  


meanest_gee_unadj <- true_set %>% rename(True = coefficient) %>% 
  filter(!(term %in% c("age", "male", "bsa"))) %>% 
  left_join(calculate_mean_estimates(gee_sim_results$ind_mdl_red) %>% 
              rename(Independence = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(gee_sim_results$ar1_mdl_red) %>% rename(`AR1` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(gee_sim_results$exch_mdl_red) %>% rename(`Exchangeable` = mean_estimate), by = "term")

meanest_gee_unadj %>% 
  mutate(bias_ind = Independence-True, 
         bias_ar1 = AR1-True, 
         bias_exch = Exchangeable-True, 
         rel_bias_ind = bias_ind/True, 
         rel_bias_ar1 = bias_ar1/True, 
         rel_bias_exch = bias_exch/True) %>% 
  relocate(c(bias_ind,rel_bias_ind), .after = Independence) %>% 
  relocate(c(bias_ar1,rel_bias_ar1), .after = AR1) %>% 
  relocate(c(bias_exch,rel_bias_exch), .after = Exchangeable)

  
# (2) Adjusted Models 
meanest_gee_adj <- true_set %>% rename(True = coefficient) %>% 
  left_join(calculate_mean_estimates(gee_sim_results$ind_mdl_full) %>% rename(`Independence` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(gee_sim_results$ar1_mdl_full) %>% rename(`AR1` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(gee_sim_results$exch_mdl_full) %>% rename(`Exchangeable` = mean_estimate), by = "term")


#5. Check the correlation parameters 
mean_rho_gee <- gee_sim_results %>%
  select(ar1_mdl_full, ar1_mdl_red, exch_mdl_full, exch_mdl_red) %>% 
  summarise(across(everything(), calculate_mean_rho, .names = "mean_{col}")) %>%
  pivot_longer(everything(), names_to = "model_type", values_to = "mean_values") %>%
  unnest(mean_values) %>% 
  mutate(corstr = c("AR1", "AR1", "Exchangeable", "Exchangeable"), 
         covariates = c("Adjusted", "Unadjusted", "Adjusted", "Unadjusted"),
         true_rho = rep(rho, 4)) %>% 
  relocate(c(corstr, covariates), .after = model_type) %>% 
  relocate(true_rho, .after = mean_rho) %>% 
  select(-model_type)


### QLS ------------------------------------------------------------------------
# 1. Convergence 
qls_convergence_df <- qls_fits_df %>%
  mutate(ind_full = map_lgl(ind_mdl_full, ~any(.x$convergence == TRUE)),
         ar1_full = map_lgl(ar1_mdl_full, ~any(.x$convergence == TRUE)),
         exch_full = map_lgl(exch_mdl_full, ~any(.x$convergence == TRUE)),
         ind_red = map_lgl(ind_mdl_red, ~any(.x$convergence == TRUE)), 
         ar1_red = map_lgl(ar1_mdl_red, ~any(.x$convergence == TRUE)), 
         exch_red = map_lgl(exch_mdl_red, ~any(.x$convergence == TRUE))) %>% 
  select(sim_id, ind_full:exch_red) 

qls_convergence <- qls_convergence_df %>% select(-sim_id) %>% colMeans()

# Filter out non-converged simulations
qls_non_converged <- qls_convergence_df %>% 
  filter(if_all(ind_full:exch_red, ~ . == FALSE)) %>% 
  select(sim_id, ind_full:exch_red) 
print(paste("There are", nrow(gee_non_converged), "simulations did not converge for GEE fit."))

qls_sim_results <- qls_fits_df %>% filter(!(sim_id %in% qls_non_converged$sim_id))


# 2. Check Extreme standard errors 
# Remove Simulations with SE > 5 
se_max <- 5 
qls_extreme_se <- qls_sim_results %>% 
  mutate(ind_full_se = map_lgl(ind_mdl_full, ~any(.x$std_error > se_max)),
         ar1_full_se = map_lgl(ar1_mdl_full, ~any(.x$std_error > se_max)),
         exch_full_se = map_lgl(exch_mdl_full, ~any(.x$std_error > se_max)),
         ind_red_se = map_lgl(ind_mdl_red, ~any(.x$std_error > se_max)), 
         ar1_red_se = map_lgl(ar1_mdl_red, ~any(.x$std_error > se_max)), 
         exch_red_se = map_lgl(exch_mdl_red, ~any(.x$std_error > se_max))) %>% 
  select(sim_id, ind_full_se:exch_red_se) %>% 
  mutate(across(everything(), ~replace_na(., TRUE))) %>% 
  filter(if_any(-sim_id, ~ . == TRUE | is.na(.))) 

qls_extreme_se %>% select(-sim_id) %>% colMeans(na.rm = TRUE)

qls_sim_results <- qls_sim_results %>% filter(!(sim_id %in% qls_extreme_se$sim_id))
print(paste(1000-nrow(qls_sim_results), "simulation(s) with extreme S.E. are removed for GEE fit."))
print(paste("Simulation removed:", qls_extreme_se$sim_id))


# 3. Mean GEE Fit 
# (1) Unadjusted Models 
meanest_qls_unadj <- true_set %>% rename(True = coefficient) %>% 
  filter(!(term %in% c("age", "male", "bsa"))) %>% 
  left_join(calculate_mean_estimates(qls_sim_results$ind_mdl_red) %>% rename(`Independence` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(qls_sim_results$ar1_mdl_red) %>% rename(`AR1` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(qls_sim_results$exch_mdl_red) %>% rename(`Exchangeable` = mean_estimate), by = "term")

# (2) Adjusted Models 
meanest_qls_adj <- true_set %>% rename(True = coefficient) %>% 
  left_join(calculate_mean_estimates(qls_sim_results$ind_mdl_full) %>% rename(`Independence` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(qls_sim_results$ar1_mdl_full) %>% rename(`AR1` = mean_estimate), by = "term") %>% 
  left_join(calculate_mean_estimates(qls_sim_results$exch_mdl_full) %>% rename(`Exchangeable` = mean_estimate), by = "term")


## Check the correlation parameters --------------------------------------------
mean_rho_qls <- qls_sim_results %>%
  select(ar1_mdl_full, ar1_mdl_red, exch_mdl_full, exch_mdl_red) %>% 
  summarise(across(everything(), calculate_mean_rho, .names = "mean_{col}")) %>%
  pivot_longer(everything(), names_to = "model_type", values_to = "mean_values") %>%
  unnest(mean_values) %>% 
  mutate(corstr = c("AR1", "AR1", "Exchangeable", "Exchangeable"), 
         covariates = c("Adjusted", "Unadjusted", "Adjusted", "Unadjusted"),
         true_rho = rep(rho, 4)) %>% 
  relocate(c(corstr, covariates), .after = model_type) %>% 
  relocate(true_rho, .after = mean_rho) %>% 
  select(-model_type)



