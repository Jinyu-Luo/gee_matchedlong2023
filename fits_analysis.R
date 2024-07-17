library(tidyverse)
library(kableExtra)
load("~/Documents/Projects/Practicum_GEE/Outputs/simfit.RData")
load("~/Documents/Projects/Practicum_GEE/data/realcoefs.RData")
load("~/Documents/Projects/Practicum_GEE/data/qlsfits.RData")

true_coefs <- red_ar1$coefficients$Estimate
names(true_coefs) <- c("g0", "bav", "visit", "bav_visit")
corr_alpha <- red_ar1$corr$Estimate
rho <- 0.3
se_max <- 5

true_set <- data.frame(term = c("(Intercept)", "visit", "bav", "bav:visit"), 
                       coefficient = true_coefs, 
                       row.names = NULL)

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
evaluate_fits <- function(df, true_set) {
  df %>%
    map_dfr(~.x %>% select(term, estimate, std_error, adj_std_error)) %>%
    group_by(term) %>%
    summarize(mean_estimate = mean(estimate, na.rm = TRUE), 
              est_std = sd(estimate, na.rm=TRUE), 
              mean_se = mean(std_error, na.rm=TRUE), 
              mean_adj_se = mean(adj_std_error, na.rm = TRUE)) %>% 
    left_join(true_set, by = "term") %>% 
    mutate(mean_bias = mean_estimate - coefficient, 
           mean_rel_bias = mean_bias/abs(coefficient)*100, 
           mean_mse = est_std^2+mean_bias^2) %>% 
    select(term, coefficient, mean_estimate, est_std, mean_se, mean_adj_se, mean_bias, mean_rel_bias, mean_mse)
}

# Function to calculate relative bias for each simulation 
cal_rel_bias <- function(df, true_set){
 df %>% 
    map_dfr(~.x %>% select(term, estimate)) %>%
    left_join(true_set, by = "term") %>% 
    group_by(term) %>% 
    mutate(bias = estimate - coefficient, 
           rel_bias = bias/abs(coefficient), 
           sim = row_number()) %>% 
   ungroup() %>% na.omit() %>% 
   select(sim, term, rel_bias) %>% 
   pivot_wider(names_from = term, values_from = rel_bias) %>% 
   rename(b0 = `(Intercept)`, 
          b1 = bav, 
          b2 = visit, 
          b3 = `bav:visit`)
}

# 3. Function to calculate mean rho and rho_se for AR1 and Exchangeable 
calculate_mean_rho <- function(df) {
  df %>%
    map_dfr(~ .x %>% select(rho) %>% distinct()) %>%
    summarise(mean_rho = mean(rho), 
              sd_rho = sd(rho)) %>% 
    mutate(bias = mean_rho - 0.3, 
           rel_bias = bias/0.3, 
           mse = sd_rho^2+bias^2)
}

# Function to extract tau from QLS estimation 
pull_tau <- function(df){
  df %>%
    map_dfr(~ .x %>% select(tau) %>% distinct())  %>% 
    pull(tau)
}


# 4. Function for calculating convergence 
calculate_convergence <- function(fits_df){
  result <- fits_df %>%
    mutate(
      ind_full = map_lgl(ind_mdl_full, ~any(.x$convergence == TRUE)),
      ar1_full = map_lgl(ar1_mdl_full, ~any(.x$convergence == TRUE)),
      exch_full = map_lgl(exch_mdl_full, ~any(.x$convergence == TRUE)),
      ind_red = map_lgl(ind_mdl_red, ~any(.x$convergence == TRUE)), 
      ar1_red = map_lgl(ar1_mdl_red, ~any(.x$convergence == TRUE)), 
      exch_red = map_lgl(exch_mdl_red, ~any(.x$convergence == TRUE))
    ) %>%
    select(sim_id, ind_full:exch_red)
  
  return(result) 
}

# 5. Function for searching simulations with extreme estimate S.E. 
extr_se_search <- function(sim_results) {
  extreme_df <- sim_results %>% 
    mutate(ind_full_se = map_lgl(ind_mdl_full, ~any(.x$std_error > se_max)),
           ar1_full_se = map_lgl(ar1_mdl_full, ~any(.x$std_error > se_max)),
           exch_full_se = map_lgl(exch_mdl_full, ~any(.x$std_error > se_max)),
           ind_red_se = map_lgl(ind_mdl_red, ~any(.x$std_error > se_max)), 
           ar1_red_se = map_lgl(ar1_mdl_red, ~any(.x$std_error > se_max)), 
           exch_red_se = map_lgl(exch_mdl_red, ~any(.x$std_error > se_max))) %>% 
    select(sim_id, ind_full_se:exch_red_se) %>% 
    mutate(across(everything(), ~replace_na(., TRUE))) %>% 
    filter(if_any(-sim_id, ~ . == TRUE | is.na(.))) 
  
  n_extreme <- nrow(extreme_df)
  if(n_extreme != 0){
    paste(print("There are"), n_extreme, "simulations have extreme values. ")
    sim_results <- sim_results %>% filter(!(sim_id %in% extreme_df$sim_id))
    print(paste(1000-nrow(sim_results), "simulations with extreme S.E. are removed for GEE fit."))
    print("Removed Simulation ID:")
    print(unlist(extreme_df$sim_id))
  } else {print("No extreme S.E.")}
  return(sim_results)
}

extract_term_data <- function(df, term, mod) {
  output <- df %>%
    select(sim_id, ends_with(mod)) %>%
    pivot_longer(cols = ends_with(mod), names_to = "model", values_to = "data") %>%
    rowwise() %>%
    mutate(data = list(filter(data, term == !!term))) %>%
    unnest(cols = c(data)) %>%
    ungroup() %>%
    mutate(sim_id = as.integer(sim_id), 
           model = case_when(str_detect(model, "ind") ~ "Independence", 
                             str_detect(model, "ar1") ~ "AR1", 
                             TRUE ~ "Exchangeable")) %>% 
    select(sim_id, model, term, estimate, std_error, lower, upper, adj_lower, adj_upper)
  return(output)
}

# Function to generate plot for each term
generate_plot <- function(term_data) {
  ggplot(term_data, aes(y = mean_estimate, x = interaction(type, adjusted), 
                        color = model)) + 
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, 
                  position = position_dodge(width = 0.5)) +
    theme_bw() +
    labs(x = "", y = expression(beta), color = "") +
    scale_color_manual(values =  c("#17365d", "#548cd4","#8cb4e2")) + # Darker sand colors
    theme(strip.background = element_rect(fill = "#17365d"),
          strip.text = element_text(color = "white"), 
          legend.position="bottom")
}


# GEE --------------------------------------------------------------------------
# gee_convergence_df <- calculate_convergence(gee_fits_df)
# gee_convergence <- gee_convergence_df %>% select(-sim_id) %>% colMeans()
# # Filter out non-converged simulations
# gee_non_converged <- gee_convergence_df %>%
#   filter(if_all(ind_full:exch_red, ~ . == FALSE)) %>%
#   select(sim_id, ind_full:exch_red)
# gee_conv_info <- paste("There are", nrow(gee_non_converged), "simulation(s) did not converge for GEE fit.")
# gee_sim_results <- gee_fits_df %>% filter(!(sim_id %in% gee_non_converged$sim_id))

# QLS --------------------------------------------------------------------------
# qls_convergence_df <- calculate_convergence(qls_fits_df)
# qls_convergence <- qls_convergence_df %>% select(-sim_id) %>% colMeans()
# # Filter out non-converged simulations
# qls_non_converged <- qls_convergence_df %>%
#   filter(if_all(ind_full:exch_red, ~ . == FALSE)) %>%
#   select(sim_id, ind_full:exch_red)
# qls_conv_info <- paste("There are", nrow(gee_non_converged), "simulations did not converge for GEE fit.")
# qls_sim_results <- qls_fits_df %>% filter(!(sim_id %in% qls_non_converged$sim_id))
# 
# gee_sim_results <- extr_se_search(gee_sim_results)
# qls_sim_results <- extr_se_search(qls_sim_results)

# GEE 
# (1) Coverage probabilities for models without adjustment by age, male, and BSA
# coverage_gee_unadj <- calculate_coverage(gee_sim_results$ind_mdl_red, true_set) %>%
#   rbind(calculate_coverage(gee_sim_results$ar1_mdl_red, true_set)) %>%
#   rbind(calculate_coverage(gee_sim_results$exch_mdl_red, true_set)) %>%
#   na.omit()
# 
# # (2) Coverage probabilities for models adjusted by age, male, and BSA
# coverage_gee_adj <- calculate_coverage(gee_sim_results$ind_mdl_full, true_set) %>%
#   rbind(calculate_coverage(gee_sim_results$ar1_mdl_full, true_set)) %>%
#   rbind(calculate_coverage(gee_sim_results$exch_mdl_full, true_set)) %>%
#   na.omit()
## Significant improvement by including covariates for adjustment 

# QLS 
# # (1) Coverage probabilities for models without adjustment by age, male, and BSA
# coverage_qls_unadj <- calculate_coverage(qls_sim_results$ind_mdl_red, true_set) %>%
#   rbind(calculate_coverage(qls_sim_results$ar1_mdl_red, true_set)) %>%
#   rbind(calculate_coverage(qls_sim_results$exch_mdl_red, true_set)) %>%
#   na.omit()
# 
# # (2) Coverage probabilities for models adjusted by age, male, and BSA
# coverage_qls_adj <- calculate_coverage(qls_sim_results$ind_mdl_full, true_set) %>%
#   rbind(calculate_coverage(qls_sim_results$ar1_mdl_full, true_set)) %>%
#   rbind(calculate_coverage(qls_sim_results$exch_mdl_full, true_set)) %>%
#   na.omit()


# bav_red <- extract_term_data(gee_sim_results, term = "bav", mod = "_red") %>% mutate(type = "GEE") %>%
#   rbind(extract_term_data(qls_sim_results, term = "bav", mod = "_red") %>% mutate(type = "QLS")) %>%
#   relocate(type, .after = sim_id)
# 
# visit_red <- extract_term_data(gee_sim_results, term = "visit", mod = "_red") %>% mutate(type = "GEE") %>%
#   rbind(extract_term_data(qls_sim_results, term = "visit", mod = "_red") %>% mutate(type = "QLS")) %>%
#   relocate(type, .after = sim_id)
# 
# bav_visit_red <- extract_term_data(gee_sim_results, term = "bav:visit", mod = "_red") %>% mutate(type = "GEE") %>%
#   rbind(extract_term_data(qls_sim_results, term = "bav:visit", mod = "_red") %>% mutate(type = "QLS")) %>%
#   relocate(type, .after = sim_id)
# 
# 
# red_mdl_df <- rbind(bav_red, visit_red, bav_visit_red)


# bav_full <- extract_term_data(gee_sim_results, term = "bav", mod = "_full") %>% mutate(type = "GEE") %>%
#   rbind(extract_term_data(qls_sim_results, term = "bav", mod = "_full") %>% mutate(type = "QLS")) %>%
#   relocate(type, .after = sim_id)
# 
# visit_full <- extract_term_data(gee_sim_results, term = "visit", mod = "_full") %>% mutate(type = "GEE") %>%
#   rbind(extract_term_data(qls_sim_results, term = "visit", mod = "_full") %>% mutate(type = "QLS")) %>%
#   relocate(type, .after = sim_id)
# 
# bav_visit_full <- extract_term_data(gee_sim_results, term = "bav:visit", mod = "_full") %>% mutate(type = "GEE") %>%
#   rbind(extract_term_data(qls_sim_results, term = "bav:visit", mod = "_full") %>% mutate(type = "QLS")) %>%
#   relocate(type, .after = sim_id)
# 
# age_df <- extract_term_data(gee_sim_results, term = "age", mod = "_full") %>% mutate(type = "GEE") %>%
#   rbind(extract_term_data(qls_sim_results, term = "age", mod = "_full") %>% mutate(type = "QLS")) %>%
#   relocate(type, .after = sim_id)
# 
# male_df <- extract_term_data(gee_sim_results, term = "male", mod = "_full") %>% mutate(type = "GEE") %>%
#   rbind(extract_term_data(qls_sim_results, term = "male", mod = "_full") %>% mutate(type = "QLS")) %>%
#   relocate(type, .after = sim_id)
# 
# bsa_df <- extract_term_data(gee_sim_results, term = "bsa", mod = "_full") %>% mutate(type = "GEE") %>%
#   rbind(extract_term_data(qls_sim_results, term = "bsa", mod = "_full") %>% mutate(type = "QLS")) %>%
#   relocate(type, .after = sim_id)


# full_mdl_df <- rbind(bav_full, visit_full, bav_visit_full, age_df, male_df, bsa_df)

# red_summary <- red_mdl_df %>%
#   mutate(term = case_when(term == "bav" ~ "BAV",
#                           term == "visit" ~ "Visit",
#                           TRUE ~ "BAV:Visit")) %>%
#   group_by(term, model, type) %>%
#   summarise(mean_estimate = mean(estimate),
#             mean_lower = mean(lower),
#             mean_upper = mean(upper),
#             mean_adj_lower = mean(adj_lower),
#             mean_adj_upper = mean(adj_upper),
#             .groups = 'drop') %>%
#   pivot_longer(cols = c(mean_lower, mean_upper, mean_adj_lower, mean_adj_upper),
#                names_to = "ci_type", values_to = "ci_value") %>%
#   mutate(ci_type = case_when(
#     ci_type == "mean_lower" ~ "Unadjusted Lower",
#     ci_type == "mean_upper" ~ "Unadjusted Upper",
#     ci_type == "mean_adj_lower" ~ "Adjusted Lower",
#     ci_type == "mean_adj_upper" ~ "Adjusted Upper")) %>%
#   pivot_wider(names_from = ci_type, values_from = ci_value) %>%
#   pivot_longer(cols = c("Unadjusted Lower", "Unadjusted Upper", "Adjusted Lower", "Adjusted Upper"),
#                names_to = "ci_type", values_to = "ci_value") %>%
#   separate(ci_type, into = c("adjusted", "boundary"), sep = " ") %>%
#   pivot_wider(names_from = boundary, values_from = ci_value)


# full_summary <- full_mdl_df %>%
#   mutate(term = case_when(term == "bav" ~ "BAV",
#                           term == "visit" ~ "Visit",
#                           term == "age" ~ "Age",
#                           term == "male" ~ "Male",
#                           term == "bsa" ~ "BSA",
#                           TRUE ~ "BAV:Visit")) %>%
#   group_by(term, model, type) %>%
#   summarise(mean_estimate = mean(estimate),
#             mean_lower = mean(lower),
#             mean_upper = mean(upper),
#             mean_adj_lower = mean(adj_lower),
#             mean_adj_upper = mean(adj_upper),
#             .groups = 'drop') %>%
#   pivot_longer(cols = c(mean_lower, mean_upper, mean_adj_lower, mean_adj_upper),
#                names_to = "ci_type", values_to = "ci_value") %>%
#   mutate(ci_type = case_when(
#     ci_type == "mean_lower" ~ "Unadjusted Lower",
#     ci_type == "mean_upper" ~ "Unadjusted Upper",
#     ci_type == "mean_adj_lower" ~ "Adjusted Lower",
#     ci_type == "mean_adj_upper" ~ "Adjusted Upper")) %>%
#   pivot_wider(names_from = ci_type, values_from = ci_value) %>%
#   pivot_longer(cols = c("Unadjusted Lower", "Unadjusted Upper", "Adjusted Lower", "Adjusted Upper"),
#                names_to = "ci_type", values_to = "ci_value") %>%
#   separate(ci_type, into = c("adjusted", "boundary"), sep = " ") %>%
#   pivot_wider(names_from = boundary, values_from = ci_value)


# Generate plots for each term
# red_mdl_plots <- red_summary %>%
#   split(.$term) %>%
#   lapply(generate_plot)
# 
# full_mdl_plots <- full_summary %>%
#   split(.$term) %>%
#   lapply(generate_plot)


# save(coverage_gee_adj, coverage_gee_unadj, coverage_qls_adj, coverage_qls_unadj,
#      gee_sim_results, qls_sim_results, true_set, red_mdl_df, full_mdl_df,
#      red_summary, full_summary, red_mdl_plots, full_mdl_plots,
#      gee_conv_info, qls_conv_info,
#      file = "Outputs/fits_analysis.RData")

