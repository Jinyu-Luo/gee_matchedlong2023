library(tidyverse)
library(kableExtra)
load("~/Documents/Projects/Practicum_GEE/Outputs/simfit.RData")
load("~/Documents/Projects/Practicum_GEE/data/realcoefs.RData")
load("~/Documents/Projects/Practicum_GEE/data/qlsfits.RData")

true_set <- red_ar1$coefficients %>% select(Estimate) %>% 
  rename(truth = Estimate) %>% 
  rownames_to_column(var = "term")

# Helper functions -------------------------------------------------------------
# 1. Function to calculate coverage for each simulation 
calculate_coverage <- function(fit_list, true_set) {
  fit_list %>%
    map_dfr(~.x %>%
              right_join(true_set, by = "term") %>%
              mutate(unadjusted = truth >= lower & truth <= upper,
                     adjusted = truth >= adj_lower & truth <= adj_upper)) %>%
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
    mutate(mean_bias = mean_estimate - truth, 
           mean_rel_bias = mean_bias/abs(truth)*100, 
           mean_mse = est_std^2+mean_bias^2) %>% 
    select(term, truth, mean_estimate, est_std, mean_se, mean_adj_se, mean_bias, mean_rel_bias, mean_mse)
}

# Function to calculate relative bias for each simulation 
cal_rel_bias <- function(df, true_set){
 df %>% 
    map_dfr(~.x %>% select(term, estimate)) %>%
    left_join(true_set, by = "term") %>% 
    group_by(term) %>% 
    mutate(bias = estimate - truth, 
           rel_bias = bias/abs(truth), 
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
  temp <- fits_df %>%
    mutate(
      ind_full = map_lgl(ind_mdl_full, ~any(.x$convergence == TRUE)),
      ar1_full = map_lgl(ar1_mdl_full, ~any(.x$convergence == TRUE)),
      exch_full = map_lgl(exch_mdl_full, ~any(.x$convergence == TRUE)),
      ind_red = map_lgl(ind_mdl_red, ~any(.x$convergence == TRUE)), 
      ar1_red = map_lgl(ar1_mdl_red, ~any(.x$convergence == TRUE)), 
      exch_red = map_lgl(exch_mdl_red, ~any(.x$convergence == TRUE))
    ) %>%
    select(sim_id, ind_full:exch_red)
  
  converged <- temp %>% select(-sim_id) %>% colMeans()
  diverged <- temp %>% filter(if_all(ind_full:exch_red, ~ . == FALSE)) 
  if(nrow(diverged) != 0){
    diverged_sim <- diverged %>% pull(sim_id)
    out_data <- fits_df %>% filter(!(sim_id %in% diverged))
  } else{
    out_data <- fits_df
    diverged_sim <- NULL 
  }
  return(list(data = out_data, diverged = diverged_sim)) 
}

# 5. Function for searching simulations with extreme estimate S.E. 
extr_se_search <- function(sim_results, se_max) {
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
  messages <- character()
  
  if(n_extreme != 0){
    messages <- c(messages, paste("There are", n_extreme, "simulations have extreme values."))
    sim_results <- sim_results %>% filter(!(sim_id %in% extreme_df$sim_id))
    messages <- c(messages, paste(1500 - nrow(sim_results), "simulations with extreme S.E. are removed for the fit."))
    messages <- c(messages, "Removed Simulation ID:")
    messages <- c(messages, paste(unlist(extreme_df$sim_id), collapse = ", "))
  } else {
    messages <- c(messages, "No extreme S.E.")
  }
  
  list(
    filtered_results = sim_results,
    messages = messages, 
    extreme_id = unlist(extreme_df$sim_id)
  )
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
