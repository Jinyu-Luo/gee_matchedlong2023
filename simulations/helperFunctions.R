# title: Simulation Helper Functions 
# author: Jinyu Luo
# version: 2024-07-18

rm(list = ls())
# Required R packages ----------------------------------------------------------
library(CorBin)
library(tidyverse)
library(geepack)
library(parallel)
library(survival)
library(simsurv)
load("data/realcoefs.RData")

# Global Variables -------------------------------------------------------------
n_patients <- 250 
n_sim <- 1500
maxT <- 6
rho <- 0.3
alpha_ci <- 0.05

true_coefs <- red_ar1$coefficients$Estimate
names(true_coefs) <- c("g0", "bav", "visit", "bav_visit")
corr_alpha <- red_ar1$corr$Estimate
surv_coefs <- surv_coefs[,-4]
rownames(surv_coefs)[4] <- "male"
rho <- 0.3
se_max <- 5

true_set <- data.frame(term = c("(Intercept)", "visit", "bav", "bav:visit"), 
                       coefficient = true_coefs, 
                       row.names = NULL)

# Functions for GEE and QLS Fit ------------------------------------------------
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

# Function to fit GEE models 
fit_gee <- function(df){
  fit_res <- list()
  for (i in 1:length(df)) {
    dat <- df[[i]]
    sim_id <- i 
    mdl_res <- list()
    for (mdl in names(model_specs)) {
      m <- model_specs[[mdl]]
      mdl_res[[mdl]] <- get_gee_results(df = dat, formula = m$formula, 
                                        corstr = m$corstr, adjusted = m$adjusted)
    }
    fit_res[[i]] <- c(list(sim_id = sim_id), mdl_res)
    print(paste("Completed simulation", i, "out of", n_sim))
  }
  
  output <- tibble(
    sim_id = map(fit_res, "sim_id"),
    ind_mdl_full = map(fit_res, "ind_mdl_full"),
    ind_mdl_red = map(fit_res, "ind_mdl_red"), 
    ar1_mdl_full = map(fit_res, "ar1_mdl_full"), 
    ar1_mdl_red = map(fit_res, "ar1_mdl_red"), 
    exch_mdl_full = map(fit_res, "exch_mdl_full"), 
    exch_mdl_red = map(fit_res, "exch_mdl_red")
  )
  return(output)
}



# Functions for Fit Analysis ---------------------------------------------------
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
    dplyr::filter(if_any(-sim_id, ~ . == TRUE | is.na(.))) 
  
  n_extreme <- nrow(extreme_df)
  if(n_extreme != 0){
    paste(print("There are"), n_extreme, "simulations have extreme values. ")
    sim_results <- sim_results %>% dplyr::filter(!(sim_id %in% extreme_df$sim_id))
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
    mutate(data = list(dplyr::filter(data, term == !!term))) %>%
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
