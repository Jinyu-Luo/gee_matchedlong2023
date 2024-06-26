
rm(list = ls())
# install.packages("simsurv")
library(mvtnorm)
library(nlme)
library(simsurv)
library(CorBin)
library(tidyverse)
library(optmatch)
library(geepack)

## ---- Step 0: setup-sim ----
n <- 100
n_bav <- n*0.2
n_tav <- n*0.8
rho <- 0.6
alpha = 0.05

# 2. Coefficients for exposure group simulation 
params_sets <- list(
  psm_coefs = c(b0 = -0.24917, age = -0.07118, male = 1.24005, bsa = 2.94631), 
  gee_coefs = c(b0=-2.41129228, bav=-1.53758808, age=-0.03742342, 
                male=1.29569822, bsa=1.40495289, visit=-0.02214736, 
                bav_visit=0.31010373), 
  surv_coefs = c(bav=-0.124, age = -0.021, male=-0.139)
)



# Data Simulation ------------------------------------------------------
# 1. Function for simulate outcomes for 100 patients 
simulate_data <- function(coeffs){
  gee_b <- coeffs$gee_coefs
  psm_b <- coeffs$psm_coefs
  surv_b <- coeffs$surv_coefs
  patients <- data.frame(# baseline information simulation
    pid = 1:n,
    age = round(rnorm(n, mean = 50, sd = 10)),
    male = rbinom(n, size = 1, prob = 0.65),
    bsa = rnorm(n, mean = 2, sd = 0.3)) %>% 
    mutate(bav = if_else(pid <= n_bav, 1, 0)) %>% 
    relocate(bav, .before = age) %>% 
    expand_grid(timepoint = 1:6) %>% 
    arrange(pid, timepoint) %>% 
    mutate(bav_timepoint = timepoint*bav, 
           logit_y = gee_b['b0'] + bav * gee_b['bav'] + age * gee_b['age'] + 
             male * gee_b['male'] + bsa * gee_b['bsa'] + timepoint * gee_b['visit'] + 
             bav_timepoint * gee_b['bav_visit'], 
           prob_y = exp(logit_y)/(1+exp(logit_y)), 
           y = round(prob_y))
  
  baseline_info <- sim_patients %>% filter(timepoint==1)
  # baseline_info %>% group_by(y) %>% summarise(n=n())
  # A tibble: 2 × 2
  #      y     n
  #  <dbl> <int>
  #      0   430
  #      1    70
  # baseline_info %>% group_by(bav) %>% summarise(n=n())
  # A tibble: 2 × 2
  #    bav     n
  #  <dbl> <int>
  #      0   400
  #      1   100
  
  ## Simulate dropout timepoints -------------------------------------------------
  Z <- baseline_info %>% select(pid, bav, age, male) 
  
  # set.seed(5207)
  dropouts <- simsurv(lambdas = exp(0.611), # lambda: baseline hazard rate
                      gammas = exp(0.260), # gamma: shape parameter 
                      x = Z,
                      betas = surv_b,
                      maxt = 5, dist = "weibull") %>% 
    mutate(drop_time = round(eventtime)+1) %>% 
    rename(pid = id) %>% select(pid, drop_time)
  
  dropped_data <- patients %>% left_join(dropouts, by = "pid") %>% 
    filter(timepoint <= drop_time) %>% 
    select(-c(drop_time, logit_y, prob_y, bav_timepoint))
  
  # Propensity Score Matching --------------------------------------------------
  ps_model <- glm(bav ~ age + male + bsa, family = binomial, data = baseline_info)
  pps_match <- pairmatch(ps_model, data = baseline_info)
  matched_df <- data.frame(baseline_info, matched = pps_match, check.rows = TRUE) %>% 
    filter(!is.na(matched))
  N_pat <- nrow(matched_df) # number of patients in the matched set 
  
  # Final simulated data -------------------------------------------------------
  dropped_data %>% filter(pid %in% matched_df$pid) %>% 
    left_join(matched_df %>% select(pid, matched), by = "pid") 
}

# tempfit <- geeglm(y ~ timepoint*bav+age+male+bsa, 
#                   family = binomial('logit'),
#                   wave = factor(timepoint), 
#                   corstr = "independence", 
#                   id = pid, data = temp, 
#                   scale.fix = TRUE) %>% summary()
# 
# tempest <- tempfit$coefficients
# rownames(tempfit$coefficients)
# tempfit$coefficients$Std.err

get_model_results <- function(df, formula, corstr) {
  model <- tryCatch({
    geeglm(formula, family = binomial('logit'), wave = factor(timepoint), corstr = corstr, id = pid, data = df)
  }, error = function(e) NULL, warning = function(w) NULL)
  
  if (is.null(model)) {
    return(data.frame(term = NA, estimate = NA, std_error = NA, 
                      lower = NA, upper = NA, 
                      adj_lower = NA, adj_upper = NA, convergence = FALSE))
  }
  
  fit <- summary(model)
  est <- fit$coefficients
  z <- qnorm(1 - alpha / 2)
  lower <- est[, "Estimate"] - z * est[, "Std.err"]
  upper <- est[, "Estimate"] + z * est[, "Std.err"]
  
  p <- nrow(est)
  N <- nrow(df)
  v_cov <- vcov(model)
  V_df <- (N / (N - p)) * v_cov
  adj_se <- sqrt(diag(V_df))
  adj_lower <- est[, "Estimate"] - z * adj_se
  adj_upper <- est[, "Estimate"] + z * adj_se
  
  result <- data.frame(term = rownames(est), 
                       estimate = est[, "Estimate"], 
                       std_error = est[, "Std.err"], 
                       lower = lower, upper = upper, 
                       adj_lower = adj_lower, adj_upper = adj_upper,
                       convergence = TRUE)
  return(result)
}


set.seed(5207)
sim_df <- expand.grid(sim_id = 1:1000) %>% 
  mutate(data= map(sim_id, function(id){simulate_data(coeffs=params_sets)})) %>% 
  mutate(
    ind_mdl_full = map(data, function(df) {
      get_model_results(df, y ~ timepoint * bav + age + male + bsa, "independence")
    }),
    ind_mdl_red = map(data, function(df) {
      get_model_results(df, y ~ timepoint * bav, "independence")
    }),
    ar1_mdl_full = map(data, function(df) {
      get_model_results(df, y ~ timepoint * bav + age + male + bsa, "ar1")
    }),
    ar1_mdl_red = map(data, function(df) {
      get_model_results(df, y ~ timepoint * bav, "ar1")
    }),
    exch_mdl_full = map(data, function(df) {
      get_model_results(df, y ~ timepoint * bav + age + male + bsa, "exchangeable")
    }),
    exch_mdl_red = map(data, function(df) {
      get_model_results(df, y ~ timepoint * bav, "exchangeable")
    })
  )


# Function to clean term names
clean_terms <- function(df) {
  df %>%
    mutate(term = gsub("true", "", tolower(term), ignore.case = TRUE),
           term = gsub("(^[[:punct:]]|[[:punct:]]$)", "", term),
           term = gsub("[[:punct:]]", "_", term))
}

# Apply `clean_terms` to each data frame in `sim_df`
sim_df <- sim_df %>%
  mutate(
    ind_mdl_full = map(ind_mdl_full, clean_terms),
    ind_mdl_red = map(ind_mdl_red, clean_terms),
    ar1_mdl_full = map(ar1_mdl_full, clean_terms),
    ar1_mdl_red = map(ar1_mdl_red, clean_terms),
    exch_mdl_full = map(exch_mdl_full, clean_terms),
    exch_mdl_red = map(exch_mdl_red, clean_terms)
  )


# Calculate the convergence proportion for each model type
convergence_summary <- sim_df %>%
  summarise(
    ind_mdl_full_convergence = mean(map_lgl(ind_mdl_full, ~ .x$convergence[1])),
    ind_mdl_red_convergence = mean(map_lgl(ind_mdl_red, ~ .x$convergence[1])),
    ar1_mdl_full_convergence = mean(map_lgl(ar1_mdl_full, ~ .x$convergence[1])),
    ar1_mdl_red_convergence = mean(map_lgl(ar1_mdl_red, ~ .x$convergence[1])),
    exch_mdl_full_convergence = mean(map_lgl(exch_mdl_full, ~ .x$convergence[1])),
    exch_mdl_red_convergence = mean(map_lgl(exch_mdl_red, ~ .x$convergence[1])))

convergence_tbl <- data.frame(
  corstr = c(rep("Independence", 2), rep("AR1",2), rep("Exchangeable", 2)),
  cov_sets = rep(c("Full", "Reduced"), 3), 
  convergence = c(
    convergence_summary$ind_mdl_full_convergence,
    convergence_summary$ind_mdl_red_convergence,
    convergence_summary$ar1_mdl_full_convergence,
    convergence_summary$ar1_mdl_red_convergence,
    convergence_summary$exch_mdl_full_convergence,
    convergence_summary$exch_mdl_red_convergence
  ))

sim_res <- bind_rows(
  sim_df %>% 
    select(sim_id, ind_mdl_full) %>% 
    unnest(cols = c(ind_mdl_full)) %>% 
    pivot_wider(id_cols = sim_id,
                names_from = "term",
                values_from = c(estimate, std_error)) %>% 
    mutate(corstr = "independence", params = "Full"),
  
  sim_df %>% 
    select(sim_id, ind_mdl_red) %>% 
    unnest(cols = c(ind_mdl_red)) %>% 
    pivot_wider(id_cols = sim_id,
                names_from = "term",
                values_from = c(estimate, std_error)) %>% 
    mutate(corstr = "independence", params = "Reduced"),
  
  sim_df %>% 
    select(sim_id, ar1_mdl_full) %>% 
    unnest(cols = c(ar1_mdl_full)) %>% 
    pivot_wider(id_cols = sim_id,
                names_from = "term",
                values_from = c(estimate, std_error)) %>% 
    mutate(corstr = "ar1", params = "Full"),
  
  sim_df %>% 
    select(sim_id, ar1_mdl_red) %>% 
    unnest(cols = c(ar1_mdl_red)) %>% 
    pivot_wider(id_cols = sim_id,
                names_from = "term",
                values_from = c(estimate, std_error)) %>% 
    mutate(corstr = "ar1", params = "Reduced"),
  
  sim_df %>% 
    select(sim_id, exch_mdl_full) %>% 
    unnest(cols = c(exch_mdl_full)) %>% 
    pivot_wider(id_cols = sim_id,
                names_from = "term",
                values_from = c(estimate, std_error)) %>% 
    mutate(corstr = "exchangeable", params = "Full"), 
  
  sim_df %>% 
    select(sim_id, exch_mdl_red) %>% 
    unnest(cols = c(exch_mdl_red)) %>% 
    pivot_wider(id_cols = sim_id,
                names_from = "term",
                values_from = c(estimate, std_error)) %>% 
    mutate(corstr = "exchangeable", params = "Reduced")
)

names(sim_res)

sim_res_summary <- sim_res %>% 
  group_by(params, corstr) %>% 
  summarize_at(vars(estimate_intercept, estimate_timepoint, estimate_bav, 
                    estimate_age, estimate_male, estimate_bsa, estimate_timepoint_bav),
               list(avg = ~ mean(., na.rm = T),
                    std = ~ sd(., na.rm = T))) %>% 
  mutate(
         # bias
         intercept_bias = estimate_intercept_avg - params_sets$gee_coefs['b0'],
         time_bias      = estimate_timepoint_avg - params_sets$gee_coefs['visit'],
         bav_bias       = estimate_bav_avg - params_sets$gee_coefs['bav'],
         age_bias       = estimate_age_avg - params_sets$gee_coefs['age'],
         bsa_bias       = estimate_bsa_avg - params_sets$gee_coefs['bsa'],
         time_bav_bias  = estimate_timepoint_bav_avg - params_sets$gee_coefs['bav_visit'],
         
         # relative bias
         intercept_bias_rel = (estimate_intercept_avg - params_sets$gee_coefs['b0'])/params_sets$gee_coefs['b0'],
         time_bias_rel  = (estimate_timepoint_avg - params_sets$gee_coefs['visit'])/params_sets$gee_coefs['visit'],
         bav_bias_rel   = (estimate_bav_avg - -params_sets$gee_coefs['bav'])/abs(-params_sets$gee_coefs['bav']),
         age_bias_rel   = (estimate_age_avg - params_sets$gee_coefs['age'])/params_sets$gee_coefs['age'],
         bsa_bias_red   = (estimate_bsa_avg - params_sets$gee_coefs['bsa'])/params_sets$gee_coefs['bsa'],
         time_bav_bias_rel = (estimate_timepoint_bav_avg - params_sets$gee_coefs['bav_visit'])/params_sets$gee_coefs['bav_visit'],
         
         # mean squared error
         intercept_mse = estimate_intercept_std^2 + intercept_bias^2,
         time_mse      = estimate_timepoint_std^2 + time_bias^2,
         bav_mse       = estimate_bav_std^2 + bav_bias^2,
         age_mse       = estimate_age_std^2 + age_bias^2,
         bsa_mse       = estimate_bsa_std^2 + bsa_bias^2,
         time_bav_mse = estimate_timepoint_bav_std^2 + time_bav_bias^2)

sim_res_summary

