# title: Cardiovascular Data Simulation 
# author: Jinyu Luo
# version: 2024-06-26
rm(list = ls())
# Required R packages ----------------------------------------------------------
library(CorBin)
library(optmatch)
library(tidyverse)
library(geepack)
library(lme4)
library(parallel)
library(survival)
library(simsurv)
library(kableExtra)
library(R.utils)
load("data/realcoefs.RData")

# Global Variables -------------------------------------------------------------
n_patients <- 250 
n_sim <- 1000
maxT <- 6
rho <- 0.3
alpha_ci <- 0.05

true_coefs <- full_ar1$coefficients$Estimate
names(true_coefs) <- c("g0", "bav", "visit", "age", "male", "bsa", "bav_visit")
corr_alpha <- full_ar1$corr$Estimate
surv_coefs <- surv_coefs[,-4]
rownames(surv_coefs)[4] <- "male"


# GEE Helper Functions ---------------------------------------------------------
# 1. Function for creating profile for one simulation 
patient_profile <- function(){
  
  full_data <- NULL 
  
  for (id in 1:n_patients){
    # Simulate covariates
    age = rnorm(1, mean = 50, sd = 10)
    male = rbinom(1, size = 1, prob = 0.65)
    bsa = rnorm(1, mean = 2, sd = 0.3)
    logit_bav <- ps_coefs[1] + ps_coefs[2]*age + ps_coefs[3]*male + 
      ps_coefs[4]*bsa
    prob_bav <- exp(logit_bav) / (1+exp(logit_bav))
    bav = rbinom(1, size = 1, prob = prob_bav)
    
    # Simulate binary outcome 
    vst <- 1:maxT
    logit_y <- true_coefs[1] + true_coefs[2]*bav + true_coefs[3]*vst +
      true_coefs[4]*age + true_coefs[5]*male + true_coefs[6]*bsa + 
      true_coefs[7]*bav*vst
    prob_y <- exp(logit_y) / (1+exp(logit_y)) 
    y <- t(cBern(n = 1, p = prob_y, rho = rho, type = "DCP"))
    
    full_data <- rbind(full_data, 
                       data.frame(id = rep(id, each = maxT), 
                                  age = rep(age, each = maxT), 
                                  male = rep(male, each = maxT), 
                                  bsa = rep(bsa, each = maxT), 
                                  bav = rep(bav, each = maxT), 
                                  visit = vst) %>% cbind(y)
                       )
  }

  full_data
}

# 2. Function to simulate dropouts 
sim_dropouts <- function(dat){
  
  temp <- dat %>% relocate(c(y, bav), .after = id) 
  
  base_df <- temp %>% filter(visit == 1) %>% select(-visit)
  
  # (1) Drop visits for patients who had only one follow up 
  pat_set1 <- simsurv(lambdas = exp(surv_coefs$fit1[6]), # scale parameter
                      gammas = exp(surv_coefs$fit1[7]), # shape parameter
                      x = base_df, maxt = 1, 
                      betas = t(surv_coefs)[1,1:5],
                      dist = "weibull") %>% filter(status == 1)  %>% pull(id)
  temp <- temp %>% 
    filter(!((id %in% pat_set1) & visit > 2)) %>% 
    group_by(id) %>% 
    mutate(total_visit = n()) %>% 
    ungroup()
  
  
  # (2) Drop visits for patients who had two follow up visits, i.e., total_visit == 3 
  Z2 <- temp %>% filter(total_visit > 2) %>% 
    group_by(id) %>% slice(1) %>% 
    select(-c(visit, total_visit))
  pat_set2 <- simsurv(lambdas = exp(surv_coefs$fit2[6]), # scale parameter
                      gammas = exp(surv_coefs$fit2[7]), # shape parameter
                      x = Z2, maxt = 1, 
                      betas = t(surv_coefs)[2,1:5],
                      dist = "weibull") %>% filter(status == 1)  %>% pull(id)
  temp <- temp %>% 
    filter(!((id %in% pat_set2) & visit > 3)) %>% 
    group_by(id) %>% 
    mutate(total_visit = n()) %>% 
    ungroup()
  
  # (3) Drop visits for patients who had three follow up visits, i.e., total_visit == 4 
  Z3 <- temp %>% group_by(id) %>% slice(1) %>% 
    filter(total_visit > 3) %>% select(-c(visit, total_visit))
  pat_set3 <- simsurv(lambdas = exp(surv_coefs$fit3[6]), 
                      gammas = exp(surv_coefs$fit3[7]), # shape parameter
                      x = Z3, maxt = 1,
                      betas = t(surv_coefs)[3,1:5],
                      dist = "weibull") %>% filter(status == 1)  %>% pull(id)
  temp <- temp %>% 
    filter(!((id %in% pat_set3) & visit > 4)) %>% 
    group_by(id) %>% 
    mutate(total_visit = n()) %>% 
    ungroup()
  
  # (4) Drop visits for patients who had four follow up visits, i.e., total_visit == 5 
  Z4 <- temp %>% group_by(id) %>% slice(1)%>% 
    filter(total_visit == 6)  %>% select(-c(visit, total_visit))
  pat_set4 <- simsurv(lambdas = exp(surv_coefs$fit5[6]), 
                      gammas = exp(surv_coefs$fit5[7]), # shape parameter
                      x = Z4, maxt = 1,
                      betas = t(surv_coefs)[4,1:5],
                      dist = "weibull") %>% filter(status == 1)  %>% pull(id)
  temp <- temp %>% 
    filter(!((id %in% pat_set4) & visit > 5)) %>% 
    # filter(!((id %in% Z4$id[pat_set4]) & visit > 5)) %>% 
    group_by(id) %>% 
    mutate(total_visit = n()) %>% 
    ungroup()
  
  temp 
}

# 3. Function for Propensity score matching 
ps_match <- function(dat){
  base_info <- dat %>% filter(visit == 1) %>% select(-total_visit)
  ps_model <- glm(bav ~ age + male + bsa, family = binomial, data = base_info)
  pps_match <- pairmatch(ps_model, data = base_info)
  matched_df <- data.frame(base_info, matched = pps_match, check.rows = TRUE) %>% 
    filter(!is.na(matched))
  matchid <- matched_df %>% select(id, matched)
  finaldata <- dat %>% right_join(matchid, by = "id")
  finaldata
}

# 4. Function for GEE fit 
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

# 5. Function to calculate coverage for each simulation 
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

# 6. Function to calculate mean estimates for each term
calculate_mean_estimates <- function(results_list) {
  results_list %>%
    map_dfr(~.x %>% 
              select(term, estimate)) %>%
    group_by(term) %>%
    summarize(mean_estimate = mean(estimate, na.rm = TRUE))
}

# 7. Function to calculate mean rho and rho_se for AR1 and Exchangeable 
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

# Data Simulation --------------------------------------------------------------
set.seed(5207)

sim_df <- expand.grid(sim_id = 1:n_sim) %>%
  mutate(full_data= map(sim_id, function(id){patient_profile()})) %>%
  mutate(dropout_data = map(full_data, ~sim_dropouts(.x))) %>%
  mutate(matched_data = map(dropout_data, ~ps_match(.x)))

matched_df <- sim_df %>% pull(matched_data)

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
# [1] 7

sim_results <- gee_fits_df %>% filter(!(sim_id %in% non_converged$sim_id))
nrow(sim_results)
# [1] 993

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
  filter(if_any(-sim_id, ~ . == TRUE))

extreme_sim %>% select(-sim_id) %>% colMeans()
# ind_full_se  ar1_full_se exch_full_se   ind_red_se   ar1_red_se   exch_red_se  
#     0.66667      0.58333      0.80556      0.00000      0.02778       0.00000

sim_results <- sim_results %>% filter(!(sim_id %in% extreme_sim$sim_id))
nrow(sim_results)
# [1] 957


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




# QLS Helper Functions ---------------------------------------------------------
alpha_stg1_ar1 <- function(mdat, Z, Qinv){
  Fa <- Fb <- 0 
  for (i in mdat$clusterID) { # for each pair
    S1_j <- S2_j <- S1_ja <- S1_jb <- 0
    
    t_i1 <- nlevels(as.factor(mdat[mdat$clusterID == i & mdat$cluster.var==1,]$visit))
    t_i2 <- nlevels(as.factor(mdat[mdat$clusterID == i & mdat$cluster.var==2,]$visit))
    t_ij <- min(c(t_i1, t_i2))
    Z_i1 <- Z[mdat$clusterID == i & mdat$cluster.var == 1]
    Z_i2 <- Z[mdat$clusterID == i & mdat$cluster.var == 2]

    # Check if the lengths match t_ij before creating the matrices
    if (length(Z_i1) >= t_ij && length(Z_i2) >= t_ij) {
      matZ_i1 <- matrix(Z_i1[1:t_ij], nrow = t_ij)
      matZ_i2 <- matrix(Z_i2[1:t_ij], nrow = t_ij)
      
      if (t_ij > 1) {
        for (k in 1:(t_ij-1)) {
          matZ1 <- matrix(c(matZ_i1[k], matZ_i2[k]), nrow = 2)
          matZ2 <- matrix(c(matZ_i1[k+1], matZ_i2[k+1]), nrow = 2)
          S2_j <- S2_j + t(matZ1) %*% Qinv %*% matZ2
        }
        if (t_ij == 2) {
          for (k in 1:t_ij) {
            matZ <- matrix(c(matZ_i1[k], matZ_i2[k]), nrow = 2)
            S1_j <- S1_j + t(matZ) %*% Qinv %*% matZ
          }
        } else {
          for (k in 1:t_ij) {
            matZ <- matrix(c(matZ_i1[k], matZ_i2[k]), nrow = 2)
            S1_ja <- S1_ja + t(matZ) %*% Qinv %*% matZ
          }
          for (k in 2:(t_ij-1)) {
            matZ <- matrix(c(matZ_i1[k], matZ_i2[k]), nrow = 2)
            S1_jb <- S1_jb + t(matZ) %*% Qinv %*% matZ
          }
          S1_j <- S1_ja + S1_jb
        }
      }
      
      Fa <- Fa + S1_j
      Fb <- Fb + S2_j
    } else {
      warning(paste("Cluster", i, "has inconsistent lengths for Z_i1 and Z_i2 with t_ij =", t_ij))
    }
  }
    
  ### stage 1 estimate of alpha
  var_discriminant <- (Fa - 2 * Fb) * (Fa + 2 * Fb)
  if (var_discriminant < 0) {
    warning("Quasi-variance discriminant is negative. Setting alpha0 to NA.")
    alpha0 <- NA
  } else {
    alpha0 <- (Fa - sqrt(var_discriminant)) / (2 * Fb)
  }
  return(alpha0)
}


alpha_stg2_ar1 <- function(alpha0){
  alpha <- as.numeric( 2 * alpha0 / ( 1 + alpha0 ^ 2 ) )
  return(alpha)
}


estalpha1_exch <- function(mdat, Z, Qinv){
  match.call()
  alphafun <- function(alpha){
    GG1 <- GG2 <- 0
    for (i in unique(mdat$clusterID)){
      GG1j <- GG2j <- 0
      
      t_i1 <- nlevels(as.factor(mdat[mdat$clusterID == i & mdat$cluster.var==1,]$visit))
      t_i2 <- nlevels(as.factor(mdat[mdat$clusterID == i & mdat$cluster.var==2,]$visit))
      t_ij <- max(c(t_i1, t_i2))
      Z_i1 <- Z[mdat$clusterID == i & mdat$cluster.var == 1]
      if (t_i1 < t_ij) {Z_i1 <- c(Z_i1, rep(0, t_ij - t_i1))}
      matZ_i1 <- matrix(Z_i1, nrow = t_ij) 
      Z_i2 <- Z[mdat$clusterID == i & mdat$cluster.var == 2]
      if (t_i2 < t_ij) {Z_i2 <- c(Z_i2, rep(0, t_ij - t_i2))} 
      matZ_i2 <- matrix(Z_i2, nrow = t_ij) 
      
      #if(t_ij > 1){
      g1 <- vector()
      for(t in 1:t_ij){
        matZ <- matrix(c(matZ_i1[t],matZ_i2[t]),nrow=2)
        g1[t] <- t(matZ) %*% Qinv %*% matZ
      } 
      G1 <- sum(g1)
      
      g2 <- vector()
      G2 <- 0 #
      if(t_ij > 1){ #
        for(t in 1:(t_ij - 1)){
          for(tt in (t+1):t_ij){
            matZ1 <- matrix(c(matZ_i1[t],matZ_i2[t]),nrow=2)
            matZ2 <- matrix(c(matZ_i1[tt],matZ_i2[tt]),nrow=2)
            g2 <- c(g2, t(matZ1) %*% Qinv %*% matZ2) 
          }
        }
        G2 <- sum(g2)
        
      }
      denom <- ( 1 + ( t_ij - 1 ) * alpha ) ^ 2
      num1 <- alpha ^ 2 * ( t_ij - 1 ) * ( t_ij - 2 ) + 2 * alpha * ( t_ij - 1 )
      num2 <- ( 1 + alpha ^ 2 * ( t_ij - 1 ) )
      
      GG1j <- GG1j + ( G1 * num1 ) / denom
      GG2j <- GG2j + ( G2 * num2 ) / denom
      
      GG1 <- GG1 + GG1j
      GG2 <- GG2 + GG2j
    }
    GG1 - 2 * GG2
  }
  ### stage 1 estimate of alpha
  alpha0 <- uniroot(alphafun, c(0,1), tol = 1e-10, extendInt = "yes")$root
  return(alpha0)
}

estalpha2_exch <- function(alpha0, mdat){
  match.call()
  alphapart1 <- alphapart2 <- 0
  
  for (i in mdat$clusterID){
    alphapart1j <- alphapart2j <- 0
    
    t_ij <- 5 #nlevels(as.factor(mdat[mdat$cluster_id == i & mdat$cluster.var == j,]$visit))
    if(t_ij > 1){
      alphapart1num <- alpha0 * ( t_ij - 1 )* ( alpha0 * (t_ij - 2) + 2 )
      alphapart2num <- ( t_ij - 1 ) * ( 1 + alpha0 ^ 2 * (t_ij - 1) )
      alphaden <- ( 1 + alpha0 * ( t_ij - 1 ) ) ^ 2
      
      alphapart1j <- alphapart1j + alphapart1num / alphaden
      alphapart2j <- alphapart2j + alphapart2num / alphaden
    }
    alphapart1 <- alphapart1 + alphapart1j
    alphapart2 <- alphapart2 + alphapart2j
  }
  alpha <- alphapart1 / alphapart2
  return(alpha)
}

tau_stg1 <- function(mdat, maxT, Z, corstr, alpha0){
  Fa <- Fb <- 0
  for (i in mdat$clusterID){
    a_1 <- a_2 <- 0
    t_i1 <- nlevels(as.factor(mdat[mdat$clusterID == i & mdat$cluster.var==1,]$visit))
    t_i2 <- nlevels(as.factor(mdat[mdat$clusterID == i & mdat$cluster.var==2,]$visit))
    if (corstr == "independence") {Rinv1 <- solve(diag(t_i1))} 
    if (corstr == "ar1") {Rinv1 <- solve(ar1_cor(alpha0, t_i1))} 
    if (corstr == "exchangeable") {Rinv1 <- solve(exch_cormat(alpha0, t_i1))}
    if (corstr == "independence") {Rinv2 <- solve(diag(t_i2))} 
    if (corstr == "ar1") {Rinv2 <- solve(ar1_cor(alpha0, t_i2))} 
    if (corstr == "exchangeable") {Rinv2 <- solve(exch_cormat(alpha0, t_i2))}
    Rinv <- solve(exch_cormat(alpha0, maxT))
    
    Z_i1 <- Z[mdat$clusterID == i & mdat$cluster.var == 1]
    matZ_i1 <- matrix(Z_i1, nrow = t_i1)
    Z_i2 <- Z[mdat$clusterID == i & mdat$cluster.var == 2]
    matZ_i2 <- matrix(Z_i2, nrow = t_i2)
    a_1 <- a_1 + t(matZ_i1) %*% Rinv1 %*% matZ_i1 + t(matZ_i2) %*% Rinv2 %*% matZ_i2
    
    if (maxT > t_i1) {matZ_i1 <- c(matZ_i1, rep(0, maxT - t_i1))}
    if (maxT > t_i2) {matZ_i2 <- c(matZ_i2, rep(0, maxT - t_i2))}
    a_2 <- a_2 + t(matZ_i1) %*% Rinv %*% matZ_i2
    
    Fa <- Fa + a_1
    Fb <- Fb + a_2
    
  }
  ### stage 1 estimate of tau
  tau0 <- ( Fa - sqrt( ( Fa - 2 * Fb ) * ( Fa + 2 * Fb ) ) ) / ( 2 * Fb )
  return(tau0)
}

tau_stg2 <- function(tau0){
  tau <- as.numeric( 2 * tau0 / ( 1 + tau0 ^ 2 ) )
  return(tau)
}

exch_cormat <- function(rho, n) {
  # Create an nxn matrix filled with tau
  cor_matrix <- matrix(rho, n, n)
  # Set the diagonal to 1
  diag(cor_matrix) <- 1
  return(cor_matrix)
}

ar1_cor <- function(rho,n) {
  rho <- as.numeric(rho)
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

Sigma <- function(data,tau, alpha, corstr, time.var){
  Sigma_list <- list()
  for (i in unique(data$clusterID)) {
    Qi <- exch_cormat(tau, 2)
    ti1 <- nlevels(as.factor(data[data$clusterID == i & data$cluster.var==1,]$visit))
    ti2 <- nlevels(as.factor(data[data$clusterID == i & data$cluster.var==2,]$visit))
    ni <- ti1+ti2
    if (ti1 >= ti2)
    {
      if (corstr == "independence") {Ri <- diag(ti1)}
      if (corstr == "ar1") {Ri <- ar1_cor(alpha, ti1)}
      if (corstr == "exchangeable") {Ri <- exch_cormat(alpha, ti1)}
      Fi <- kronecker(Qi,Ri)
      Sigma_i <- Fi[1:ni, 1:ni]
      Sigma_list <- c(Sigma_list, list(Sigma_i))
    }
    else {
      if (corstr == "independence") {Ri <- diag(ti2)}
      if (corstr == "ar1") {Ri <- ar1_cor(alpha, ti2)}
      if (corstr == "exchangeable") {Ri <- exch_cormat(alpha, ti2)}
      Fi <- kronecker(Qi,Ri)
      Sigma_i <- Fi[1:ni, 1:ni]
      Sigma_list <- c(Sigma_list, list(Sigma_i))
    }
  }
  Sigma_list
}

beta_hat <- function(formula,data, time.var, corstr, tau, alpha) {
  X <- model.matrix(object=formula, data = data) #design matrix
  y <- as.matrix(data$y) #response variable
  Sigma_list <- list()
  Xt_Sigma_inv_X <- list()
  Xt_Sigma_inv_y <- list()
  S <- Sigma(data=data, tau=tau, alpha=alpha, corstr=corstr, time.var=time.var)
  for (i in 1:length(S)) {
    ti1 <- nlevels(as.factor(data[data$clusterID == i & data$cluster.var==1,]$visit))
    ti2 <- nlevels(as.factor(data[data$clusterID == i & data$cluster.var==2,]$visit))
    if (ti1 >= ti2){
      Xi <- rbind(X[data$clusterID==i & data$cluster.var==1,], 
                  X[data$clusterID==i & data$cluster.var==2,])
      yi <- rbind(as.matrix(y[data$clusterID==i & data$cluster.var==1,]), 
                  as.matrix(y[data$clusterID==i & data$cluster.var==2]))
    }
    else {
      Xi <- rbind(X[data$clusterID==i & data$cluster.var==2,], 
                  X[data$clusterID==i & data$cluster.var==1,])
      yi <- rbind(as.matrix(y[data$clusterID==i & data$cluster.var==2,]), 
                  as.matrix(y[data$clusterID==i & data$cluster.var==1]))
    }
    Sigma_inv <- solve(S[[i]])
    Xt_Sigma_inv_X_i <- t(Xi) %*% Sigma_inv %*% Xi
    Xt_Sigma_inv_X[[i]] <- Xt_Sigma_inv_X_i
    Xt_Sigma_inv_y_i <- t(Xi) %*% Sigma_inv %*% yi
    Xt_Sigma_inv_y[[i]] <- Xt_Sigma_inv_y_i
  }
  return(solve(Reduce("+", Xt_Sigma_inv_X)) %*% Reduce("+",Xt_Sigma_inv_y))
}


sandwich <- function(formula,data,beta_hat,alpha, corstr){
  X <- model.matrix(object=formula, data = data)
  y <- as.matrix(data$y)
  W <- list()
  mid <- list()
  for (i in unique(data$clusterID)) {
    Xi <- X[data$clusterID==i,]
    yi <- y[data$clusterID==i]
    Zi <- yi - (Xi %*% as.matrix(beta_hat))
    ti1 <- nlevels(as.factor(data[data$clusterID == i & data$cluster.var==1,]$visit))
    ti2 <- nlevels(as.factor(data[data$clusterID == i & data$cluster.var==2,]$visit))
    ni <- ti1 + ti2
    Ai <- diag(ni) ^ (1/2)
    if (corstr == "independence") {Ri <- diag(ni)}
    if (corstr == "ar1") {Ri <- ar1_cor(alpha, ni)}
    if (corstr == "exchangeable") {Ri <- exch_cormat(alpha, ni)}
    Wi <- t(Xi) %*% Ai %*% solve(Ri) %*% Ai %*% Xi
    W[[i]] <- Wi
    mid_i <- t(Xi) %*% Ai %*% solve(Ri) %*% Zi %*% t(Zi) %*% solve(Ri) %*% Ai %*% Xi
    mid[[i]] <- mid_i
  }
  Wn_inv <- solve(Reduce("+", W))
  mid_n <- Reduce("+", mid)
  out <- list()
  out$vcov <- Wn_inv %*% mid_n %*% Wn_inv
  out$se <- sqrt(diag(out$vcov))
  return(out)
}


qls <- function(formula, data, corstr, maxT, time.var){
  iter <- 0
  bdiff <- c(1,1,1,1)
  alpha0  <- 0.1 # initial alpha estimate
  
  # use independent GEE to get initial beta estimates 
  init_mod <- geeglm(formula, data = data, family = binomial('logit'), 
                     id = id, waves = factor(visit), 
                     corstr = "independence", scale.fix = TRUE) 
  #summary(init_mod)
  beta0 <- as.vector(coef(init_mod))
  Z0 <- residuals(init_mod,"pearson") #init_mod$residuals Z0[1:5][1,]
  
  # compute initial tau estimate
  tau0 <- tau_stg1(mdat=data, maxT=maxT, Z = Z0, corstr = corstr,alpha0=alpha0) 
  
  while(max(abs(bdiff)) > .00000001){
    betahat <- beta_hat(formula=formula,data=data, time.var=time.var, 
                        corstr=corstr, tau=tau0, alpha=alpha0)
    beta1 <- as.vector(betahat)
    if (all(!is.na(betahat))){bdiff <- beta1 - beta0} #***
    
    # update tau0
    Z1 <- as.matrix(df$y) - 
      model.matrix(object=formula, data = data) %*% as.matrix(betahat)
    
    tau00 <- tau_stg1(mdat=data, maxT=maxT, Z = Z1, corstr = corstr,alpha0=alpha0) 
    # update alpha0 (initial alpha0 for the next iteration)
    if (!is.na(tau00)) {tau0 <- tau00}
    #print(tau0)
    
    Qinv <- solve(exch_cormat(tau0, 2))
    if (corstr == "independence") {alpha0 <- 0}
    if (corstr == "ar1") {alpha0 <- alpha_stg1_ar1(mdat=data, Z=Z1, Qinv=Qinv)}
    if (corstr == "exchangeable") {alpha0 <- estalpha1_exch(mdat=data, Z=Z1, Qinv=Qinv)}
    
    iter <- iter + 1
    beta0 <- beta1
    # print(paste("iter:", iter, sep = " "))
    # print(paste("alpha0:",alpha0, sep = " "))
    # print(paste("tau0:",as.numeric(tau0), sep = " "))
    # print(paste("bdiff:",max(abs(bdiff)), sep = " "))
  }
  
  # after converge, get stage 2 estimates
  tau2 <- tau_stg2(tau0)
  if (corstr == "independence") {alpha2 <- alpha0}
  if (corstr == "ar1") {alpha2 <- alpha_stg2_ar1(alpha0)}
  if (corstr == "exchangeable") {alpha2 <- estalpha2_exch(alpha0, mdat = data)}
  
  betahat1 <- beta_hat(formula=formula,data=data, time.var=time.var, 
                       corstr=corstr, tau=tau2, alpha=alpha2)
  beta <- as.vector(betahat1)
  sandwich_out <- sandwich(formula = formula, data = data, beta_hat = betahat1, 
                           alpha = alpha2, corstr = corstr)
  se <- sandwich_out$se
  vcov <- sandwich_out$vcov
  
  fit <- list()
  fit$call <- match.call()
  fit$coefficients <- beta
  fit$se <- se
  fit$alpha <- alpha2
  fit$tau <- tau2
  fit$niter <- iter
  fit$vcov <- vcov
  fit
}

qls_ci <- function(model, level = 0.95) {
  # Calculate the lower and upper bounds of the confidence interval for each parameter
  lower <- model$coefficients - qnorm(1-(1-0.95)/2) * model$se
  upper <- model$coefficients + qnorm(1-(1-0.95)/2) * model$se
  results <- data.frame(lower = lower, upper = upper)
  return(results)
}

adj_qls_ci <- function(model, N, level = 0.95) {
  # Calculate the adjusted standard errors using DF-corrected sandwich estimator
  adj_se <- sqrt(diag((N/(N-p))*model$vcov))
  # Calculate the lower and upper bounds of the confidence interval for each parameter
  lower <- model$coefficients - qnorm(1-(1-0.95)/2) * adj_se
  upper <- model$coefficients + qnorm(1-(1-0.95)/2) * adj_se
  results <- data.frame(lower = lower, upper = upper)
  return(results)
}


get_qls_results <- function(df, formula, corstr, type) {
  mdl_name <- paste(type, corstr, sep = " ")
  print(paste("Starting model fitting with correlation structure:", mdl_name))
  
  model <- tryCatch({
    withTimeout({
      qls(formula, data = df, time.var = df$visit, maxT = 6, corstr = corstr)
    }, timeout = 60, onTimeout = "warning")
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
  
  result <- data.frame(term = names(model$se), 
                       estimate = model$coefficients, 
                       std_error = model$se, 
                       convergence = TRUE, 
                       rho = model$alpha,
                       tau = model$tau) %>% cbind(model$vcov)
  
  return(result)
}


# QLS Fit ----------------------------------------------------------------------
# QLS data 
qls_df <- sim_df %>% select(sim_id, matched_data) %>% 
  mutate(matched_data = map(matched_data, ~ .x %>% 
                              arrange(matched) %>% 
                              mutate(clusterID = as.integer(factor(matched))) %>% 
                              select(-matched) %>% 
                              group_by(clusterID) %>% 
                              mutate(cluster.var = ifelse(bav == 0, 1, 2), 
                                     order = row_number()) %>% 
                              relocate(c(clusterID, cluster.var, order), .after = id)), 
         n_pairs = map_int(matched_data, ~ n_distinct(.x$clusterID)))


qls_mdls <- list(
  ind_mdl_full = list(formula = y ~ visit * bav + age + male + bsa, 
                      corstr = "independence", type = "Adjusted"),
  ind_mdl_red = list(formula = y ~ visit * bav, 
                     corstr = "independence", type = "Unadjusted"),
  ar1_mdl_full = list(formula = y ~ visit * bav + age + male + bsa, 
                      corstr = "ar1", type = "Adjusted"),
  ar1_mdl_red = list(formula = y ~ visit * bav, 
                     corstr = "ar1", type = "Unadjusted"),
  exch_mdl_full = list(formula = y ~ visit * bav + age + male + bsa, 
                       corstr = "exchangeable",  type = "Adjusted"),
  exch_mdl_red = list(formula = y ~ visit * bav, 
                      corstr = "exchangeable", type = "Unadjusted")
)


# Initialize a list to store the results
qls_results <- list()

# Loop through each dataset and fit all models
for (i in seq_len(nrow(qls_df))) {
  df <- qls_df$matched_data[[i]]
  sim_id <- qls_df$sim_id[i]
  
  model_results <- list()
  
  for (mdl in names(qls_mdls)) {
    spec <- qls_mdls[[mdl]]
    model_results[[mdl]] <- get_qls_results(df, spec$formula, spec$corstr, spec$type)
  }
  
  qls_results[[i]] <- c(list(sim_id = sim_id), model_results)
  
  print(paste("Completed simulation", i, "out of", nrow(qls_df)))
}


qls_fits_df <- tibble(
  sim_id = map(qls_results, "sim_id"),
  ind_mdl_full = map(qls_results, "ind_mdl_full"), 
  ind_mdl_red = map(qls_results, "ind_mdl_red"), 
  ar1_mdl_full = map(qls_results, "ar1_mdl_full"), 
  ar1_mdl_red = map(qls_results, "ar1_mdl_red"), 
  exch_mdl_full = map(qls_results, "exch_mdl_full"), 
  exch_mdl_red = map(qls_results, "exch_mdl_red")
)

## QLS Convergence ------------------------------------------------------------
qls_convergence <- qls_fits_df %>%
  mutate(ind_full = map_lgl(ind_mdl_full, ~any(.x$convergence == TRUE)),
         ar1_full = map_lgl(ar1_mdl_full, ~any(.x$convergence == TRUE)),
         exch_full = map_lgl(exch_mdl_full, ~any(.x$convergence == TRUE)),
         ind_red = map_lgl(ind_mdl_red, ~any(.x$convergence == TRUE)), 
         ar1_red = map_lgl(ar1_mdl_red, ~any(.x$convergence == TRUE)), 
         exch_red = map_lgl(exch_mdl_red, ~any(.x$convergence == TRUE))) %>% 
  select(sim_id, ind_full:exch_red) 

qls_convergence_res <- qls_convergence %>% select(-sim_id) %>% colMeans()

# Filter out non-converged simulations
qls_non_converged <- qls_convergence %>% 
  filter(if_all(ind_full:exch_red, ~ . == FALSE)) %>% 
  select(sim_id, ind_full:exch_red) 
nrow(qls_non_converged)

qls_sim_results <- qls_fits_df %>% filter(!(sim_id %in% qls_non_converged$sim_id))
nrow(sim_results)
# [1] 993

# Remove Simulations with SE > 5 
qls_extreme_sim <- qls_sim_results %>% 
  mutate(ind_full_se = map_lgl(ind_mdl_full, ~any(.x$std_error > se_max)),
         ar1_full_se = map_lgl(ar1_mdl_full, ~any(.x$std_error > se_max)),
         exch_full_se = map_lgl(exch_mdl_full, ~any(.x$std_error > se_max)),
         ind_red_se = map_lgl(ind_mdl_red, ~any(.x$std_error > se_max)), 
         ar1_red_se = map_lgl(ar1_mdl_red, ~any(.x$std_error > se_max)), 
         exch_red_se = map_lgl(exch_mdl_red, ~any(.x$std_error > se_max))) %>% 
  select(sim_id, ind_full_se:exch_red_se) %>% 
  filter(if_any(-sim_id, ~ . == TRUE))

qls_extreme_sim %>% select(-sim_id) %>% colMeans(na.rm = TRUE) 

qls_sim_results <- qls_sim_results %>% filter(!(sim_id %in% qls_extreme_sim$sim_id))
nrow(sim_results)




# 
# qls_fit_results <- qls_df %>% 
#   mutate(
#     ind_mdl_full = map(matched_data, function(df) {
#       get_qls_results(df, y ~ visit * bav + age + male + bsa, "ind")
#     }),
#     ind_mdl_red = map(matched_data, function(df) {
#       get_qls_results(df, y ~ visit * bav, "ind")
#     }),
#     ar1_mdl_full = map(matched_data, function(df) {
#       get_qls_results(df, y ~ visit * bav + age + male + bsa, "ar1")
#     }),
#     ar1_mdl_red = map(matched_data, function(df) {
#       get_qls_results(df, y ~ visit * bav, "ar1")
#     }),
#     exch_mdl_full = map(matched_data, function(df) {
#       get_qls_results(df, y ~ visit * bav + age + male + bsa, "exch")
#     }),
#     exch_mdl_red = map(matched_data, function(df) {
#       get_qls_results(df, y ~ visit * bav, "exch")
#     }))

# save(sim_results, convergence, non_converged, extreme_sim, 
#      unadjusted_coverage, adjusted_coverage, 
#      mean_est_adj, mean_est_unadj, true_set,
#      matched_info, mean_rho_values, qls_df,
#      file = "outputs/sim_results.RData")

