
library(CorBin)
library(optmatch)
library(tidyverse)
library(geepack)
library(lme4)
library(parallel)
library(survival)
library(simsurv)
library(kableExtra)
library(survival)

# Initialize Global Variables
simnum <- 1050 # 1050 # the total number of simulation
N <- 500       # the total number of patients for each simulation 
maxT <- 5      # maximum number of visit 
visits <- rep(1:maxT, N)     # vector of visits for all patients
pid <- rep(1:N, each = maxT) # vector of patient IDs
rho <- 0.3

# Coefficients for exposure group simulation 
a_bav <- -0.4  # intercept 
b1_bav <- -0.1 # age 
b2_bav <- 1.2  # male
b3_bav <- 3    # BSA

# Coefficients for outcome simulation
a_oc <- -1     # intercept 
b1_oc <- -0.05 # visit 
b2_oc <- -2    # BAV 
b3_oc <- -0.05 # age 
b4_oc <- 1.4   # male 
b5_oc <- 0.5   # BSA
b6_oc <- 0.5   # BAV * visit

true_par_full <- c(a_oc, b1_oc, b2_oc, b3_oc, b4_oc, b5_oc, b6_oc)
par_name_full <- c("Intercept", "age", "male", "bsa", "bav", "visit", "bav:visit")
num_coeffs_full <- length(true_par_full)

true_par_red <- c(true_par_full[1:3], b6_oc)
par_name_red <- c("Intercept", "bav", "visit", "bav:visit")
num_coeffs_red <- length(true_par_red)

# Survival variables 
beta_age <- -0.03810  # Effect of age on hazard rate
beta_BAV <- 0.16574  # Effect of BAV on hazard rate
beta_y <- 0.93280
beta_male <- -0.67235
beta_bsa <- -0.52135
# wb_age <- -0.019
# wb_bav <- 0.159
# wb_male <- -0.462
# wb_bsa <- -0.456
# wb_y <- 0.654



# Function to calculate CIs with unadjusted and DF-adjusted sandwich estimator
geeci <- function(model, N, alpha = 0.05) {
  fit <- summary(model)
  est <- fit$coefficients[, "Estimate"]
  se <- fit$coefficients[, "Std.err"] 
  
  p <- length(est) # Number of parameters
  v_cov <- vcov(model) # Sandwich estimator
  V_df <- (N / (N - p)) * v_cov # DF correction
  adj_se <- sqrt(diag(V_df)) # Standard errors from corrected matrix
  
  z <- qnorm(1 - alpha/2)
  lower <- est - z * se; adj_lower <- est - z*adj_se
  upper <- est + z * se; adj_upper <- est + z*adj_se
  
  result <- data.frame(estimate = est, se = se, 
                       lower = lower, upper = upper, 
                       adj_lower = adj_lower, adj_upper = adj_upper)
  return(result)
}

# Function for Fitting GEE model and Extract Estimates 
gee.output <- function(type = "full", dat, corstr, N) {
  if (type == "full"){
    true_p = true_par_full
    fmod = y ~ visit + bav + age + male + bsa + bav:visit
    fit <- geeglm(fmod, family = binomial("logit"), id = dat$id, data = dat, 
                  corstr = corstr, scale.fix = TRUE)
  } else{
    true_p = true_par_red
    rmod = y ~ bav+visit+bav:visit
    fit <- geeglm(rmod, family = binomial("logit"), id = dat$id, data = dat, 
                  corstr = corstr, scale.fix = TRUE)
  }
  
  est <- geeci(fit, N)
  cov_prob <- c(ifelse((est$lower <= true_p & est$upper >= true_p), 1, 0),
                ifelse((est$adj_lower <= true_p & est$adj_upper >= true_p), 1, 0))
  est_set <- est[,1:2] %>% mutate(bias = (estimate - true_p)/true_p, mse = (estimate - true_p)^2)
  out <- list("Estimates" = est_set, "Cov.Prob" = cov_prob)
  
  return(out)
}

corstrs <- c("independence", "exchangeable", "ar1")
mod_types <- c("full", "reduced")

est_res <- list(); se_res <- list(); bias_res <- list(); mse_res <- list(); cov_prob <- list()
est_do <- list(); se_do <- list(); bias_do <- list(); mse_do <- list(); Cprob_do <- list()

for (i in 1:2) {
  p1 <- ifelse(i == 2, 4, 7)
  p2 <- ifelse(i == 2, 8, 14)
  est_res[[i]] <- list("independence" = matrix(NA, simnum, p1), 
                       "exchangeable" = matrix(NA, simnum, p1), 
                       "ar1" = matrix(NA, simnum, p1))
  est_do[[i]] <- list("independence" = matrix(NA, simnum, p1), 
                      "exchangeable" = matrix(NA, simnum, p1), 
                      "ar1" = matrix(NA, simnum, p1))
  se_res[[i]] <- list("independence" = matrix(NA, simnum, p1), 
                      "exchangeable" = matrix(NA, simnum, p1), 
                      "AR1" = matrix(NA, simnum, p1))
  se_do[[i]] <- list("independence" = matrix(NA, simnum, p1), 
                     "exchangeable" = matrix(NA, simnum, p1), 
                     "ar1" = matrix(NA, simnum, p1))
  bias_res[[i]] <- list("independence" = matrix(NA, simnum, p1), 
                        "exchangeable" = matrix(NA, simnum, p1), 
                        "ar1" = matrix(NA, simnum, p1))
  bias_do[[i]] <- list("independence" = matrix(NA, simnum, p1), 
                       "exchangeable" = matrix(NA, simnum, p1), 
                       "ar1" = matrix(NA, simnum, p1))
  mse_res[[i]] <- list("independence" = matrix(NA, simnum, p1), 
                       "exchangeable" = matrix(NA, simnum, p1), 
                       "ar1" = matrix(NA, simnum, p1))
  mse_do[[i]] <- list("independence" = matrix(NA, simnum, p1), 
                      "exchangeable" = matrix(NA, simnum, p1), 
                      "ar1" = matrix(NA, simnum, p1))
  cov_prob[[i]] <- list("independence" = matrix(NA, simnum, p2), 
                        "exchangeable" = matrix(NA, simnum, p2), 
                        "ar1" = matrix(NA, simnum, p2))
  Cprob_do[[i]] <- list("independence" = matrix(NA, simnum, p2), 
                        "exchangeable" = matrix(NA, simnum, p2), 
                        "ar1" = matrix(NA, simnum, p2))
}
names(est_res) <- mod_types
names(se_res) <- mod_types
names(bias_res) <- mod_types
names(mse_res) <- mod_types
names(est_do) <- mod_types
names(se_do) <- mod_types
names(bias_do) <- mod_types
names(mse_do) <- mod_types

# =================== Simulation ====================== # 
set.seed(5207)
for (s in 1:simnum) {
  tryCatch({ # catch any error without breaking the simulation 
    simlist <- list()
    for(i in 1:N){
      # Baseline Covariates 
      male <- rbinom(1, size = 1, prob = 0.65) 
      age <- round(rnorm(1, mean = 50, sd = 10))
      BSA <- rnorm(1, mean = 2, sd = 0.3)
      
      # Generate exposure groups 
      logit_bav <- a_bav + b1_bav*age + b2_bav*male + b3_bav*BSA # logit of Pr(BAV = 1)
      prob_bav <- exp(logit_bav) / (1+exp(logit_bav))
      BAV <- rbinom(1, size = 1, prob = prob_bav) 
      
      age_5 <- c(1:5)+age
      
      visit <- 1:maxT
      logit_y <- a_oc + b1_oc*visit + b2_oc*BAV + b3_oc*age_5 + b4_oc*male + b5_oc*BSA + b6_oc*BAV*visit
      prob_y <- exp(logit_y) / (1+exp(logit_y)) # a vector of  marginal probabilities with dimension maxT = 5
      y <- cBern(n = 1, p = prob_y, rho = rho, type = "DCP") 
      
      records <- data.frame(id = rep(i, 5), visit = visit, y = y[1,], age = age_5, 
                            bsa = rep(BSA, 5), male = rep(male, 5), bav = rep(BAV, 5))
      
      simlist <- rbind(simlist, records)
    }
    
    uni_data <- simlist %>% group_by(id) %>% slice(1)
    
    covs <- uni_data %>% select(id, y, age, bsa, male, bav) 
    
    # Simulate time to dropout for subject i 
    drop_sim <- simsurv(lambdas = 0.499, # log(scale) = exp(-0.695)
                        gammas = 1, # shape parameter
                        x = covs, maxt = 5, 
                        betas = c(age = wb_age, bav = wb_bav, y = wb_y, 
                                  male = wb_male, bsa = wb_bsa), 
                        dist = "weibull")  %>% # Assuming exponential distribution for simplicity
      mutate(drop_vst = round(eventtime), 
             to_drop = case_when(drop_vst == 0 ~ 1, TRUE ~ drop_vst)) 
    
    drop_sim %>% group_by(to_drop) %>% summarise(n=n())
    
    data_dropped <- simlist %>% 
      left_join(drop_sim, by = "id") %>% 
      filter(visit <= to_drop) %>% 
      select(-c(eventtime, status, drop_vst, to_drop))
    
    ####################### Propensity Score Matching ##########################
    base_dat_drop <- data_dropped %>%  group_by(id) %>% slice(1)
    # Step1: Estimate propensity scores with logistic regression
    ps_model <- glm(y ~ age + male + bsa, family = binomial, data = base_dat_drop)
    
    # Step2: Perform propensity score matching 
    pps_match <- pairmatch(ps_model, data = base_dat_drop)
    
    # Step3: Get the matched data 
    matched_df <- data.frame(base_dat_drop, matched = pps_match, check.rows = TRUE) %>% 
      filter(!is.na(matched))
    
    # The number of patients in the matched dataset 
    N_patients <- nrow(matched_df) 
    
    # Step4: Get the matched data in long format 
    matched_long <- data_dropped %>% filter(id %in% matched_df$id) %>% 
      left_join(matched_df %>% select(id, matched), by = "id") 
    
    ############################# Fit GEE #################################
    for (m in 1:3) {
      str <- corstrs[m]
      # When there is no drop out
      # (1) Full covariate set 
      gee_full <- gee.output(type = "full", dat = simlist, corstr = str, N = N_patients)
      est_res[[1]][[m]][s, ] <- gee_full$Estimates[,"estimate"]
      se_res[[1]][[m]][s, ] <- gee_full$Estimates[,"se"]
      bias_res[[1]][[m]][s, ] <- gee_full$Estimates[,"bias"]
      mse_res[[1]][[m]][s, ] <- gee_full$Estimates[,"mse"]
      cov_prob[[1]][[m]][s, ] <- gee_full$Cov.Prob
      
      # (2) Reduced covariate set 
      gee_red <- gee.output(type = "reduced", dat = simlist, corstr = str, N = N_patients)
      est_res[[2]][[m]][s, ] <- gee_red$Estimates[,"estimate"]
      se_res[[2]][[m]][s, ] <- gee_red$Estimates[,"se"]
      bias_res[[2]][[m]][s, ] <- gee_red$Estimates[,"bias"]
      mse_res[[2]][[m]][s, ] <- gee_red$Estimates[,"mse"]
      cov_prob[[2]][[m]][s, ] <- gee_red$Cov.Prob
      
      # When there is drop out
      # (1) Full covariate set 
      gee_full_do <- gee.output(type = "full", dat = matched_long, corstr = str, N = N_patients)
      est_do[[1]][[m]][s, ] <- gee_full_do$Estimates[,"estimate"]
      se_do[[1]][[m]][s, ] <- gee_full_do$Estimates[,"se"]
      bias_do[[1]][[m]][s, ] <- gee_full_do$Estimates[,"bias"]
      mse_do[[1]][[m]][s, ] <- gee_full_do$Estimates[,"mse"]
      Cprob_do[[1]][[m]][s, ] <- gee_full_do$Cov.Prob
      
      # (2) Reduced covariate set 
      gee_red_do <- gee.output(type = "reduced", dat = matched_long, corstr = str, N = N_patients)
      est_do[[2]][[m]][s, ] <- gee_red_do$Estimates[,"estimate"]
      se_do[[2]][[m]][s, ] <- gee_red_do$Estimates[,"se"]
      bias_do[[2]][[m]][s, ] <- gee_red_do$Estimates[,"bias"]
      mse_do[[2]][[m]][s, ] <- gee_red_do$Estimates[,"mse"]
      Cprob_do[[2]][[m]][s, ] <- gee_red_do$Cov.Prob
    }
  }, error = function(e){
    cat("Error in iteration", s, ": ", e$message, "\n")
  })
  print(s)
}

# Table for coverage probabilities from full covariates 
full_cp <- rbind(colMeans(cov_prob[[1]]$independence, na.rm = TRUE),
                 colMeans(cov_prob[[1]]$exchangeable, na.rm = TRUE),
                 colMeans(cov_prob[[1]]$ar1, na.rm = TRUE), 
                 colMeans(Cprob_do[[1]]$independence, na.rm = TRUE),
                 colMeans(Cprob_do[[1]]$exchangeable, na.rm = TRUE),
                 colMeans(Cprob_do[[1]]$ar1, na.rm = TRUE)) 
param_full <- c("Intercept", "Age", "Male", "BSA", "BAV", "Visit", "BAV:Visit")
colnames(full_cp) <- rep(param_full, 2)
rownames(full_cp) <- rep(c("Independent", "Exchangeable", "AR1"), 2)
library(kableExtra)
kableExtra::kable(full_cp, escape = FALSE, digits = 3, 
                  caption = "True Parameter Coverage Rate from Full Model with 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  add_header_above(c(" ", "Unadjusted" = 7, "Adjusted" = 7)) %>% 
  column_spec(8, bold = T) %>% 
  column_spec(15, bold = T) %>% 
  pack_rows("With Drop out", 1, 3) %>% 
  pack_rows("Without Drop outs", 4, 6)

# Table for coverage probabilities from reduced covariates 
red_cp <- rbind(colMeans(cov_prob[[2]]$independence, na.rm = TRUE),
                colMeans(cov_prob[[2]]$exchangeable, na.rm = TRUE),
                colMeans(cov_prob[[2]]$ar1, na.rm = TRUE), 
                colMeans(Cprob_do[[2]]$independence, na.rm = TRUE),
                colMeans(Cprob_do[[2]]$exchangeable, na.rm = TRUE),
                colMeans(Cprob_do[[2]]$ar1, na.rm = TRUE)) 
red_param <- c("Intercept", "BAV", "Visit", "BAV:Visit")
colnames(red_cp) <- rep(red_param, 2)
rownames(red_cp) <- rep(c("Independent", "Exchangeable", "AR1"), 2)
kableExtra::kable(red_cp, escape = FALSE, digits = 3, 
                  caption = "True Parameter Coverage Rate from Reduced Model with 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  add_header_above(c(" ", "Unadjusted" = 4, "Adjusted" = 4)) %>% 
  column_spec(5, bold = T) %>% 
  column_spec(9, bold = T) %>% 
  pack_rows("With Drop out", 1, 3) %>% 
  pack_rows("Without Drop outs", 4, 6)


# Results for Full Set of Parameter Estimates
res_cols <- c("Mean Estimate", "SD(Estimate)", "Mean Std. Error", "SD(Std. Error)", "Mean Bias", "Mean MSE")
# 1. Full set of covariates
full_results <- cbind(est.ind = colMeans(est_res[[1]]$independence, na.rm = TRUE),
                      sd.ind = apply(est_res[[1]]$independence, 2, sd, na.rm = TRUE), 
                      se.ind = colMeans(se_res[[1]]$independence, na.rm = TRUE),
                      se.sd.ind = apply(se_res[[1]]$independence, 2, sd, na.rm = TRUE),
                      bias.ind = colMeans(bias_res[[1]]$independence, na.rm = TRUE),
                      mse.ind = colMeans(mse_res[[1]]$independence, na.rm = TRUE), 
                      est.ex = colMeans(est_res[[1]]$exchangeable, na.rm = TRUE),
                      sd.ex = apply(est_res[[1]]$exchangeable, 2, sd, na.rm = TRUE), 
                      se.ex = colMeans(se_res[[1]]$exchangeable, na.rm = TRUE),
                      se.sd.ex = apply(se_res[[1]]$exchangeable, 2, sd, na.rm = TRUE),
                      bias.ex = colMeans(bias_res[[1]]$exchangeable, na.rm = TRUE),
                      mse.ex = colMeans(mse_res[[1]]$exchangeable, na.rm = TRUE), 
                      est.ar1 = colMeans(est_res[[1]]$ar1, na.rm = TRUE),
                      sd.ar1 = apply(est_res[[1]]$ar1, 2, sd, na.rm = TRUE),
                      se.ar1 = colMeans(se_res[[1]]$AR1, na.rm = TRUE),
                      se.sd.ar1 = apply(se_res[[1]]$AR1, 2, sd, na.rm = TRUE),
                      bias.ar1 = colMeans(bias_res[[1]]$ar1, na.rm = TRUE),
                      mse.ar1 = colMeans(mse_res[[1]]$ar1, na.rm = TRUE)) %>% t()
colnames(full_results) <- param_full
rownames(full_results) <- rep(res_cols, 3)
kableExtra::kable(full_results, escape = FALSE, digits = 3, 
                  caption = "Estimation Results for Full Model without dropout from 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  column_spec(8, bold = T) %>% 
  pack_rows("Independent", 1, 6) %>% 
  pack_rows("Exchangeable", 7, 12) %>% 
  pack_rows("AR1", 13, 18)

# Results for Full Set of Parameter Estimates with Drop outs 
res_cols <- c("Mean Estimate", "SD(Estimate)", "Mean Std. Error", "SD(Std. Error)", "Mean Bias", "Mean MSE")
# 1. Full set of covariates
full_results_do <- cbind(est.ind = colMeans(est_do[[1]]$independence, na.rm = TRUE),
                      sd.ind = apply(est_do[[1]]$independence, 2, sd, na.rm = TRUE), 
                      se.ind = colMeans(se_do[[1]]$independence, na.rm = TRUE),
                      se.sd.ind = apply(se_do[[1]]$independence, 2, sd, na.rm = TRUE),
                      bias.ind = colMeans(bias_do[[1]]$independence, na.rm = TRUE),
                      mse.ind = colMeans(mse_do[[1]]$independence, na.rm = TRUE), 
                      est.ex = colMeans(est_do[[1]]$exchangeable, na.rm = TRUE),
                      sd.ex = apply(est_do[[1]]$exchangeable, 2, sd, na.rm = TRUE), 
                      se.ex = colMeans(se_do[[1]]$exchangeable, na.rm = TRUE),
                      se.sd.ex = apply(se_do[[1]]$exchangeable, 2, sd, na.rm = TRUE),
                      bias.ex = colMeans(bias_do[[1]]$exchangeable, na.rm = TRUE),
                      mse.ex = colMeans(mse_do[[1]]$exchangeable, na.rm = TRUE), 
                      est.ar1 = colMeans(est_do[[1]]$ar1, na.rm = TRUE),
                      sd.ar1 = apply(est_do[[1]]$ar1, 2, sd, na.rm = TRUE),
                      se.ar1 = colMeans(se_do[[1]]$ar1, na.rm = TRUE),
                      se.sd.ar1 = apply(se_do[[1]]$ar1, 2, sd, na.rm = TRUE),
                      bias.ar1 = colMeans(bias_do[[1]]$ar1, na.rm = TRUE),
                      mse.ar1 = colMeans(mse_do[[1]]$ar1, na.rm = TRUE)) %>% t()
colnames(full_results_do) <- param_full
rownames(full_results_do) <- rep(res_cols, 3)
kableExtra::kable(full_results_do, escape = FALSE, digits = 3, 
                  caption = "Estimation Results for Full Model with Drop outs from 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  column_spec(8, bold = T) %>% 
  pack_rows("Independent", 1, 6) %>% 
  pack_rows("Exchangeable", 7, 12) %>% 
  pack_rows("AR1", 13, 18)


# 2. Reduced Covariate set 
reduced_result <- cbind(est.ind = colMeans(est_res[[2]]$independence, na.rm = TRUE),
                        sd.ind = apply(est_res[[2]]$independence, 2, sd, na.rm = TRUE), 
                        se.ind = colMeans(se_res[[2]]$independence, na.rm = TRUE),
                        se.sd.ind = apply(se_res[[2]]$independence, 2, sd, na.rm = TRUE),
                        bias.ind = colMeans(bias_res[[2]]$independence, na.rm = TRUE),
                        mse.ind = colMeans(mse_res[[2]]$independence, na.rm = TRUE), 
                        est.ex = colMeans(est_res[[2]]$exchangeable, na.rm = TRUE),
                        sd.ex = apply(est_res[[2]]$exchangeable, 2, sd, na.rm = TRUE), 
                        se.ex = colMeans(se_res[[2]]$exchangeable, na.rm = TRUE),
                        se.sd.ex = apply(se_res[[2]]$exchangeable, 2, sd, na.rm = TRUE),
                        bias.ex = colMeans(bias_res[[2]]$exchangeable, na.rm = TRUE),
                        mse.ex = colMeans(mse_res[[2]]$exchangeable, na.rm = TRUE), 
                        est.ar1 = colMeans(est_res[[2]]$ar1, na.rm = TRUE),
                        sd.ar1 = apply(est_res[[2]]$ar1, 2, sd, na.rm = TRUE),
                        se.ar1 = colMeans(se_res[[2]]$AR1, na.rm = TRUE),
                        se.sd.ar1 = apply(se_res[[2]]$AR1, 2, sd, na.rm = TRUE),
                        bias.ar1 = colMeans(bias_res[[2]]$ar1, na.rm = TRUE),
                        mse.ar1 = colMeans(mse_res[[2]]$ar1, na.rm = TRUE)) %>% t()
colnames(reduced_result) <- red_param
rownames(reduced_result) <- rep(res_cols, 3)
kableExtra::kable(reduced_result, escape = FALSE, digits = 3, 
                  caption = "Estimate Results for Reduced Model without dropouts from 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  column_spec(5, bold = T) %>% 
  pack_rows("Independent", 1, 6) %>% 
  pack_rows("Exchangeable", 7, 12) %>% 
  pack_rows("AR1", 13, 18)

# Reduced model with drop outs 
reduced_result_do <- cbind(est.ind = colMeans(est_do[[2]]$independence, na.rm = TRUE),
                        sd.ind = apply(est_do[[2]]$independence, 2, sd, na.rm = TRUE), 
                        se.ind = colMeans(se_do[[2]]$independence, na.rm = TRUE),
                        se.sd.ind = apply(se_do[[2]]$independence, 2, sd, na.rm = TRUE),
                        bias.ind = colMeans(bias_do[[2]]$independence, na.rm = TRUE),
                        mse.ind = colMeans(mse_do[[2]]$independence, na.rm = TRUE), 
                        est.ex = colMeans(est_do[[2]]$exchangeable, na.rm = TRUE),
                        sd.ex = apply(est_do[[2]]$exchangeable, 2, sd, na.rm = TRUE), 
                        se.ex = colMeans(se_do[[2]]$exchangeable, na.rm = TRUE),
                        se.sd.ex = apply(se_do[[2]]$exchangeable, 2, sd, na.rm = TRUE),
                        bias.ex = colMeans(bias_do[[2]]$exchangeable, na.rm = TRUE),
                        mse.ex = colMeans(mse_do[[2]]$exchangeable, na.rm = TRUE), 
                        est.ar1 = colMeans(est_do[[2]]$ar1, na.rm = TRUE),
                        sd.ar1 = apply(est_do[[2]]$ar1, 2, sd, na.rm = TRUE),
                        se.ar1 = colMeans(se_do[[2]]$ar1, na.rm = TRUE),
                        se.sd.ar1 = apply(se_do[[2]]$ar1, 2, sd, na.rm = TRUE),
                        bias.ar1 = colMeans(bias_do[[2]]$ar1, na.rm = TRUE),
                        mse.ar1 = colMeans(mse_do[[2]]$ar1, na.rm = TRUE)) %>% t()
colnames(reduced_result_do) <- red_param
rownames(reduced_result_do) <- rep(res_cols, 3)
kableExtra::kable(reduced_result_do, escape = FALSE, digits = 3, 
                  caption = "Estimate Results for Reduced Model with dropout from 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  column_spec(5, bold = T) %>% 
  pack_rows("Independent", 1, 6) %>% 
  pack_rows("Exchangeable", 7, 12) %>% 
  pack_rows("AR1", 13, 18)


