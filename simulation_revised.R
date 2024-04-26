# title: Simulation for Cardiovascular Data
# author: Jinyu Luo
# date: "2024-04-24"

# Required R packages
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

##---------------------------------------------------------------
##                       Global Variables                      --
##---------------------------------------------------------------
N.sim <- 1050  # the total number of simulation
N <- 500       # the total number of patients in each simulation 
maxT <- 5      # maximum number of visit 
pid <- rep(1:N, each = maxT) # vector of patient IDs
rho <- 0.1

# =================== TRUE Parameter Coefficients ===================== #
# Coefficients for exposure group simulation
glm_coef <- c(b0=-0.25676, age = -0.06684, sex = 1.30776, bsa = 2.72470)

# Coefficients for outcome simulation
full_coef <- c(b0=-2.0127, bav=-2.9662, visit=0.0501, age=-0.0523, 
               sex=1.6753, bsa=1.4444, bav_visit=0.6532)

red_coef <- c(b0=-2.0127, bav=-2.9662, visit=0.0501, bav_visit=0.6532)

# Survival coefficients for mortality
death_coef <- c(scale = 4.718, age = -0.002, bav=-1.359, y=1.688)

# Survival coefficients for dropouts 
dropout_coef <- c(age = -0.019, male = -0.462, bsa = -0.456, bav=0.159, y=0.654)

# ====================== Helper Functions ======================= #
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
    true_p = full_coef
    mod = y ~ bav * visit + age + male + bsa 
    
  } else{
    true_p = red_coef
    mod = y ~ bav * visit
  }
  fit <- geeglm(mod, family = binomial("logit"), id = dat$id, data = dat, 
                corstr = corstr, scale.fix = TRUE)
  est <- geeci(fit, N)
  coveraged <- c(ifelse((est$lower <= true_p & true_p <= est$upper), 1, 0),
                 ifelse((est$adj_lower <= true_p & true_p <= est$adj_upper), 1, 0))
  est_set <- est[,1:2] %>% mutate(bias = (estimate - true_p)/true_p, mse = (estimate - true_p)^2)
  out <- list("Estimates" = est_set, "Cov.Prob" = coveraged)
  
  return(out)
}

corstrs <- c("independence", "exchangeable", "ar1")
mod_types <- c("full", "reduced")

# Function to create a list of matrices
create_matrices <- function(simnum, p) {
  list("independence" = matrix(NA, simnum, p),
       "exchangeable" = matrix(NA, simnum, p),
       "ar1" = matrix(NA, simnum, p))
}

# Initialize lists to store results
est_res <- list(); se_res <- list(); bias_res <- list(); mse_res <- list(); coverage <- list()
est_do <- list(); se_do <- list(); bias_do <- list(); mse_do <- list(); coverage_do <- list()

# Create list of matrices for above containers 
for (i in 1:2) {
  p1 <- ifelse(i == 2, 4, 7)
  p2 <- ifelse(i == 2, 8, 14)
  
  # Create matrices for each type and each scenario
  est_res[[i]] <- create_matrices(simnum, p1)
  se_res[[i]] <- create_matrices(simnum, p1)
  bias_res[[i]] <- create_matrices(simnum, p1)
  mse_res[[i]] <- create_matrices(simnum, p1)
  coverage[[i]] <- create_matrices(simnum, p2)

  est_do[[i]] <- create_matrices(simnum, p1)
  se_do[[i]] <- create_matrices(simnum, p1)
  bias_do[[i]] <- create_matrices(simnum, p1)
  mse_do[[i]] <- create_matrices(simnum, p1)
  coverage_do[[i]] <- create_matrices(simnum, p2)
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
      base_covs <- c(a = 1,
                     age = round(rnorm(1, mean = 50, sd = 10)), 
                     male = rbinom(1, size = 1, prob = 0.65) , 
                     bsa = rnorm(1, mean = 2, sd = 0.3))
      logit_bav <- sum(base_covs*glm_coef)
      prob_bav <- exp(logit_bav) / (1+exp(logit_bav))
      
      # Generate Exposure Group
      BAV <- rbinom(1, size = 1, prob = prob_bav) 
      
      # Simulate the initial outcome 
      vst <- 1:maxT
      logit_y <- full_coef[1]+BAV*full_coef[2]+vst*full_coef[3]+
        sum(base_covs[2:4]*full_coef[4:6])+BAV*vst*full_coef[7]
      prob_y <- exp(logit_y) / (1+exp(logit_y)) 
      y <- as.vector(cBernEx(n = 1, p = prob_y, rho = rho))
      
      records <- data.frame(id = rep(i, 5), y=y, bav = rep(BAV, 5), 
                            visit = 1:maxT, age = rep(base_covs[2], 5), 
                            male = rep(base_covs[3], 5), 
                            bsa = rep(base_covs[4], 5)) 
      
      simlist <- rbind(simlist, records)
    }
    
    base_data <- simlist %>% group_by(id) %>% filter(visit == 1)
    
    # Simulate death  
    # death_covs <- base_data %>% select(id, age,bav,y) 
    # death_sim <- simsurv(lambda = death_coef[1], # scale parameter
    #                      gammas = 1, # shape parameter 
    #                      x = death_covs, maxt = 5, 
    #                      betas = death_coef[2:4], 
    #                      dist = "weibull") %>% 
    #     mutate(drop_vst = round(eventtime))
    # 
    # dropped_death <- simlist %>% left_join(death_sim %>% select(id, drop_vst)) %>% 
    #   filter(drop_vst!=0) %>% 
    #   group_by(id) %>% 
    #   filter(visit < drop_vst) %>% 
    #   ungroup()
    
    # ============ Simulate time to dropout for subject i ==================== #
    # drop_sim <- simsurv(lambdas = 0.499, # log(scale) = exp(-0.695)
    #                     gammas = 1, # shape parameter
    #                     x = covs, maxt = 5, 
    #                     betas = c(age = wb_age, bav = wb_bav, y = wb_y, 
    #                               male = wb_male, bsa = wb_bsa), 
    #                     dist = "weibull")  %>% # Assuming exponential distribution for simplicity
    #   mutate(drop_vst = round(eventtime), 
    #          to_drop = case_when(drop_vst == 0 ~ 1, TRUE ~ drop_vst)) 
    # 
    # drop_sim %>% group_by(to_drop) %>% summarise(n=n())
    
    # data_dropped <- simlist %>% 
    #   left_join(drop_sim, by = "id") %>% 
    #   filter(visit <= to_drop) %>% 
    #   select(-c(eventtime, status, drop_vst, to_drop))
    
    dropouts_X <- base_data %>% select(id, age, male,bsa,bav,y)
    dropouts_sim <- simsurv(lambda = death_coef[1], # scale parameter
                            gammas = 1, # shape parameter 
                            x = dropouts_X, maxt = 5, 
                            betas = dropout_coef, 
                            dist = "weibull") %>% 
      mutate(drop_vst = round(eventtime))
    
    dropped_df <- simlist %>% left_join(dropouts_sim %>% select(id, drop_vst)) %>% 
      filter(drop_vst!=0) %>% 
      group_by(id) %>% 
      filter(visit < drop_vst) %>% 
      ungroup()
    
    ####################### Propensity Score Matching ##########################
    match_function <- function(dropped = TRUE){
      if (dropped){ data = dropped_df} else {data = simlist}
      base_dat <- data %>% group_by(id) %>% slice(1)
      ps_model <- glm(bav ~ age + male + bsa, family = binomial, data = base_dat)
      pps_match <- pairmatch(ps_model, data = base_dat)
      matched_df <- data.frame(base_dat, matched = pps_match, check.rows = TRUE) %>% 
        filter(!is.na(matched))
      Npatients <- nrow(matched_df)
      matched_long <- simlist %>% filter(id %in% matched_df$id) %>% 
        left_join(matched_df %>% select(id, matched), by = "id") 
      return(list(N = Npatients, data = matched_long))
    }
    
    full_data <- match_function(dropped = FALSE)
    dropouts <- match_function()
    # base_dat <- dropped_death %>% group_by(id) %>% slice(1)
    # # Step1: Estimate propensity scores with logistic regression
    # ps_model <- glm(y ~ age + male + bsa, family = binomial, data = base_dat)
    # 
    # # Step2: Perform propensity score matching 
    # pps_match <- pairmatch(ps_model, data = base_dat)
    # 
    # # Step3: Get the matched data 
    # matched_df <- data.frame(base_dat, matched = pps_match, check.rows = TRUE) %>% 
    #   filter(!is.na(matched))
    # 
    # # The number of patients in the matched dataset 
    # N_patients <- nrow(matched_df) 
    # 
    # # Step4: Get the matched data in long format 
    # matched_long <- simlist %>% filter(id %in% matched_df$id) %>% 
    #   left_join(matched_df %>% select(id, matched), by = "id") 
    
    ############################# Fit GEE #################################
    for (m in 1:3) {
      str <- corstrs[m]
      # When there is no drop out
      # (1) Full covariate set 
      gee_full <- gee.output(type = "full", dat = full_data$data, corstr = str, N = full_data$N)
      est_res[[1]][[m]][s, ] <- gee_full$Estimates[,"estimate"]
      se_res[[1]][[m]][s, ] <- gee_full$Estimates[,"se"]
      bias_res[[1]][[m]][s, ] <- gee_full$Estimates[,"bias"]
      mse_res[[1]][[m]][s, ] <- gee_full$Estimates[,"mse"]
      coverage[[1]][[m]][s, ] <- gee_full$Cov.Prob
      
      # (2) Reduced covariate set 
      gee_red <- gee.output(type = "reduced", dat = full_data$data, corstr = str, N = full_data$N)
      est_res[[2]][[m]][s, ] <- gee_red$Estimates[,"estimate"]
      se_res[[2]][[m]][s, ] <- gee_red$Estimates[,"se"]
      bias_res[[2]][[m]][s, ] <- gee_red$Estimates[,"bias"]
      mse_res[[2]][[m]][s, ] <- gee_red$Estimates[,"mse"]
      coverage[[2]][[m]][s, ] <- gee_red$Cov.Prob
      
      # When there is drop out
      # (1) Full covariate set 
      gee_full_do <- gee.output(type = "full", dat = dropouts$data, corstr = str, N = dropouts$N)
      est_do[[1]][[m]][s, ] <- gee_full_do$Estimates[,"estimate"]
      se_do[[1]][[m]][s, ] <- gee_full_do$Estimates[,"se"]
      bias_do[[1]][[m]][s, ] <- gee_full_do$Estimates[,"bias"]
      mse_do[[1]][[m]][s, ] <- gee_full_do$Estimates[,"mse"]
      coverage_do[[1]][[m]][s, ] <- gee_full_do$Cov.Prob
      
      # (2) Reduced covariate set 
      gee_red_do <- gee.output(type = "reduced", dat = dropouts$data, corstr = str, N = dropouts$N)
      est_do[[2]][[m]][s, ] <- gee_red_do$Estimates[,"estimate"]
      se_do[[2]][[m]][s, ] <- gee_red_do$Estimates[,"se"]
      bias_do[[2]][[m]][s, ] <- gee_red_do$Estimates[,"bias"]
      mse_do[[2]][[m]][s, ] <- gee_red_do$Estimates[,"mse"]
      coverage_do[[2]][[m]][s, ] <- gee_red_do$Cov.Prob
    }
  }, error = function(e){
    cat("Error in iteration", s, ": ", e$message, "\n")
  })
  print(s)
}

# Indicator of convergence 
# Full Model 
threshold_value <- 10
ind_conv_full <- which(sapply(se_res[[1]][[1]], function(x) x > threshold_value))
ind_conv_red <- which(sapply(se_res[[2]][[1]], function(x) x > threshold_value))
ex_conv_full <- which(sapply(se_res[[1]][[2]], function(x) x > threshold_value))
ex_conv_red <- which(sapply(se_res[[2]][[2]], function(x) x > threshold_value))
ar1_conv_full <- which(sapply(se_res[[1]][[3]], function(x) x > threshold_value))
ar1_conv_red <- which(sapply(se_res[[2]][[3]], function(x) x > threshold_value))

non_converged <- list(full = list(ind_conv_full, ex_conv_full, ar1_conv_full), 
                      red = list(ind_conv_red, ex_conv_red, ar1_conv_red))

est_res$full$exchangeable <- est_res$full$exchangeable[-ex_conv_full, ]
est_res$reduced$exchangeable <- est_res$reduced$exchangeable[-ex_conv_red, ]

# Table for coverage probabilities from full covariates 
full_cp <- rbind(colMeans(coverage[[1]]$independence, na.rm = TRUE),
                 colMeans(coverage[[1]]$exchangeable, na.rm = TRUE),
                 colMeans(coverage[[1]]$ar1, na.rm = TRUE), 
                 colMeans(coverage_do[[1]]$independence, na.rm = TRUE),
                 colMeans(coverage_do[[1]]$exchangeable, na.rm = TRUE),
                 colMeans(coverage_do[[1]]$ar1, na.rm = TRUE)) 
param_full <- c("Intercept","BAV","Visit", "Age", "Male", "BSA","BAV:Visit")
colnames(full_cp) <- rep(param_full, 2)
rownames(full_cp) <- rep(c("Independent", "Exchangeable", "AR1"), 2)
library(kableExtra)
kableExtra::kable(full_cp, escape = FALSE, digits = 3, 
                  caption = "True Parameter Coverage Rate from Full Model with 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  add_header_above(c(" ", "Unadjusted Standard Error" = 7, "Adjusted Standard Error" = 7)) %>% 
  column_spec(8, bold = T) %>% 
  column_spec(15, bold = T) %>% 
  pack_rows("Without Drop out", 1, 3) %>% 
  pack_rows("With Drop outs", 4, 6)

# Table for coverage probabilities from reduced covariates 
red_cp <- rbind(colMeans(coverage[[2]]$independence, na.rm = TRUE),
                colMeans(coverage[[2]]$exchangeable, na.rm = TRUE),
                colMeans(coverage[[2]]$ar1, na.rm = TRUE), 
                colMeans(coverage_do[[2]]$independence, na.rm = TRUE),
                colMeans(coverage_do[[2]]$exchangeable, na.rm = TRUE),
                colMeans(coverage_do[[2]]$ar1, na.rm = TRUE)) 
red_param <- c("Intercept", "BAV", "Visit",  "BAV:Visit")
colnames(red_cp) <- rep(red_param, 2)
rownames(red_cp) <- rep(c("Independent", "Exchangeable", "AR1"), 2)
kableExtra::kable(red_cp, escape = FALSE, digits = 3, 
                  caption = "True Parameter Coverage Rate from Reduced Model with 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  add_header_above(c(" ", "Unadjusted Standard Error" = 4, "Adjusted Standard Error" = 4)) %>% 
  column_spec(5, bold = T) %>% 
  column_spec(9, bold = T) %>% 
  pack_rows("Without Drop out", 1, 3) %>% 
  pack_rows("With Drop outs", 4, 6)


# Results for Full Set of Parameter Estimates
res_cols <- c("Mean Estimate", "SD(Estimate)", "Mean Std. Error", "SD(Std. Error)", "Mean Bias", "Mean MSE")
# 1. Full set of covariates
full_results <- cbind(`True Values` = full_coef, 
                      est.ind = colMeans(est_res[[1]]$independence, na.rm = TRUE),
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
                      se.ar1 = colMeans(se_res[[1]]$ar1, na.rm = TRUE),
                      se.sd.ar1 = apply(se_res[[1]]$ar1, 2, sd, na.rm = TRUE),
                      bias.ar1 = colMeans(bias_res[[1]]$ar1, na.rm = TRUE),
                      mse.ar1 = colMeans(mse_res[[1]]$ar1, na.rm = TRUE)) %>% t()
colnames(full_results) <- param_full
rownames(full_results) <- c("True Value", rep(res_cols, 3))
kableExtra::kable(full_results, escape = FALSE, digits = 3, 
                  caption = "Estimation Results for Full Model without dropout from 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  column_spec(8, bold = T) %>% 
  pack_rows("Independent", 2, 7) %>% 
  pack_rows("Exchangeable", 8, 13) %>% 
  pack_rows("AR1", 14, 19)



# Results for Full Set of Parameter Estimates with Drop outs 
res_cols <- c("Mean Estimate", "SD(Estimate)", "Mean Std. Error", "SD(Std. Error)", "Mean Bias", "Mean MSE")
# 1. Full set of covariates
full_results_do <- cbind(`True values` = full_coef, 
                         est.ind = colMeans(est_do[[1]]$independence, na.rm = TRUE),
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
rownames(full_results_do) <- c("True Value", rep(res_cols, 3))
kableExtra::kable(full_results_do, escape = FALSE, digits = 3, 
                  caption = "Estimation Results for Full Model with Drop outs from 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  column_spec(8, bold = T) %>% 
  pack_rows("Independent", 2, 7) %>% 
  pack_rows("Exchangeable", 8, 13) %>% 
  pack_rows("AR1", 14, 19)


# 2. Reduced Covariate set 
reduced_result <- cbind(`True values` = red_coef, 
                        est.ind = colMeans(est_res[[2]]$independence, na.rm = TRUE),
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
                        se.ar1 = colMeans(se_res[[2]]$ar1, na.rm = TRUE),
                        se.sd.ar1 = apply(se_res[[2]]$ar1, 2, sd, na.rm = TRUE),
                        bias.ar1 = colMeans(bias_res[[2]]$ar1, na.rm = TRUE),
                        mse.ar1 = colMeans(mse_res[[2]]$ar1, na.rm = TRUE)) %>% t()
colnames(reduced_result) <- red_param
rownames(reduced_result) <- c("True Value", rep(res_cols, 3))
kableExtra::kable(reduced_result, escape = FALSE, digits = 3, 
                  caption = "Estimate Results for Reduced Model without dropouts from 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  column_spec(5, bold = T) %>% 
  pack_rows("Independent", 2, 7) %>% 
  pack_rows("Exchangeable", 8, 13) %>% 
  pack_rows("AR1", 14, 19)

# Reduced model with drop outs 
reduced_result_do <- cbind(`True values` = red_coef, 
                           est.ind = colMeans(est_do[[2]]$independence, na.rm = TRUE),
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
rownames(reduced_result_do) <- c("True Value", rep(res_cols, 3))
kableExtra::kable(reduced_result_do, escape = FALSE, digits = 3, 
                  caption = "Estimate Results for Reduced Model with dropout from 1000 Simulations") %>%
  kable_styling(full_width = F, position = "center") %>% 
  column_spec(5, bold = T) %>% 
  pack_rows("Independent", 2, 7) %>% 
  pack_rows("Exchangeable", 8, 13) %>% 
  pack_rows("AR1", 14, 19)

