# title: Simulation for Cardiovascular Data
# author: Jinyu Luo
# version: "2024-04-25 after the meeting on 24th"

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

## Global Variables ---------------------------------------------
# 1. Simulation Variables 
N.sim <- 1030  # the total number of simulation
N <- 500       # the total number of patients in each simulation 
maxT <- 5      # maximum number of visit 
pid <- rep(1:N, each = maxT) # vector of patient IDs
rho <- 0.1

# 2. Coefficients for exposure group simulation 
glm_coef <- c(b0=1.600, # intercept
              age = -0.081, male = 1.289, bsa = 2.315)

# 3. Coefficients for outcome simulation
# Full model 
gee_coef <- c(b0=-3.815, bav=-1.938, visit=0.114, age=-0.036, 
              male=1.196, bsa=1.955, bav_visit=0.345)
param_full <- c("Intercept","BAV","Visit", "Age", "Male", "BSA","BAV:Visit")

# Reduced model 
red_gee_coef <- c(gee_coef[1:3], gee_coef[7])
red_param <- c("Intercept", "BAV", "Visit",  "BAV:Visit")

# 4. Coefficients for dropouts (Survival outcome)
surv_coef1 <- c(scale=exp(0.485), age=-0.021, bav=0.482,male=-0.261,bsa=-0.275,y=0.276)
surv_coef2 <- c(scale=exp(3.321),age=-0.007,bav=-0.032,male=-0.401,bsa=-0.619,y=-0.038)
surv_coef3 <- c(scale = exp(1.342), age =-0.022,bav = 0.293,male = -0.347,bsa = -0.093,y = 0.889)

## Helper Functions ---------------------------------------------
# 1. Function to calculate CIs with 
#    unadjusted and DF-adjusted sandwich estimator
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


# 2. Function for Fitting GEE model and Extract Estimates 
gee.output <- function(type = "full", dat, corstr, N) {
  if (type == "full"){
    true_p = gee_coef
    mod = y ~ bav * visit + age + male + bsa 
    
  } else{
    true_p = red_gee_coef
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

# Function to create a list of matrices
create_matrices <- function(N.sim, p1, p2) {
  list("Estimation" = matrix(NA, N.sim, p1),
       "SE" = matrix(NA, N.sim, p1),
       "Bias" = matrix(NA, N.sim, p1), 
       "MSE" = matrix(NA, N.sim, p1), 
       "Coverage" = matrix(NA, N.sim, p2))
}

##  Result Containers -------------------------------------------
# Correlation structures
corstrs <- c("independence", "exchangeable", "ar1") 

### 1. Without Dropouts -----------------------------------------
# Full model
fullmod_res <- list(independent = create_matrices(N.sim, 7, 14),   
                    exchangeable = create_matrices(N.sim, 7, 14),
                    ar1 = create_matrices(N.sim, 7, 14))
# Reduced model
redmod_res <- list(independent = create_matrices(N.sim, 4, 8),   
                   exchangeable = create_matrices(N.sim, 4, 8),
                   ar1 = create_matrices(N.sim, 4, 8))

### 2. With Dropouts --------------------------------------------
# Full model
fullmod_resdo <- list(independent = create_matrices(N.sim, 7, 14),   
                      exchangeable = create_matrices(N.sim, 7, 14),
                      ar1 = create_matrices(N.sim, 7, 14))
# Reduced model
redmod_resdo <- list(independent = create_matrices(N.sim, 4, 8),   
                     exchangeable = create_matrices(N.sim, 4, 8),
                     ar1 = create_matrices(N.sim, 4, 8))

# Dropout container
N_dropouts <- list()

# Effective sample size 
eff_sizes <- NULL

## Simulation -------------------------------------------------
set.seed(0426)
for (s in 1:N.sim) {
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
      logit_y <- gee_coef[1]+BAV*gee_coef[2]+vst*gee_coef[3]+
        sum(base_covs[2:4]*gee_coef[4:6])+BAV*vst*gee_coef[7]
      prob_y <- exp(logit_y) / (1+exp(logit_y)) 
      y <- as.vector(cBernEx(n = 1, p = prob_y, rho = rho))
      
      records <- data.frame(id = rep(i, maxT),  age = rep(base_covs[2], maxT),
                            bav = rep(BAV, maxT),  
                            male = rep(base_covs[3], maxT),
                            bsa = rep(base_covs[4], maxT),y=y, 
                            visit = 1:maxT) 
      
      simlist <- rbind(simlist, records)
    }
    
    ### Dropouts Simulation -------------------------------------
    # 1. Simulate Patients who dropped out after the second measurement
    visit2_df <- simlist %>% filter(visit == 2)
    Z1 <- visit2_df %>% select(-visit)
    # Assuming dropout time follows exponential distribution 
    surv_prob1 <- simsurv(lambdas = surv_coef1[1], 
                          gammas = 1, # shape parameter
                          x = Z1, maxt = 1,
                          betas = surv_coef1[2:6],
                          dist = "weibull") 
    
    # record dropout numbers
    N_dropouts <- rbind(N_dropouts, table(surv_prob1$status))
    
    # remove dropouts
    to_drop <- which(surv_prob1$status == 1)
    dropped1 <- simlist %>% filter(!(visit >= 3 & id %in% to_drop))
    
    # 2. Simulate Patients who dropped out after the third measurement
    visit3_df <- dropped1 %>% filter(visit == 3 & !(id %in% to_drop))
    Z2 <- visit3_df %>% select(-visit)
    # Assuming dropout time follows exponential distribution 
    surv_prob2 <- simsurv(lambdas = surv_coef2[1], 
                          gammas = 1, # shape parameter
                          x = Z2, maxt = 1,
                          betas = surv_coef2[2:6],
                          dist = "weibull") 
    
    # record dropout numbers
    N_dropouts <- rbind(N_dropouts, table(surv_prob2$status))
    
    # remove dropouts
    to_drop <- which(surv_prob2$status == 1)
    dropped2 <- dropped1 %>% filter(!(visit >= 4 & id %in% to_drop))
    
    # 3. Simulate Patients who dropped out after the fourth measurement
    visit4_df <- dropped2 %>% filter(visit == 4 & !(id %in% to_drop))
    Z3 <- visit4_df %>% select(-visit)
    # Assuming dropout time follows exponential distribution 
    surv_prob3 <- simsurv(lambdas = surv_coef3[1], 
                          gammas = 1, # shape parameter
                          x = Z3, maxt = 1,
                          betas = surv_coef3[2:6],
                          dist = "weibull") 
    
    # record dropout numbers
    N_dropouts <- rbind(N_dropouts, table(surv_prob3$status))
    
    # remove dropouts
    to_drop <- which(surv_prob3$status == 1)
    dropped3 <- dropped2 %>% filter(!(visit > 4 & id %in% to_drop))
    
    ### Propensity Score Matching ------------------------------------
    match_function <- function(dropped = TRUE){
      if (dropped){ data = dropped3} else {data = simlist}
      base_dat <- data %>% filter(visit == 1)
      ps_model <- glm(bav ~ age + male + bsa, family = binomial, data = base_dat)
      pps_match <- pairmatch(ps_model, data = base_dat)
      matched_df <- data.frame(base_dat, matched = pps_match, check.rows = TRUE) %>% 
        filter(!is.na(matched))
      Npatients <- nrow(matched_df)
      matched_long <- simlist %>% filter(id %in% matched_df$id) %>% 
        left_join(matched_df %>% select(id, matched), by = "id") 
      return(list(N = Npatients, data = matched_long))
    }
    
    no_dropouts <- match_function(dropped = FALSE)
    dropouts <- match_function()
    
    ### GEE -----------------------------------------------------------
    for (m in 1:3) {
      str <- corstrs[m]
      # When there is no drop out
      # (1) Full covariate set 
      gee_full <- gee.output(type = "full", dat = no_dropouts$data, 
                             corstr = str, N = no_dropouts$N)
      fullmod_res[[m]]$Estimation[s,] <- gee_full$Estimates[,"estimate"]
      fullmod_res[[m]]$SE[s,] <- gee_full$Estimates[,"se"]
      fullmod_res[[m]]$Bias[s,] <- gee_full$Estimates[,"bias"]
      fullmod_res[[m]]$MSE[s,] <- gee_full$Estimates[,"mse"]
      fullmod_res[[m]]$Coverage[s,] <- gee_full$Cov.Prob
      
      # (2) Reduced covariate set 
      gee_red <- gee.output(type = "reduced", dat = no_dropouts$data, 
                            corstr = str, N = no_dropouts$N)
      redmod_res[[m]]$Estimation[s,] <- gee_red$Estimates[,"estimate"]
      redmod_res[[m]]$SE[s,] <- gee_red$Estimates[,"se"]
      redmod_res[[m]]$Bias[s,] <- gee_red$Estimates[,"bias"]
      redmod_res[[m]]$MSE[s,]<- gee_red$Estimates[,"mse"]
      redmod_res[[m]]$Coverage[s,]<- gee_red$Cov.Prob
      
      # When there is drop out
      # (1) Full covariate set 
      gee_full_do <- gee.output(type = "full", dat = dropouts$data, corstr = str, N = dropouts$N)
      fullmod_resdo[[m]]$Estimation[s,] <- gee_full$Estimates[,"estimate"]
      fullmod_resdo[[m]]$SE[s,] <- gee_full$Estimates[,"se"]
      fullmod_resdo[[m]]$Bias[s,] <- gee_full$Estimates[,"bias"]
      fullmod_resdo[[m]]$MSE[s,] <- gee_full$Estimates[,"mse"]
      fullmod_resdo[[m]]$Coverage[s,] <- gee_full$Cov.Prob
      
      # (2) Reduced covariate set 
      gee_red_do <- gee.output(type = "reduced", dat = dropouts$data, corstr = str, N = dropouts$N)
      redmod_resdo[[m]]$Estimation[s,] <- gee_red$Estimates[,"estimate"]
      redmod_resdo[[m]]$SE[s,] <- gee_red$Estimates[,"se"]
      redmod_resdo[[m]]$Bias[s,] <- gee_red$Estimates[,"bias"]
      redmod_resdo[[m]]$MSE[s,]<- gee_red$Estimates[,"mse"]
      redmod_resdo[[m]]$Coverage[s,]<- gee_red$Cov.Prob
    }
  }, error = function(e){
    cat("Error in iteration", s, ": ", e$message, "\n")
  })
  print(s)
}


## Coverage Probability -----------------------------------------
### 1. Full Model -----------------------------------------------
full_cp <- rbind(
  # Without dropouts 
  colMeans(fullmod_res$independent$Coverage, na.rm = TRUE),
  colMeans(fullmod_res$exchangeable$Coverage,  na.rm = TRUE),
  colMeans(fullmod_res$ar1$Coverage,  na.rm = TRUE), 
  # With dropouts
  colMeans(fullmod_resdo$independent$Coverage, na.rm = TRUE),
  colMeans(fullmod_resdo$exchangeable$Coverage, na.rm = TRUE),
  colMeans(fullmod_resdo$ar1$Coverage, na.rm = TRUE)) 
colnames(full_cp) <- rep(param_full, 2)
rownames(full_cp) <- rep(c("Independent", "Exchangeable", "AR1"), 2)


### 2. Reduced Model ---------------------------------------------
red_cp <- rbind(
  # Without dropouts
  colMeans(redmod_res$independent$Coverage, na.rm = TRUE),
  colMeans(redmod_res$exchangeable$Coverage, na.rm = TRUE),
  colMeans(redmod_res$ar1$Coverage,  na.rm = TRUE), 
  # With dropouts 
  colMeans(redmod_resdo$independent$Coverage, na.rm = TRUE),
  colMeans(redmod_resdo$exchangeable$Coverage, na.rm = TRUE),
  colMeans(redmod_resdo$ar1$Coverage, na.rm = TRUE)) 
colnames(red_cp) <- rep(red_param, 2)
rownames(red_cp) <- rep(c("Independent", "Exchangeable", "AR1"), 2)

save(fullmod_res, redmod_res,fullmod_resdo, redmod_resdo, N_dropouts,
     file = "~/Documents/GEE-Practicum/Outputs/raw_res.rds")
