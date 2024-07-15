# title: Full Data Simulation 
# author: Jinyu Luo
# version: 2024-06-30
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
load("data/realcoefs.RData")

# Global Variables -------------------------------------------------------------
n_patients <- 250 
n_sim <- 1000
maxT <- 6
rho <- 0.3
alpha_ci <- 0.05

true_coefs <- red_ar1$coefficients$Estimate
names(true_coefs) <- c("g0", "bav", "visit", "bav_visit")
corr_alpha <- red_ar1$corr$Estimate
surv_coefs <- surv_coefs[,-4]
rownames(surv_coefs)[4] <- "male"

# Helper Functions ---------------------------------------------------------
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
      true_coefs[4]*bav*vst
      # true_coefs[4]*age + true_coefs[5]*male + true_coefs[6]*bsa + 
      
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

# Data Simulation --------------------------------------------------------------
set.seed(5207)

sim_df <- expand.grid(sim_id = 1:n_sim) %>%
  mutate(full_data= map(sim_id, function(id){patient_profile()})) %>%
  mutate(dropout_data = map(full_data, ~sim_dropouts(.x))) %>%
  mutate(matched_data = map(dropout_data, ~ps_match(.x)))

save(sim_df, file = "Outputs/sim_data.RData")

