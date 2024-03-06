set.seed(5207)

library(CorBin)
library(optmatch)
library(tidyverse)
library(geepack)
library(lme4)

simnum <- 1000 # the total number of simulation
N <- 500       # the total number of patients for each simulation 
maxT <- 5      # maximum number of visit 
visits <- rep(1:maxT, N)     # vector of visits for all patients
pid <- rep(1:N, each = maxT) # vector of patient IDs
rho <- 0.3

# Coefficients for exposure group simulation 
a_ex <- -0.4  # intercept 
b1_ex <- -0.1 # age 
b2_ex <- 1.2  # male
b3_ex <- 3    # BSA

# Coefficients for outcome simulation
a_oc <- -1     # intercept 
b1_oc <- -0.05 # visit 
b2_oc <- -2    # BAV 
b3_oc <- -0.05 # age 
b4_oc <- 1.4   # male 
b5_oc <- 0.5   # BSA
b6_oc <- 0.5   # BAV * visit

# Containers for simulation results 
simdata_list <- vector("list", simnum) 
adj_results <- matrix(NA, simnum, 6)   # estimations for adjusted GEEs
unadj_results <- matrix(NA, simnum, 6) # estimations for un-adjusted GEEs
glmm_coef <- matrix(NA, simnum, 4)  # coefficient estimates from GLMM
glmm_se <- matrix(NA, simnum, 4)   # SE estimates from GLMM


# =================== Simulation ====================== # 
for (s in 1:simnum) {
  simlist <- list()
  for(i in 1:N){
   # Baseline Covariates 
   male <- rbinom(1, size = 1, prob = 0.65) 
   age <- rnorm(1, mean = 50, sd = 10)
   BSA <- rnorm(1, mean = 2, sd = 0.3)
   
   # Generate exposure groups 
   logit_bav <- a_ex + b1_ex*age + b2_ex*male + b3_ex*BSA # logit of Pr(BAV = 1)
   prob_bav <- exp(logit_bav) / (1+exp(logit_bav))
   BAV <- rbinom(1, size = 1, prob = prob_bav) 
   
   visit <- 1:maxT
   logit_y <- a_oc + b1_oc*visit + b2_oc*BAV + b3_oc*age + b4_oc*male + b5_oc*BSA + b6_oc*BAV*visit
   prob_y <- exp(logit_y) / (1+exp(logit_y))
   y <- cBern(n = maxT, p = prob_y, rho = rho, type = "DCP") 
  
   datai <- cbind(rep(i, maxT), c(1:maxT), as.numeric(t(y)), rep(BAV, maxT), 
                  rep(male, maxT), rep(age, maxT), rep(BSA, maxT))
   simlist[[i]] <- datai
  }
  simdata <- as.data.frame(as.matrix(do.call("rbind", simlist)))
  names(simdata) <- c("id", "visit", "y", "BAV", "male", "age", "bsa")
  
  ####################### Propensity Score Matching ###########################
  match_df <- simdata %>% group_by(id) %>% slice(1) # prep data for matching
  
  # Step1: Estimate propensity scores with logistic regression
  ps_model <- glm(y ~ age + male + bsa, family = binomial, data = match_df)
  
  # Step2: Perform propensity score matching 
  pps_match <- pairmatch(ps_model, data = match_df)
  
  # Step3: Get the matched data 
  matched_df <- data.frame(match_df, matched = pps_match, check.rows = TRUE) %>% 
     filter(!is.na(matched))
  
  K <- nrow(matched_df)
  
  # Step4: Get the matched data in long format 
  matched_long <- simdata %>% filter(id %in% matched_df$id) %>% 
    left_join(matched_df %>% select(id, matched), by = "id") %>% 
    mutate(mid = matched)
  
  # # Step5: Define the interaction explicitly and add to the data 
  # #        to avoid singular design matrix 
  # matched_long$bav.visit <- interaction(matched_long$BAV,matched_long$visit)
  # 
  # # Step6: Remove the levels in interaction for which we have no observations
  # matched_long$bav.visit <- droplevels(matched_long$bav.visit)
  # 
  # # Step7: Reorder observations by id because the documentary of geeglm 
  # #        saids the function assumes that the data is sorted 
  # matched_long <- matched_long[order(matched_long$id,matched_long$visit),]
  
  ######################## Fit GEE adjusted models #############################
  ## (1) Independence Correlation Structure
  gee.ind.adj <- geeglm(y ~ age+male+bsa+BAV+visit+BAV*visit, 
                        family = binomial("logit"), id = id, 
                        data = matched_long, corstr = "independence", scale.fix = TRUE)
  estbeta.ind <- gee.ind.adj$coef[7] # get the interaction effect 
  estse.ind <- sqrt(vcov(gee.ind.adj)[7,7])
  
  ### (2) Exchangeable Correlation Structure
  gee.exch.adj <- geeglm(y ~ age+male+bsa+BAV+visit+BAV*visit,  
                         family = binomial("logit"), id = id, 
                         data = matched_long, corstr = "exchangeable", scale.fix=TRUE)
  estbeta.exch <- gee.exch.adj$coef[7]
  estse.exch <- sqrt(vcov(gee.exch.adj)[7, 7])
  
  ### (3) AR1 Correlation Structure 
  gee.ar1.adj <- geeglm(y ~ age+male+bsa+BAV+visit+BAV*visit,  
                        family = binomial("logit"), id = id, 
                        data = matched_long, corstr = "ar1", scale.fix=TRUE)
  estbeta.ar1 <- gee.ar1.adj$coef[7]
  estse.ar1 <- sqrt(vcov(gee.ar1.adj)[7,7])
  
  resultvec <- c(estbeta.ind, estse.ind, estbeta.exch, estse.exch, estbeta.ar1, estse.ar1)
  names(resultvec) <- c("estbeta.ind", "estse.ind", "estbeta.exch", "estse.exch", "estbeta.ar1", "estse.ar1")
  adj_results[s, ] <- resultvec
  
  ######################## Fit GEE Unadjusted models ###########################
  ## (1) Independence Correlation Structure
  gee.ind <- geeglm(y ~ BAV+visit+BAV*visit, family = binomial("logit"), 
                        id = id, data = matched_long, corstr = "independence", scale.fix = TRUE)
  estbeta.ind <- gee.ind$coef[3] # get the interaction effect 
  estse.ind <- sqrt(vcov(gee.ind)[3, 3])
  
  ### (2) Exchangeable Correlation Structure
  gee.exch <- geeglm(y ~ BAV+visit+BAV*visit, family = binomial("logit"), 
                         id = id, data = matched_long, corstr = "exchangeable", scale.fix=TRUE)
  estbeta.exch <- gee.exch$coef[3]
  estse.exch <- sqrt(vcov(gee.exch)[3,3])
  
  ### (3) AR1 Correlation Structure 
  gee.ar1 <- geeglm(y ~ BAV+visit+BAV*visit, family = binomial("logit"), 
                        id = id, data = matched_long, corstr = "ar1", scale.fix=TRUE)
  estbeta.ar1 <- gee.ar1$coef[3]
  estse.ar1 <- sqrt(vcov(gee.ar1)[3,3])
  
  resultvec <- c(estbeta.ind, estse.ind, estbeta.exch, estse.exch, estbeta.ar1, estse.ar1)
  names(resultvec) <- c("estbeta.ind", "estse.ind", "estbeta.exch", "estse.exch", "estbeta.ar1", "estse.ar1")
  unadj_results[s, ] <- resultvec
  
  ######################## Generalized Linear Mixed Model ######################
  glmm <- glmer(y ~ BAV+visit+BAV*visit+ (1|id), data = matched_long, nAGQ = 50,
                family = binomial, control = glmerControl(optimizer = "bobyqa"))
  glmm_coef[s, ] <- unname(coef(summary(glmm))[,1])
  glmm_se[s, ] <- unname(coef(summary(glmm))[,2])
  print(s)
}
