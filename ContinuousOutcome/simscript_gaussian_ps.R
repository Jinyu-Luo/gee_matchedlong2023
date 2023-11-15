### generate big cohort 
set.seed(2)
simnum <- 1000
# total sample size
N <- 250
# total number of visits
maxT <- 5
# vector of visits for all patients
visit <- rep(1:maxT, N)
# vector of patient ids
id <- rep(1:N, each = maxT)

# proportion of bav = 0 in hospital cohort
pbav0 = 0.1

# coef for baseline covariates (High overlap)
# b01 <- 2.2
# b02 <- -0.002
# b03 <- -0.4
# b04 <- 0.1
# # for bav = 1
# b11 <- 2.2
# b12 <- -0.004
# b13 <- -0.2
# b14 <- 0.1

# coef for baseline covariates for bav = 0 (Moderate)
# b01 <- 6 # 4-->6
# b02 <- -0.1
# b03 <- -1.3
# b04 <- 1.5
# # coef for baseline covariates for bav = 1
# b11 <- 6 # 4-->6
# b12 <- -0.1
# b13 <- -1.3
# b14 <- 1.5

# coef for baseline covariates for bav = 0 (Minimal)
b01 <- 2.7
b02 <- -0.25
b03 <- -1.2
b04 <- 6
# coef for baseline covariates for bav = 1
b11 <- 20
b12 <- -0.2
b13 <- -6
b14 <- 0.3


# coef for generating outcome values
b1_r <- 20 #intercept
b2_r <- -2 #bav
b3_r <- -0.5 #visit
b4_r <- -0.1 #age
b5_r <- -2 #sex_female
b6_r <- 2 #bsa baseline
b7_r <- 0.5 #visit:bav
beta_true <- cbind(b1_r, b2_r, b3_r, b4_r, b5_r, b6_r, b7_r)
beta_true_s <- cbind(b1_r,b3_r,b7_r)
p <- length(beta_true_s) # number of regression parameters

# create lists for output(cohort)
estbeta_ind <- vector("list", length = simnum)
estbeta_exch <- vector("list", length = simnum)
estbeta_ar1 <- vector("list", length = simnum)
estbeta_unstr <- vector("list", length = simnum)
se_ind <- vector("list", length = simnum)
se_exch <- vector("list", length = simnum)
se_ar1 <- vector("list", length = simnum)
se_unstr <- vector("list", length = simnum)

relbias_ind <- vector("list", length = simnum) # Relative bias
relbias_exch <- vector("list", length = simnum)
relbias_ar1 <- vector("list", length = simnum)
relbias_unstr <- vector("list", length = simnum)
mse_ind <- vector("list", length = simnum) # MSE
mse_exch <- vector("list", length = simnum)
mse_ar1 <- vector("list", length = simnum)
mse_unstr <- vector("list", length = simnum)
cp_ind <- vector("list", length = simnum) # list for for calculating coverage probability
cp_exch <- vector("list", length = simnum)
cp_ar1 <- vector("list", length = simnum)
cp_unstr <- vector("list", length = simnum)
outvec_ind <- vector("list", length = simnum)
outvec_exch <- vector("list", length = simnum)
outvec_ar1 <- vector("list", length = simnum)
outvec_unstr <- vector("list", length = simnum)

# create lists for output(sample)
estbeta_inds <- vector("list", length = simnum)
estbeta_exchs <- vector("list", length = simnum)
estbeta_ar1s <- vector("list", length = simnum)
estbeta_unstrs <- vector("list", length = simnum)
se_inds <- vector("list", length = simnum)
se_exchs <- vector("list", length = simnum)
se_ar1s <- vector("list", length = simnum)
se_unstrs <- vector("list", length = simnum)

cse_inds <- vector("list", length = simnum)
cse_exchs <- vector("list", length = simnum)
cse_ar1s <- vector("list", length = simnum)
cse_unstrs <- vector("list", length = simnum)

relbias_inds <- vector("list", length = simnum) # Relative bias
relbias_exchs <- vector("list", length = simnum)
relbias_ar1s <- vector("list", length = simnum)
relbias_unstrs <- vector("list", length = simnum)
mse_inds <- vector("list", length = simnum) # MSE
mse_exchs <- vector("list", length = simnum)
mse_ar1s <- vector("list", length = simnum)
mse_unstrs <- vector("list", length = simnum)
cp_inds <- vector("list", length = simnum) # list for for calculating coverage probability
cp_exchs <- vector("list", length = simnum)
cp_ar1s <- vector("list", length = simnum)
cp_unstrs <- vector("list", length = simnum)

ccp_inds <- vector("list", length = simnum) # adjusted coverage prob
ccp_exchs <- vector("list", length = simnum)
ccp_ar1s <- vector("list", length = simnum)
ccp_unstrs <- vector("list", length = simnum)

outvec_inds <- vector("list", length = simnum)
outvec_exchs <- vector("list", length = simnum)
outvec_ar1s <- vector("list", length = simnum)
outvec_unstrs <- vector("list", length = simnum)

# create lists for lme model
estbeta_lme <- vector("list", length = simnum)
se_lme <- vector("list", length = simnum)
relbias_lme <- vector("list", length = simnum)
mse_lme <- vector("list", length = simnum) 
cp_lme<- vector("list", length = simnum)
outvec_lme <- vector("list", length = simnum)

geeglm_ci <- function(model, level = 0.95) {
  # Calculate the lower and upper bounds of the confidence interval for each parameter
  lower <- coef(model) - qnorm(1-(1-0.95)/2) * summary(model)$coefficients[, "Std.err"]
  upper <- coef(model) + qnorm(1-(1-0.95)/2) * summary(model)$coefficients[, "Std.err"]
  results <- data.frame(lower = lower, upper = upper)
  return(results)
}

adj_geeglm_ci <- function(model, N, level = 0.95) {
  # Calculate the adjusted standard errors using DF-corrected sandwich estimator
  adj_se <- sqrt(diag((N/(N-p))*vcov(model)))
  # Calculate the lower and upper bounds of the confidence interval for each parameter
  lower <- coef(model) - qnorm(1-(1-0.95)/2) * adj_se
  upper <- coef(model) + qnorm(1-(1-0.95)/2) * adj_se
  results <- data.frame(lower = lower, upper = upper)
  return(results)
}

#Simulation----
for (s in 1:simnum) {
# baseline covariates
age0i <- rnorm(pbav0 * N, mean = 60, sd = 10)
female0i <- rbinom(pbav0 * N, size = 1, prob = 0.3)
bsa_bl0i <- rnorm(pbav0 * N, mean = 2, sd = 0.2)

age1i <- rnorm((1 - pbav0) * N, mean = 60, sd = 10)
female1i <- rbinom((1 - pbav0) * N, size = 1, prob = 0.3)
bsa_bl1i <- rnorm((1 - pbav0) * N, mean = 2, sd = 0.2)

ps0_xbeta <- b01 + b02 * age0i + b03 * female0i + b04 * bsa_bl0i
pscore0 <- exp(ps0_xbeta)/ (1 + exp(ps0_xbeta))
ps1_xbeta <- b11 + b12 * age1i + b13 * female1i + b14 * bsa_bl1i
pscore1 <- exp(ps1_xbeta)/ (1 + exp(ps1_xbeta))
ps_xbeta <- c(ps0_xbeta,ps1_xbeta)
pscore <- c(pscore0, pscore1)
bavi <- rbinom(N, size = 1, prob = pscore)

#check proportion of overlap
# x <- list(ps_tav = ps0_xbeta,ps_bav=ps1_xbeta)
# overlap(x,type="2",plot=TRUE)



# repeat baseline covariates maxT times
agei <- c(age0i, age1i)
femalei <- c(female0i, female1i)
bsa_bli <- c(bsa_bl0i, bsa_bl1i)
age <- rep(agei, each = maxT)
female <- rep(femalei, each = maxT)
bsa_bl <- rep(bsa_bli, each = maxT)
bav <- rep(bavi, each = maxT)

# generate time to death
table(rexp(10000, rate = 0.4) > 4)
covs <- data.frame(agei, bavi)
simsurv(dist = "exponential", lambdas = 0.4, betas = c(agei = 0.01, bavi = -0.5), x = covs, maxt = 5)

# generate outcome values
bi <- rnorm(N, mean = 0, sd = 5)
b <- rep(bi, each = maxT)
e <- rnorm(N*maxT, mean = 0, sd = 1)
# add the confounders (age, female, bsa_bl)
root <- b1_r + b2_r*bav +b3_r * visit + b4_r * age + b5_r * female + b6_r * bsa_bl + b7_r * bav * visit + b + e

simdat <- as.data.frame(cbind(id, visit, age, female, bsa_bl, bav, root))

# * Matching----
# create sample data by matching patients based on ps
simdat_base <- simdat %>% group_by(id) %>% slice(1)
ps <- glm(bav ~ age + female + bsa_bl, family = binomial, data = simdat_base)
#summary(ps)
ps.pm <- pairmatch(ps, data = simdat_base)
#summary(ps.pm) 
matched_base <- data.frame(simdat_base, matches = ps.pm, check.rows = TRUE) %>%
  filter(!is.na(matches))

K <- nrow(matched_base)

#bal.tab(ps.pm, covs = subset(simdat_base, select = c(age, female, bsa_bl)),distance = ps$fitted.values) 
matched_long <- simdat[simdat$id %in% matched_base$id,] # sample data

# join the matches column from matched_base to matched_long
matched_long <- matched_long %>% 
  left_join(matched_base %>% select(id, matches), by = "id") %>%
  mutate(matchid = matches)


# fit gee with indep, exch, ar1, unstr to entire cohort and sample data and save estimates and ses 
# repeat 1000 times and get mean coefs and ses

# ** a) entire cohort----
  #1) independece
  gee_ind <- geeglm(root ~ bav + visit + age + female + bsa_bl + bav:visit, data = simdat, family = gaussian,
                    id = id, waves = visit, corstr = "independence")
  estbeta_ind[[s]] <- coef(gee_ind) # beta estimates
  se_ind[[s]] <- summary(gee_ind)$coefficient[,2] # standard errors
  ci_ind <- geeglm_ci(gee_ind) # 95% CI
  for (i in 1:length(beta_true)) {
    relbias_ind[[s]][i] <- (estbeta_ind[[s]][i] - beta_true[i]) / beta_true[i]
    mse_ind[[s]][i] <- (estbeta_ind[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_ind[i,1]&beta_true[i]<=ci_ind[i,2], cp_ind[[s]][i]<-1, cp_ind[[s]][i]<-0)
  }
  outvec_ind[[s]] <- c(s, N, unlist(estbeta_ind[[s]]), unlist(se_ind[[s]]), 
                       unlist(relbias_ind[[s]]), unlist(mse_ind[[s]]), unlist(cp_ind[[s]]))
                    
  
  #2) exchangeable
  gee_exch <- geeglm(root ~ bav + visit + age + female + bsa_bl + bav:visit, data = simdat, family = gaussian,
                    id = id, waves = visit, corstr = "exchangeable")
  estbeta_exch[[s]] <- coef(gee_exch)
  se_exch[[s]] <- summary(gee_exch)$coefficient[,2]
  ci_exch <- geeglm_ci(gee_exch)
  for (i in 1:length(beta_true)) {
    relbias_exch[[s]][i] <- (estbeta_exch[[s]][i] - beta_true[i]) / beta_true[i]
    mse_exch[[s]][i] <- (estbeta_exch[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_exch[i,1]&beta_true[i]<=ci_exch[i,2], cp_exch[[s]][i]<-1, cp_exch[[s]][i]<-0)
  }
  outvec_exch[[s]] <- c(s, N, unlist(estbeta_exch[[s]]), unlist(se_exch[[s]]), 
                        unlist(relbias_exch[[s]]), unlist(mse_exch[[s]]), unlist(cp_exch[[s]]))
  
  #3) AR1
  gee_ar1 <- geeglm(root ~ bav + visit + age + female + bsa_bl + bav:visit, data = simdat, family = gaussian,
                     id = id, waves = visit, corstr = "ar1")
  estbeta_ar1[[s]] <- coef(gee_ar1)
  se_ar1[[s]] <- summary(gee_ar1)$coefficient[,2] 
  ci_ar1 <- geeglm_ci(gee_ar1)
  for (i in 1:length(beta_true)) {
    relbias_ar1[[s]][i] <- (estbeta_ar1[[s]][i] - beta_true[i]) / beta_true[i]
    mse_ar1[[s]][i] <- (estbeta_ar1[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_ar1[i,1]&beta_true[i]<=ci_ar1[i,2], cp_ar1[[s]][i]<-1, cp_ar1[[s]][i]<-0)
  }
  outvec_ar1[[s]] <- c(s, N, unlist(estbeta_ar1[[s]]), unlist(se_ar1[[s]]), 
                       unlist(relbias_ar1[[s]]), unlist(mse_ar1[[s]]), unlist(cp_ar1[[s]]))
  
  #4) unstructured
  gee_unstr <- geeglm(root ~ bav + visit + age + female + bsa_bl + bav:visit, data = simdat, family = gaussian,
                    id = id, waves = visit, corstr = "unstructured")
  estbeta_unstr[[s]] <- coef(gee_unstr)
  se_unstr[[s]] <- summary(gee_unstr)$coefficient[,2] 
  ci_unstr <- geeglm_ci(gee_unstr)
  for (i in 1:length(beta_true)) {
    relbias_unstr[[s]][i] <- (estbeta_unstr[[s]][i] - beta_true[i]) / beta_true[i]
    mse_unstr[[s]][i] <- (estbeta_unstr[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_unstr[i,1]&beta_true[i]<=ci_unstr[i,2], cp_unstr[[s]][i]<-1, cp_unstr[[s]][i]<-0)
    
  }
  outvec_unstr[[s]] <- c(s, N, unlist(estbeta_unstr[[s]]), unlist(se_unstr[[s]]), 
                         unlist(relbias_unstr[[s]]), unlist(mse_unstr[[s]]), unlist(cp_unstr[[s]]))
  
  # ** b) sample data----
  #1) independece
  gee_inds <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                     id = id, waves = visit, corstr = "independence")
  estbeta_inds[[s]] <- coef(gee_inds) # parameter estimates
  se_inds[[s]] <- summary(gee_inds)$coefficient[,2] # standard errors
  cse_inds[[s]] <- sqrt(diag((K/(K-p))*vcov(gee_inds)))
  ci_inds <- geeglm_ci(gee_inds) # 95% CI
  cci_inds <- adj_geeglm_ci(gee_inds,K)
  for (i in 1:length(beta_true_s)) {
    relbias_inds[[s]][i] <- (estbeta_inds[[s]][i] - beta_true_s[i]) / beta_true_s[i]
    mse_inds[[s]][i] <- (estbeta_inds[[s]][i] - beta_true_s[i])^2
    ifelse(beta_true_s[i]>=ci_inds[i,1]&beta_true_s[i]<=ci_inds[i,2], cp_inds[[s]][i]<-1, cp_inds[[s]][i]<-0)
    ifelse(beta_true_s[i]>=cci_inds[i,1]&beta_true_s[i]<=cci_inds[i,2], ccp_inds[[s]][i]<-1, ccp_inds[[s]][i]<-0)
  }
  outvec_inds[[s]] <- c(s, K, unlist(estbeta_inds[[s]]), unlist(se_inds[[s]]), unlist(cse_inds[[s]]),
                        unlist(relbias_inds[[s]]), unlist(mse_inds[[s]]), unlist(cp_inds[[s]]), unlist(ccp_inds[[s]]))

  #2) exchangeable
  gee_exchs <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                      id = id, waves = visit, corstr = "exchangeable")
  estbeta_exchs[[s]] <- coef(gee_exchs)
  se_exchs[[s]] <- summary(gee_exchs)$coefficient[,2]
  cse_exchs[[s]] <- sqrt(diag((K/(K-p))*vcov(gee_exchs)))
  ci_exchs <- geeglm_ci(gee_exchs)
  cci_exchs <- adj_geeglm_ci(gee_exchs,K)
  for (i in 1:length(beta_true_s)) {
    relbias_exchs[[s]][i] <- (estbeta_exchs[[s]][i] - beta_true_s[i]) / beta_true_s[i]
    mse_exchs[[s]][i] <- (estbeta_exchs[[s]][i] - beta_true_s[i])^2
    ifelse(beta_true_s[i]>=ci_exchs[i,1]&beta_true_s[i]<=ci_exchs[i,2], cp_exchs[[s]][i]<-1, cp_exchs[[s]][i]<-0)
    ifelse(beta_true_s[i]>=cci_exchs[i,1]&beta_true_s[i]<=cci_exchs[i,2], ccp_exchs[[s]][i]<-1, ccp_exchs[[s]][i]<-0)
  }
  outvec_exchs[[s]] <- c(s, K, unlist(estbeta_exchs[[s]]), unlist(se_exchs[[s]]), unlist(cse_exchs[[s]]),
                         unlist(relbias_exchs[[s]]), unlist(mse_exchs[[s]]), unlist(cp_exchs[[s]]), unlist(ccp_exchs[[s]]))
  #3) AR1
  gee_ar1s <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                     id = id, waves = visit, corstr = "ar1")
  estbeta_ar1s[[s]] <- coef(gee_ar1s)
  se_ar1s[[s]] <- summary(gee_ar1s)$coefficient[,2] 
  cse_ar1s[[s]] <- sqrt(diag((K/(K-p))*vcov(gee_ar1s)))
  ci_ar1s <- geeglm_ci(gee_ar1s)
  cci_ar1s <- adj_geeglm_ci(gee_ar1s,K)
  for (i in 1:length(beta_true_s)) {
    relbias_ar1s[[s]][i] <- (estbeta_ar1s[[s]][i] - beta_true_s[i]) / beta_true_s[i]
    mse_ar1s[[s]][i] <- (estbeta_ar1s[[s]][i] - beta_true_s[i])^2
    ifelse(beta_true_s[i]>=ci_ar1s[i,1]&beta_true_s[i]<=ci_ar1s[i,2], cp_ar1s[[s]][i]<-1, cp_ar1s[[s]][i]<-0)
    ifelse(beta_true_s[i]>=cci_ar1s[i,1]&beta_true_s[i]<=cci_ar1s[i,2], ccp_ar1s[[s]][i]<-1, ccp_ar1s[[s]][i]<-0)
  }
  outvec_ar1s[[s]] <- c(s, K, unlist(estbeta_ar1s[[s]]), unlist(se_ar1s[[s]]), unlist(cse_ar1s[[s]]),
                        unlist(relbias_ar1s[[s]]), unlist(mse_ar1s[[s]]), unlist(cp_ar1s[[s]]), unlist(ccp_ar1s[[s]]))
  
  #4) Unstructur
  gee_unstrs <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                       id = id, waves = visit, corstr = "unstructured")
  estbeta_unstrs[[s]] <- coef(gee_unstrs)
  se_unstrs[[s]] <- summary(gee_unstrs)$coefficient[,2] 
  cse_unstrs[[s]] <- sqrt(diag((K/(K-p))*vcov(gee_unstrs)))
  ci_unstrs <- geeglm_ci(gee_unstrs)
  cci_unstrs <- adj_geeglm_ci(gee_unstrs,K)
  for (i in 1:length(beta_true_s)) {
    relbias_unstrs[[s]][i] <- (estbeta_unstrs[[s]][i] - beta_true_s[i]) / beta_true_s[i]
    mse_unstrs[[s]][i] <- (estbeta_unstrs[[s]][i] - beta_true_s[i])^2
    ifelse(beta_true_s[i]>=ci_unstrs[i,1]&beta_true_s[i]<=ci_unstrs[i,2], cp_unstrs[[s]][i]<-1, cp_unstrs[[s]][i]<-0)
    ifelse(beta_true_s[i]>=cci_unstrs[i,1]&beta_true_s[i]<=cci_unstrs[i,2], ccp_unstrs[[s]][i]<-1, ccp_unstrs[[s]][i]<-0)
  }
  outvec_unstrs[[s]] <- c(s, K, unlist(estbeta_unstrs[[s]]), unlist(se_unstrs[[s]]), unlist(cse_unstrs[[s]]),
                          unlist(relbias_unstrs[[s]]), unlist(mse_unstrs[[s]]), unlist(cp_unstrs[[s]]), unlist(ccp_unstrs[[s]]))
  
  #linear mixed effect model
  lme <- lmer(root ~  visit + bav:visit + (1|id), data = matched_long)
  estbeta_lme[[s]] <- summary(lme)$coef[,1]
  se_lme[[s]] <- summary(lme)$coef[,2]
  ci_lme <- confint(lme)[3:5,]
  for (i in 1:length(beta_true_s)) {
    relbias_lme[[s]][i] <- (estbeta_lme[[s]][i] - beta_true_s[i]) / beta_true_s[i]
    mse_lme[[s]][i] <- (estbeta_lme[[s]][i] - beta_true_s[i])^2
    ifelse(beta_true_s[i]>=ci_lme[i,1]&beta_true_s[i]<=ci_lme[i,2], cp_lme[[s]][i]<-1, cp_lme[[s]][i]<-0)
  }
  outvec_lme[[s]] <- c(s,K, unlist(estbeta_lme[[s]]),unlist(se_lme[[s]]),NA,NA,NA,
                       unlist(relbias_lme[[s]]),unlist(mse_lme[[s]]), unlist(cp_lme[[s]]),NA,NA,NA)
    
  print(s)
}

# output ----
# 1)  entire cohort outcome
names <- c("simnum", "sample_size","(intercept)", "bav","visit", "age","female", "bsa_bl", "visit:bav",
           "SE(beta1)", "SE(beta2)", "SE(beta3)","SE(beta4)", "SE(beta5)", "SE(beta6)", "SE(beta7)",
           "RelBias(b1)","RelBias(b2)","RelBias(b3)", "RelBias(b4)","RelBias(b5)","RelBias(b6)", "RelBias(b7)",
           "MSE(b1)","MSE(b2)","MSE(b3)","MSE(b4)","MSE(b5)","MSE(b6)","MSE(b7)",
           "CovProb(b1)","CovProb(b2)","CovProb(b3)","CovProb(b4)","CovProb(b5)","CovProb(b6)","CovProb(b7)")

outvec_ind <- do.call("rbind", outvec_ind)
colnames(outvec_ind) <- names
outvec_exch <- do.call("rbind",outvec_exch)
colnames(outvec_exch) <- names
outvec_ar1 <- do.call("rbind", outvec_ar1)
colnames(outvec_ar1) <- names
outvec_unstr <- do.call("rbind", outvec_unstr)
colnames(outvec_unstr) <- names

outmean_ind <- c(simnum, colMeans(outvec_ind[,-1]), sd(outvec_ind[,3]),sd(outvec_ind[,4]),
                 sd(outvec_ind[,5]),sd(outvec_ind[,6]),sd(outvec_ind[,7]),sd(outvec_ind[,8]),sd(outvec_ind[,9]))
outmean_exch <- c(simnum,colMeans(outvec_exch[,-1]), sd(outvec_exch[,3]),sd(outvec_exch[,4]),
                  sd(outvec_exch[,5]),sd(outvec_exch[,6]),sd(outvec_exch[,7]),sd(outvec_exch[,8]),sd(outvec_exch[,9]))
outmean_ar1 <- c(simnum,colMeans(outvec_ar1[,-1]), sd(outvec_ar1[,3]),sd(outvec_ar1[,4]),
                 sd(outvec_ar1[,5]),sd(outvec_ar1[,6]),sd(outvec_ar1[,7]),sd(outvec_ar1[,8]),sd(outvec_ar1[,9]))
outmean_unstr <- c(simnum,colMeans(outvec_unstr[,-1]), sd(outvec_unstr[,3]),sd(outvec_unstr[,4]),
                   sd(outvec_unstr[,5]),sd(outvec_unstr[,6]),sd(outvec_unstr[,7]),sd(outvec_unstr[,8]),sd(outvec_unstr[,9]))

Model <- c("Independence",rep("",6),"Exchangeable",rep("",6), "AR(1)",rep("",6),"Unstructured",rep("",6))
#Parameters <- rep(c("$\\beta_1$", "$\\beta_2$", "$\\beta_3$","$\\beta_4$","$\\beta_5$","$\\beta_6$"),4)
Parameters <- rep(c("(intercept)","BAV","visit","age","female","bsa","visit:bav"),4)
MeanEstimates <- c(outmean_ind[3:9], outmean_exch[3:9],outmean_ar1[3:9],outmean_unstr[3:9])
TrueValues <- rep(beta_true,4)
MeanSE <- c(outmean_ind[10:16], outmean_exch[10:16],outmean_ar1[10:16],outmean_unstr[10:16])
SD <- c(outmean_ind[38:44], outmean_exch[38:44],outmean_ar1[38:44],outmean_unstr[38:44])
MeanRelBias <- c(outmean_ind[17:23], outmean_exch[17:23],outmean_ar1[17:23],outmean_unstr[17:23])
MSE <- c(outmean_ind[24:30], outmean_exch[24:30],outmean_ar1[24:30],outmean_unstr[24:30])
CovProb <- c(outmean_ind[31:37], outmean_exch[31:37],outmean_ar1[31:37],outmean_unstr[31:37])
numout <- round(cbind(TrueValues, MeanEstimates, MeanSE, SD, MeanRelBias, MSE, CovProb),3)
SampleSize <- rep(N, 28)
t <- data.frame(Model, SampleSize, Parameters, numout)

out_t1 <- kableExtra::kable(t, row.names=FALSE, escape = FALSE,
                  caption = "GEE models comparison for entire cohort with 1000 simulations",
                  col.names = c("Model","Sample size","Parameters","True values","Mean estimates","Mean SE",
                                "SD","Mean relative bias","MSE","Coverage prob")) %>% 
  kable_styling(full_width = F, position = "center") %>%
  #collapse_rows(columns = 1, valign = "middle") %>%
  row_spec(c(6,12,18,24),  background = "lightgrey") %>%
  column_spec(1:2,background = "transparent") 


# 2) sample data outcome
names1 <- c("simnum", "sample_size", "(intercept)", "visit", "vist:bav",
            "SE(beta1)", "SE(beta2)", "SE(beta3)","AdjSE1","AdjSE2",
            "AdjSE3","RelBias(b1)","RelBias(b2)",
            "RelBias(b3)", "MSE(b1)","MSE(b2)","MSE(b3)","CovProb(b1)",
            "CovProb(b2)","CovProb(b3)", "AdjCovProv(b1)","AdjCovProv(b2)","AdjCovProv(b3)")
outvec_inds <- do.call("rbind", outvec_inds)
colnames(outvec_inds) <- names1
outvec_exchs <- do.call("rbind",outvec_exchs)
colnames(outvec_exchs) <- names1
outvec_ar1s <- do.call("rbind", outvec_ar1s)
colnames(outvec_ar1s) <- names1
outvec_unstrs <- do.call("rbind", outvec_unstrs)
colnames(outvec_unstrs) <- names1
outvec_lme <- do.call("rbind", outvec_lme)
colnames(outvec_lme) <- names1

outmean_inds <- c(simnum,median(outvec_inds[,2]), colMeans(outvec_inds[,c(-1,-2)]), 
                  sd(outvec_inds[,3]), sd(outvec_inds[,4]),sd(outvec_inds[,5]))
outmean_exchs <- c(simnum,median(outvec_exchs[,2]),colMeans(outvec_exchs[,c(-1,-2)]),sd(outvec_exchs[,3]),
                   sd(outvec_exchs[,4]),sd(outvec_exchs[,5]))
outmean_ar1s <- c(simnum,median(outvec_ar1s[,2]),colMeans(outvec_ar1s[,c(-1,-2)]),sd(outvec_ar1s[,3]),
                  sd(outvec_ar1s[,4]),sd(outvec_ar1s[,5]))
outmean_unstrs <- c(simnum,median(outvec_unstrs[,2]),colMeans(outvec_unstrs[,c(-1,-2)]),sd(outvec_unstrs[,3]),
                    sd(outvec_unstrs[,4]),sd(outvec_unstrs[,5]))
outmean_lme <- c(simnum,median(outvec_lme[,2]), colMeans(outvec_lme[,c(-1,-2)]), sd(outvec_lme[,3]),
                 sd(outvec_lme[,4]),sd(outvec_lme[,5]))

Model <- c("Independence","","","Exchangeable","","", "AR(1)","","","Unstructured","","",
           "Linear mixed effect","","")
Parameters <- rep(c("(intercept)", "visit", "visit:bav"),5)
SampleSize1 <- rep(outmean_inds[2],15)
MeanEstimates <- c(outmean_inds[3:5], outmean_exchs[3:5],outmean_ar1s[3:5],outmean_unstrs[3:5],outmean_lme[3:5])
TrueValues <- rep(beta_true_s,5)
MeanSE <- c(outmean_inds[6:8], outmean_exchs[6:8],outmean_ar1s[6:8],outmean_unstrs[6:8], outmean_lme[6:8])
MeanSE_adj <- c(outmean_inds[9:11], outmean_exchs[9:11],outmean_ar1s[9:11],outmean_unstrs[9:11],outmean_lme[9:11])
SD <- c(outmean_inds[24:26], outmean_exchs[24:26],outmean_ar1s[24:26],outmean_unstrs[24:26],outmean_lme[24:26])
MeanRelBias <- c(outmean_inds[12:14], outmean_exchs[12:14],outmean_ar1s[12:14],outmean_unstrs[12:14],outmean_lme[12:14])
MSE <- c(outmean_inds[15:17], outmean_exchs[15:17],outmean_ar1s[15:17],outmean_unstrs[15:17],outmean_lme[15:17])
CovProb <- c(outmean_inds[18:20], outmean_exchs[18:20],outmean_ar1s[18:20],outmean_unstrs[18:20],outmean_lme[18:20])
AdjCovProb <- c(outmean_inds[21:23], outmean_exchs[21:23],outmean_ar1s[21:23],outmean_unstrs[21:23],outmean_lme[21:23])

numout <- round(cbind(TrueValues, MeanEstimates, MeanSE, MeanSE_adj, SD, MeanRelBias, MSE, CovProb,AdjCovProb),3)
t2 <- data.frame(Model, SampleSize1, Parameters, numout)

out_t2 <- kableExtra::kable(t2, row.names=FALSE, escape = FALSE,
                            caption = "GEE models comparison for matched sample with 1000 simulations",
                            col.names = c("Model","Sample size","Parameters","True values","Mean estimates",
                                          "Mean SE","Adjusted SE","SD","Mean relative bias","MSE","Coverage prob","Adj.CovProb")) %>% 
  kable_styling(full_width = F, position = "center") %>%
  #collapse_rows(columns = 1, valign = "middle") %>%
  row_spec(c(3,6,9,12,15), background = "lightgrey") %>%
  column_spec(1:2, background = "transparent") 

result <- t2[c(3,6,9,12),c(-1,-3,-6,-9,-10,-11)]
row.names(result) <- c("Independence","Exchangeable","AR(1)","Unstructured")
init_out <- kableExtra::kable(result,escape = FALSE,
                              caption = "GEE models comparison for matched sample from 1000 simulations (BAV:Visit)",
                              col.names = c("Median sample size","True Value","Mean Est","Mean SE",
                                            "Emp SE","Cov Prob")) %>%
  kable_styling(full_width = F, position = "center")

init_out


# out_t1
# out_t2

