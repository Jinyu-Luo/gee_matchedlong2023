# title: QLS Fit 
# author: Jinyu Luo
# version: 2024-06-30
rm(list = ls())
# Required R packages ----------------------------------------------------------
library(tidyverse)
library(geepack)
library(parallel)
library(survival)
library(simsurv)
library(doParallel) 
load("data/realcoefs.RData")
load("Outputs/sim_data.RData")

# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less core than available
cl <- makeCluster(num_cores)
registerDoParallel(cl)


# Global Variables -------------------------------------------------------------
maxT <- 6
rho <- 0.3
alpha_ci <- 0.05

true_coefs <- full_ar1$coefficients$Estimate
names(true_coefs) <- c("g0", "bav", "visit", "age", "male", "bsa", "bav_visit")
corr_alpha <- full_ar1$corr$Estimate
surv_coefs <- surv_coefs[,-4]
rownames(surv_coefs)[4] <- "male"


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

# Helper Functions ---------------------------------------------------------
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
                        corstr=corstr, tau=sqtau0, alpha=alpha0)
    beta1 <- as.vector(betahat)
    if (all(!is.na(betahat))){bdiff <- beta1 - beta0} #***
    
    # update tau0
    Z1 <- as.matrix(data$y) - 
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


get_qls_results <- function(df, formula, corstr, adjusted) {
  
  if(adjusted){
    print(paste("Fitting with adjusted", corstr, "correlation structure"))
  }
  
  print(paste("Fitting with unadjusted", corstr, "correlation structure"))
  
  model <- tryCatch({
    qls(formula, data = df, time.var = df$visit, maxT = 6, corstr = corstr)
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

# temp <- qls_df[[2]][[1]]
# temp_res <- get_qls_results(temp,
#                             model_specs$ar1_mdl_full$formula,
#                             model_specs$ar1_mdl_full$corstr,
#                             model_specs$ar1_mdl_full$adjusted)

# Parallel processing using foreach
qls_results <- foreach(i = seq_len(nrow(qls_df)), .packages = c('tidyverse', 'geepack', 'parallel', 'survival', 'simsurv')) %dopar% {
  df <- qls_df$matched_data[[i]]
  sim_id <- qls_df$sim_id[i]
  
  model_results <- list()
  for (mdl in names(model_specs)) {
    spec <- model_specs[[mdl]]
    model_results[[mdl]] <- get_qls_results(df, spec$formula, spec$corstr, spec$adjusted)
  }
  
  c(list(sim_id = sim_id), model_results)
}

# Stop the parallel backend
stopCluster(cl)

# Convert the results to a tibble
qls_fits_df <- tibble(
  sim_id = map(qls_results, "sim_id"),
  ind_mdl_full = map(qls_results, "ind_mdl_full"), 
  ind_mdl_red = map(qls_results, "ind_mdl_red"), 
  ar1_mdl_full = map(qls_results, "ar1_mdl_full"), 
  ar1_mdl_red = map(qls_results, "ar1_mdl_red"), 
  exch_mdl_full = map(qls_results, "exch_mdl_full"), 
  exch_mdl_red = map(qls_results, "exch_mdl_red")
)

save(qls_fits_df, file = "qls_results.RData")

# 
# # Initialize a list to store the results
# qls_results <- list()
# 
# # Loop through each dataset and fit all models
# for (i in seq_len(nrow(qls_df))) {
#   df <- qls_df$matched_data[[i]]
#   sim_id <- qls_df$sim_id[i]
#   
#   model_results <- list()
#   
#   for (mdl in names(model_specs)) {
#     spec <- model_specs[[mdl]]
#     model_results[[mdl]] <- get_qls_results(df, spec$formula, spec$corstr, spec$adjusted)
#   }
#   
#   qls_results[[i]] <- c(list(sim_id = sim_id), model_results)
#   
#   print(paste("Completed simulation", i, "out of", nrow(qls_df)))
# }

