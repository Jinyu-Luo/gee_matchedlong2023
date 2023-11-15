#____________________________________________________________
# function to estimate stage 1 alpha for ar1 structure
#____________________________________________________________
# `mdat`: a data frame containing the longitudinal data, with columns for the pair ID, 
#     cluster variable, time variable, and outcome variable.
# `id`: a vector of pair IDs.
# `Z`: a matrix (Y-mu)/V^1/2
# `Qinv`: the inverse of the correlation(between cluster) matrix.

estalpha1_ar1 <- function(mdat, Z, Qinv) {
  Fa <- Fb <- 0
  for (i in mdat$cluster_id) { # for each pair
    S1_j <- S2_j <- S1_ja <- S1_jb <- 0
    
    t_i1 <- nlevels(as.factor(mdat[mdat$cluster_id == i & mdat$cluster.var==1,]$visit))
    t_i2 <- nlevels(as.factor(mdat[mdat$cluster_id == i & mdat$cluster.var==2,]$visit))
    t_ij <- min(c(t_i1, t_i2))
    Z_i1 <- Z[mdat$cluster_id == i & mdat$cluster.var == 1]
    matZ_i1 <- matrix(Z_i1, nrow = t_ij)
    Z_i2 <- Z[mdat$cluster_id == i & mdat$cluster.var == 2]
    matZ_i2 <- matrix(Z_i2, nrow = t_ij)
    
    if (t_ij > 1) {
      for (k in 1:(t_ij-1)) {
        matZ1 <- matrix(c(matZ_i1[k],matZ_i2[k]),nrow=2)
        matZ2 <- matrix(c(matZ_i1[k+1],matZ_i2[k+1]),nrow=2)
        S2_j <- S2_j + t(matZ1) %*% Qinv %*% matZ2
      }
      if (t_ij ==2){
        for(k in 1:t_ij){
          matZ <- matrix(c(matZ_i1[k],matZ_i2[k]),nrow=2)
          S1_j <- S1_j + t(matZ) %*% Qinv %*% matZ
        }
      }else{
        for(k in 1:t_ij){
          matZ <- matrix(c(matZ_i1[k],matZ_i2[k]),nrow=2)
          S1_ja <- S1_ja + t(matZ) %*% Qinv %*% matZ
        }
        for(k in 2:(t_ij-1)){
          matZ <- matrix(c(matZ_i1[k],matZ_i2[k]),nrow=2)
          S1_jb <- S1_jb + t(matZ) %*% Qinv %*% matZ
        }
        S1_j <- S1_ja + S1_jb
      }
    }
    
    Fa <- Fa + S1_j
    Fb <- Fb + S2_j
  }
  ### stage 1 estimate of alpha
  alpha0 <- ( Fa - sqrt( ( Fa - 2 * Fb ) * ( Fa + 2 * Fb ) ) ) / ( 2 * Fb )
  return(alpha0)
}



#____________________________________________________________
# function to estimate stage 2 alpha for ar1 structure
#____________________________________________________________

estalpha2_ar1 <- function(alpha0){
  alpha <- as.numeric( 2 * alpha0 / ( 1 + alpha0 ^ 2 ) )
  return(alpha)
}




#____________________________________________________________
# function to estimate stage 1 alpha for exchangeable structure
#____________________________________________________________

estalpha1_exch <- function(mdat, Z, Qinv){
  match.call()
  alphafun <- function(alpha){
    GG1 <- GG2 <- 0
    for (i in unique(mdat$cluster_id)){
      GG1j <- GG2j <- 0
      
      t_i1 <- nlevels(as.factor(mdat[mdat$cluster_id == i & mdat$cluster.var==1,]$visit))
      t_i2 <- nlevels(as.factor(mdat[mdat$cluster_id == i & mdat$cluster.var==2,]$visit))
      t_ij <- max(c(t_i1, t_i2))
      Z_i1 <- Z[mdat$cluster_id == i & mdat$cluster.var == 1]
      if (t_i1 < t_ij) {Z_i1 <- c(Z_i1, rep(0, t_ij - t_i1))}
      matZ_i1 <- matrix(Z_i1, nrow = t_ij) 
      Z_i2 <- Z[mdat$cluster_id == i & mdat$cluster.var == 2]
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


#____________________________________________________________
# function to estimate stage 2 alpha for exch structure
#____________________________________________________________

estalpha2_exch <- function(alpha0, mdat){
  match.call()
  alphapart1 <- alphapart2 <- 0
  
  for (i in mdat$cluster_id){
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



#____________________________________________________________
# function to estimate stage 1 tau for exchangeable structure
#____________________________________________________________
esttau1_exch <- function(mdat, maxT, Z, time.str,alpha0){
  Fa <- Fb <- 0
  for (i in mdat$cluster_id){
    a_1 <- a_2 <- 0
    t_i1 <- nlevels(as.factor(mdat[mdat$cluster_id == i & mdat$cluster.var==1,]$visit))
    t_i2 <- nlevels(as.factor(mdat[mdat$cluster_id == i & mdat$cluster.var==2,]$visit))
    if (time.str == "ind") {Rinv1 <- solve(diag(t_i1))} 
    if (time.str == "ar1") {Rinv1 <- solve(ar1_cor(alpha0, t_i1))} 
    if (time.str == "exch") {Rinv1 <- solve(exch_cor(alpha0, t_i1))}
    if (time.str == "ind") {Rinv2 <- solve(diag(t_i2))} 
    if (time.str == "ar1") {Rinv2 <- solve(ar1_cor(alpha0, t_i2))} 
    if (time.str == "exch") {Rinv2 <- solve(exch_cor(alpha0, t_i2))}
    Rinv <- solve(exch_cor(alpha0, maxT))
    
    Z_i1 <- Z[mdat$cluster_id == i & mdat$cluster.var == 1]
    matZ_i1 <- matrix(Z_i1, nrow = t_i1)
    Z_i2 <- Z[mdat$cluster_id == i & mdat$cluster.var == 2]
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

#____________________________________________________________
# function to estimate stage 2 tau for exchangeable structure
#____________________________________________________________
esttau2_exch <- function(tau0){
  tau <- as.numeric( 2 * tau0 / ( 1 + tau0 ^ 2 ) )
  return(tau)
}



#____________________________________________________________
# function for exchangeable correlaiton structure
#____________________________________________________________
exch_cor <- function(rho, n) {
  # Create an nxn matrix filled with tau
  cor_matrix <- matrix(rho, n, n)
  # Set the diagonal to 1
  diag(cor_matrix) <- 1
  return(cor_matrix)
}

#____________________________________________________________
# function for ar1 correlaiton structure
#____________________________________________________________

ar1_cor <- function(rho,n) {
  rho <- as.numeric(rho)
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}



#____________________________________________________________
# function for generating Sigma
#____________________________________________________________

Sigma <- function(data,tau, alpha, time.str, time.var){
  Sigma_list <- list()
  for (i in unique(data$cluster_id)) {
    Qi <- exch_cor(tau, 2)
    ti1 <- nlevels(as.factor(data[data$cluster_id == i & data$cluster.var==1,]$visit))
    ti2 <- nlevels(as.factor(data[data$cluster_id == i & data$cluster.var==2,]$visit))
    ni <- ti1+ti2
    if (ti1 >= ti2)
    {
      if (time.str == "ind") {Ri <- diag(ti1)}
      if (time.str == "ar1") {Ri <- ar1_cor(alpha, ti1)}
      if (time.str == "exch") {Ri <- exch_cor(alpha, ti1)}
      Fi <- kronecker(Qi,Ri)
      Sigma_i <- Fi[1:ni, 1:ni]
      Sigma_list <- c(Sigma_list, list(Sigma_i))
    }
    else {
      if (time.str == "ind") {Ri <- diag(ti2)}
      if (time.str == "ar1") {Ri <- ar1_cor(alpha, ti2)}
      if (time.str == "exch") {Ri <- exch_cor(alpha, ti2)}
      Fi <- kronecker(Qi,Ri)
      Sigma_i <- Fi[1:ni, 1:ni]
      Sigma_list <- c(Sigma_list, list(Sigma_i))
    }
  }
  Sigma_list
}

#____________________________________________________________
# function for estimating beta
#____________________________________________________________

beta_hat <- function(formula,data, time.var, time.str, tau, alpha) {
  X <- model.matrix(object=formula, data = data) #design matrix
  y <- as.matrix(data$root) #response variable
  Sigma_list <- list()
  Xt_Sigma_inv_X <- list()
  Xt_Sigma_inv_y <- list()
  S <- Sigma(data=data, tau=tau, alpha=alpha, time.str=time.str, time.var=time.var)
  for (i in 1:length(S)) {
    ti1 <- nlevels(as.factor(data[data$cluster_id == i & data$cluster.var==1,]$visit))
    ti2 <- nlevels(as.factor(data[data$cluster_id == i & data$cluster.var==2,]$visit))
    if (ti1 >= ti2){
      Xi <- rbind(X[data$cluster_id==i & data$cluster.var==1,], 
                  X[data$cluster_id==i & data$cluster.var==2,])
      yi <- rbind(as.matrix(y[data$cluster_id==i & data$cluster.var==1,]), 
                  as.matrix(y[data$cluster_id==i & data$cluster.var==2]))
    }
    else {
      Xi <- rbind(X[data$cluster_id==i & data$cluster.var==2,], 
                  X[data$cluster_id==i & data$cluster.var==1,])
      yi <- rbind(as.matrix(y[data$cluster_id==i & data$cluster.var==2,]), 
                  as.matrix(y[data$cluster_id==i & data$cluster.var==1]))
    }
    Sigma_inv <- solve(S[[i]])
    Xt_Sigma_inv_X_i <- t(Xi) %*% Sigma_inv %*% Xi
    Xt_Sigma_inv_X[[i]] <- Xt_Sigma_inv_X_i
    Xt_Sigma_inv_y_i <- t(Xi) %*% Sigma_inv %*% yi
    Xt_Sigma_inv_y[[i]] <- Xt_Sigma_inv_y_i
  }
  return(solve(Reduce("+", Xt_Sigma_inv_X)) %*% Reduce("+",Xt_Sigma_inv_y))
}


#____________________________________________________________
# function for sandwich estimator
#____________________________________________________________
sandwich <- function(formula,data,beta_hat,alpha, time.str){
  X <- model.matrix(object=formula, data = data)
  y <- as.matrix(data$root)
  W <- list()
  mid <- list()
  for (i in unique(data$cluster_id)) {
    Xi <- X[data$cluster_id==i,]
    yi <- y[data$cluster_id==i]
    Zi <- yi - (Xi %*% as.matrix(beta_hat))
    ti1 <- nlevels(as.factor(data[data$cluster_id == i & data$cluster.var==1,]$visit))
    ti2 <- nlevels(as.factor(data[data$cluster_id == i & data$cluster.var==2,]$visit))
    ni <- ti1 + ti2
    Ai <- diag(ni) ^ (1/2)
    if (time.str == "ind") {Ri <- diag(ni)}
    if (time.str == "ar1") {Ri <- ar1_cor(alpha, ni)}
    if (time.str == "exch") {Ri <- exch_cor(alpha, ni)}
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


#____________________________________________________________
# function to fit QLS
#____________________________________________________________
qls <- function(formula, data, subject_id, cluster_id , cluster.var, time.var, maxT, time.str){
  iter <- 0
  bdiff <- c(1,1,1,1)
  alpha0  <- 0.1 # initial alpha estimate
  
  # use independent GEE to get initial beta estimates 
  init_mod <- geeglm(formula = formula, data = data, family = gaussian,
                     id = subject_id, waves = visit,corstr = "independence") 
  #summary(init_mod)
  beta0 <- as.vector(coef(init_mod))
  Z0 <- residuals(init_mod,"pearson") #init_mod$residuals Z0[1:5][1,]
  
  # compute initial tau estimate
  tau0 <- esttau1_exch(mdat=data, maxT=maxT, Z = Z0, time.str = time.str,alpha0=alpha0) 
  
  while(max(abs(bdiff)) > .00000001){
    betahat <- beta_hat(formula=formula,data=data, time.var=time.var, 
                        time.str=time.str, tau=tau0, alpha=alpha0)
    beta1 <- as.vector(betahat)
    if (all(!is.na(betahat))){bdiff <- beta1 - beta0} #***
    
    # update tau0
    Z1 <- as.matrix(data$root) - model.matrix(object=formula, data = data) %*% as.matrix(betahat)
    
    tau00 <- esttau1_exch(mdat=data, maxT=maxT, Z = Z1, time.str = time.str,alpha0=alpha0) 
    # update alpha0 (initial alpha0 for the next iteration)
    if (!is.na(tau00)) {tau0 <- tau00}
    #print(tau0)
    
    if (time.str == "ind") {alpha0 <- 0}
    if (time.str == "ar1") {alpha0 <- estalpha1_ar1(mdat=data, Z=Z1, Qinv=Qinv)}
    if (time.str == "exch") {alpha0 <- estalpha1_exch(mdat=data, Z=Z1, Qinv=Qinv)}
    
    iter <- iter + 1
    beta0 <- beta1
    # print(paste("iter:", iter, sep = " "))
    # print(paste("alpha0:",alpha0, sep = " "))
    # print(paste("tau0:",as.numeric(tau0), sep = " "))
    # print(paste("bdiff:",max(abs(bdiff)), sep = " "))
  }
  
  # after converge, get stage 2 estimates
  tau2 <- esttau2_exch(tau0)
  if (time.str == "ind") {alpha2 <- alpha0}
  if (time.str == "ar1") {alpha2 <- estalpha2_ar1(alpha0)}
  if (time.str == "exch") {alpha2 <- estalpha2_exch(alpha0, mdat = data)}
  
  betahat1 <- beta_hat(formula=formula,data=data, time.var=time.var, 
                       time.str=time.str, tau=tau2, alpha=alpha2)
  beta <- as.vector(betahat1)
  sandwich_out <- sandwich(formula = formula, data = data, beta_hat = betahat1, 
                           alpha = alpha2, time.str = time.str)
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



#____________________________________________________________
# function to calculate 95% CI for qls estiamtes
#____________________________________________________________
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


#____________________________________________________________

### generate big cohort 
set.seed(2)
simnum <- 1000
# total sample size
N <- 250
# total number of visits
maxT <- 5
# vector of visits for all patients
visit <- rep(0:(maxT-1), N)
# vector of patient ids
id <- rep(1:N, each = maxT)

# proportion of bav = 0 in hospital cohort
pbav0 = 0.1

# coef for baseline covariates for bav = 0 (Minimal overlap)
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
beta_true_s <- cbind(b1_r,b3_r,b2_r,b7_r)
p <- length(beta_true_s) # number of regression parameters


estbeta_ind <- vector("list", length = simnum)
estbeta_exch <- vector("list", length = simnum)
estbeta_ar1 <- vector("list", length = simnum)
se_ind <- vector("list", length = simnum)
se_exch <- vector("list", length = simnum)
se_ar1 <- vector("list", length = simnum)
dfse_ind <- vector("list", length = simnum)
dfse_exch <- vector("list", length = simnum)
dfse_ar1 <- vector("list", length = simnum)
relbias_ind <- vector("list", length = simnum) # Relative bias
relbias_exch <- vector("list", length = simnum)
relbias_ar1 <- vector("list", length = simnum)
mse_ind <- vector("list", length = simnum) # MSE
mse_exch <- vector("list", length = simnum)
mse_ar1 <- vector("list", length = simnum)
cp_ind <- vector("list", length = simnum) # list for for calculating coverage probability
cp_exch <- vector("list", length = simnum)
cp_ar1 <- vector("list", length = simnum)
ccp_ind <- vector("list", length = simnum)
ccp_exch <- vector("list", length = simnum)
ccp_ar1 <- vector("list", length = simnum)
outvec_ind <- vector("list", length = simnum)
outvec_exch <- vector("list", length = simnum)
outvec_ar1 <- vector("list", length = simnum)

alpha_ind <- rep(NA, simnum)
tau_ind <- rep(NA, simnum)
alpha_ar1 <- rep(NA, simnum)
tau_ar1 <- rep(NA, simnum)
alpha_exch <- rep(NA, simnum)
tau_exch <- rep(NA, simnum)

set.seed(2)
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
  
  # repeat baseline covariates maxT times
  agei <- c(age0i, age1i)
  femalei <- c(female0i, female1i)
  bsa_bli <- c(bsa_bl0i, bsa_bl1i)
  age <- rep(agei, each = maxT)
  female <- rep(femalei, each = maxT)
  bsa_bl <- rep(bsa_bli, each = maxT)
  bav <- rep(bavi, each = maxT)
  
  # generate time to death
  covs <- data.frame(agei, bavi)
  
  #---------------survtimes (20%/80%)---------
  # find hazrate that has mortality rate 0% (lambda=0), 20%(lambda=0.05) and 80%(lambda=0.35)
  survtimes <- simsurv(dist = "exponential", lambdas = 0.35, 
                       betas = c(agei = 0.01, bavi = -0.5), x = covs, maxt = 4)
  
  
  # generate outcome values
  bi <- rnorm(N, mean = 0, sd = 5) #----bi sd=0/0.5/1/2/5--------
  b <- rep(bi, each = maxT)
  e <- rnorm(N*maxT, mean = 0, sd = 1)
  # add the confounders (age, female, bsa_bl)
  root <- b1_r + b2_r*bav +b3_r * visit + b4_r * age + b5_r * female + b6_r * bsa_bl + b7_r * bav * visit + b + e
  
  simdat0 <- as.data.frame(cbind(id, visit, age, female, bsa_bl, bav, root))
  
  simdat <- left_join(simdat0, survtimes, by = "id") |>
    mutate(root = ifelse(visit > eventtime, NA, root)) 
  
  
  # Matching----
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
    left_join(matched_base %>% dplyr::select(id, matches), by = "id") %>%
    mutate(matchid = matches)
  
  # Re-arrange the dataset with cluster id and subject id
  matched_pair <- matched_long[order(matched_long$matches),] %>%
    mutate(cluster_id = rep(1:(K/2), each=(maxT*2)),
           subject_id = rep(1:K,each=maxT),
           cluster.var = rep(rep(1:2,each = maxT), (K/2)),
           order = rep(seq(1:10),(K/2))) %>%
    dplyr::select(-matches) 
  
  matched_pair <- filter(matched_pair, !is.na(root))
  
  # fit qls with indep, exch, ar1 to sample data and save estimates and ses 
  # 1) independence
  qls_ind <- qls(root ~  visit + bav + bav:visit, data = matched_pair, subject_id = subject_id, 
                 cluster_id = cluster_id, cluster.var = cluster.var, time.var=visit, maxT=maxT,time.str="ind")
  estbeta_ind[[s]] <- qls_ind$coefficients
  se_ind[[s]] <- qls_ind$se
  dfse_ind[[s]] <- sqrt(diag((K/(K-p))*qls_ind$vcov))
  ci_ind <- qls_ci(qls_ind)
  cci_ind <- adj_qls_ci(qls_ind,K)
  for (i in 1:length(beta_true_s)){
    relbias_ind[[s]][i] <- (estbeta_ind[[s]][i] - beta_true_s[i]) / beta_true_s[i]
    mse_ind[[s]][i] <- (estbeta_ind[[s]][i] - beta_true_s[i])^2
    ifelse(beta_true_s[i]>=ci_ind[i,1]&beta_true_s[i]<=ci_ind[i,2], cp_ind[[s]][i]<-1, cp_ind[[s]][i]<-0)
    ifelse(beta_true_s[i]>=cci_ind[i,1]&beta_true_s[i]<=cci_ind[i,2], ccp_ind[[s]][i]<-1, ccp_ind[[s]][i]<-0)
  }
  alpha_ind[s] <- qls_ind$alpha
  tau_ind[s] <- qls_ind$tau
  outvec_ind[[s]] <- c(s, K, unlist(estbeta_ind[[s]]), unlist(se_ind[[s]]), unlist(dfse_ind[[s]]),
                       unlist(relbias_ind[[s]]), unlist(mse_ind[[s]]), 
                       unlist(cp_ind[[s]]),unlist(ccp_ind[[s]]),alpha_ind[s],tau_ind[s])
  
  
  # 2) AR1
  qls_ar1 <- qls(root ~  visit + bav + bav:visit, data = matched_pair, subject_id = subject_id, 
                 cluster_id = cluster_id, cluster.var = cluster.var, time.var=visit, maxT=maxT,time.str="ar1")
  estbeta_ar1[[s]] <- qls_ar1$coefficients
  se_ar1[[s]] <- qls_ar1$se
  dfse_ar1[[s]] <- sqrt(diag((K/(K-p))*qls_ar1$vcov))
  ci_ar1 <- qls_ci(qls_ar1)
  cci_ar1 <- adj_qls_ci(qls_ar1,K)
  for (i in 1:length(beta_true_s)){
    relbias_ar1[[s]][i] <- (estbeta_ar1[[s]][i] - beta_true_s[i]) / beta_true_s[i]
    mse_ar1[[s]][i] <- (estbeta_ar1[[s]][i] - beta_true_s[i])^2
    ifelse(beta_true_s[i]>=ci_ar1[i,1]&beta_true_s[i]<=ci_ar1[i,2], cp_ar1[[s]][i]<-1, cp_ar1[[s]][i]<-0)
    ifelse(beta_true_s[i]>=cci_ar1[i,1]&beta_true_s[i]<=cci_ar1[i,2], ccp_ar1[[s]][i]<-1, ccp_ar1[[s]][i]<-0)
  }
  alpha_ar1[s] <- qls_ar1$alpha
  tau_ar1[s] <- qls_ar1$tau
  outvec_ar1[[s]] <- c(s, K, unlist(estbeta_ar1[[s]]), unlist(se_ar1[[s]]), unlist(dfse_ar1[[s]]),
                       unlist(relbias_ar1[[s]]), unlist(mse_ar1[[s]]), 
                       unlist(cp_ar1[[s]]),unlist(ccp_ar1[[s]]),alpha_ar1[s],tau_ar1[s])
  
  
  # 3) Exchangeable
  qls_exch <- qls(root ~  visit + bav + bav:visit, data = matched_pair, subject_id = subject_id, 
                  cluster_id = cluster_id, cluster.var = cluster.var, time.var=visit, maxT=maxT,time.str="exch")
  estbeta_exch[[s]] <- qls_exch$coefficients
  se_exch[[s]] <- qls_exch$se
  dfse_exch[[s]] <- sqrt(diag((K/(K-p))*qls_exch$vcov))
  ci_exch <- qls_ci(qls_exch)
  cci_exch <- adj_qls_ci(qls_exch,K)
  for (i in 1:length(beta_true_s)){
    relbias_exch[[s]][i] <- (estbeta_exch[[s]][i] - beta_true_s[i]) / beta_true_s[i]
    mse_exch[[s]][i] <- (estbeta_exch[[s]][i] - beta_true_s[i])^2
    ifelse(beta_true_s[i]>=ci_exch[i,1]&beta_true_s[i]<=ci_exch[i,2], cp_exch[[s]][i]<-1, cp_exch[[s]][i]<-0)
    ifelse(beta_true_s[i]>=cci_exch[i,1]&beta_true_s[i]<=cci_exch[i,2], ccp_exch[[s]][i]<-1, ccp_exch[[s]][i]<-0)
  }
  alpha_exch[s] <- qls_exch$alpha
  tau_exch[s] <- qls_exch$tau
  outvec_exch[[s]] <- c(s, K, unlist(estbeta_exch[[s]]), unlist(se_exch[[s]]), unlist(dfse_ar1[[s]]),
                        unlist(relbias_exch[[s]]), unlist(mse_exch[[s]]), 
                        unlist(cp_exch[[s]]),unlist(ccp_exch[[s]]),alpha_exch[s],tau_exch[s])
  print(s)
}

# output
names <- c("simnum", "sample_size", "(intercept)", "visit", "bav", "vist:bav",
           "SE(beta1)", "SE(beta2)", "SE(beta3)","SE(beta4)",
           "DF-CorrectedSE(1)", "DF-CorrectedSE(beta2)", "DF-CorrectedSE(beta3)","DF-CorrectedSE(beta4)",
           "RelBias(b1)","RelBias(b2)","RelBias(b3)", "RelBias(b4)",
           "MSE(b1)","MSE(b2)","MSE(b3)","MSE(b4)" ,
           "CovProb(b1)", "CovProb(b2)","CovProb(b3)","CovProb(b4)",
           "AdjCovProb(b1)","AdjCovProb(b2)","AdjCovProb(b3)","AdjCovProb(b4)","alpha","tau")

outvec_ind <- do.call("rbind", outvec_ind) 
colnames(outvec_ind) <- names
outvec_exch <- do.call("rbind",outvec_exch)
colnames(outvec_exch) <- names
outvec_ar1 <- do.call("rbind", outvec_ar1)
colnames(outvec_ar1) <- names


# save output
outvec <- rbind(outvec_ind,outvec_exch,outvec_ar1)
corstr <- rep(c("Independence","Exchangeable","AR(1)"), each = simnum)
sample_out <- as.data.frame(cbind(corstr, outvec))
#saveRDS(sample_out,"QLS_R8S5") # Rate0/2/8 Sig0/1/5--------



outmean_ind <- c(simnum, median(outvec_ind[,2]),colMeans(outvec_ind[,c(-1,-2)]),
                 sd(outvec_ind[,3]),sd(outvec_ind[,4]),sd(outvec_ind[,5]),sd(outvec_ind[,24]),sd(outvec_ind[,25]))
outmean_exch <- c(simnum,median(outvec_exch[,2]),colMeans(outvec_exch[,c(-1,-2)]), 
                  sd(outvec_exch[,3]),sd(outvec_exch[,4]),sd(outvec_exch[,5]),sd(outvec_exch[,24]),sd(outvec_exch[,25]))
outmean_ar1 <- c(simnum,median(outvec_ar1[,2]),colMeans(outvec_ar1[,c(-1,-2)]), 
                 sd(outvec_ar1[,3]),sd(outvec_ar1[,4]),sd(outvec_ar1[,5]),sd(outvec_ar1[,24]),sd(outvec_ar1[,25]))

Model <- c("Independence","","","Exchangeable","","", "AR(1)","","")
Parameters <- rep(c("(intercept)", "visit", "visit:bav"),3)
SampleSize <- rep(outmean_ind[2],9)
MeanEstimates <- c(outmean_ind[3:5], outmean_exch[3:5],outmean_ar1[3:5])
TrueValues <- rep(beta_true_s,3)
MeanSE <- c(outmean_ind[6:8], outmean_exch[6:8],outmean_ar1[6:8])
MeanSE_adj <- c(outmean_ind[9:11], outmean_exch[9:11],outmean_ar1[9:11])
MeanRelBias <- c(outmean_ind[12:14], outmean_exch[12:14],outmean_ar1[12:14])
MSE <- c(outmean_ind[15:17], outmean_exch[15:17],outmean_ar1[15:17])
CovProb <- c(outmean_ind[18:20], outmean_exch[18:20],outmean_ar1[18:20])
CovProb_adj <- c(outmean_ind[21:23], outmean_exch[21:23],outmean_ar1[21:23])
Meanalpha <- c(rep(outmean_ind[24],3),rep(outmean_exch[24],3), rep(outmean_ar1[24],3))
Meantau <- c(rep(outmean_ind[25],3),rep(outmean_exch[25],3), rep(outmean_ar1[25],3))
SD <- c(outmean_ind[26:28], outmean_exch[26:28],outmean_ar1[26:28])
SDalpha <- c(rep(outmean_ind[29],3),rep(outmean_exch[29],3), rep(outmean_ar1[29],3))
SDtau <- c(rep(outmean_ind[30],3),rep(outmean_exch[30],3), rep(outmean_ar1[30],3))

numout <- round(cbind(TrueValues, MeanEstimates, MeanSE, MeanSE_adj,SD, 
                      MeanRelBias, MSE, CovProb,CovProb_adj,Meanalpha,SDalpha, Meantau,SDtau),3)
tab <- data.frame(Model, SampleSize, Parameters, numout)

out_tab <- kableExtra::kable(tab, row.names=FALSE, escape = FALSE,
                             caption = "QLS models comparison for matched sample with 1000 simulations",
                             col.names = c("Model","Median sample size","Parameters","True Value","Mean Est",
                                           "Mean SE","DF-corrected SE","Emp SE","Mean Rel Bias","MSE",
                                           "Cov Prob","Adj.Cov Prob","Mean alpha", "Emp SE alpha","Mean tau","Emp SE tau")) %>% 
  kable_styling(full_width = F, position = "center") %>%
  row_spec(c(3,6,9), background = "lightgrey") %>%
  column_spec(1:2, background = "transparent") 



result <- tab[c(3,6,9),c(-1,-3,-6,-9,-10,-11)]
row.names(result) <- c("Independence","Exchangeable","AR(1)")
partial_out <- kableExtra::kable(result,escape = FALSE,
                              caption = "QLS models comparison for matched sample from 1000 simulations (BAV:Visit)",
                              col.names = c("Median sample size","True Value","Mean Est","Mean SE",
                                            "Emp SE","Cov Prob","Mean alpha","Emp SE alpha","Mean tau","Emp SE tau")) %>%
  kable_styling(full_width = F, position = "center")

