# title: QLS Functions 
# author: Jinyu Luo
# version: 2024-06-30

# 1. Function to create exchangeable correlation matrix
exch_cormat <- function(rho, n) {
  cor_matrix <- matrix(rho, n, n)
  diag(cor_matrix) <- 1
  return(cor_matrix)
}

# 2. Function to create AR(1) correlation matrix
ar1_cormat <- function(rho, n) {
  rho <- as.numeric(rho)
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  return(rho^exponent)
}


# 3. Function to estimate Stage 1 tau 
tau_stg1 <- function(mdat, maxT, Z, corstr, alpha0){
  Fa <- Fb <- 0
  for (i in mdat$clusterID){
    a_1 <- a_2 <- 0
    t_i1 <- nlevels(as.factor(mdat[mdat$clusterID == i & mdat$cluster.var==1,]$visit))
    t_i2 <- nlevels(as.factor(mdat[mdat$clusterID == i & mdat$cluster.var==2,]$visit))
    if (corstr == "independence") {Rinv1 <- solve(diag(t_i1))} 
    if (corstr == "ar1") {Rinv1 <- solve(ar1_cormat(alpha0, t_i1))} 
    if (corstr == "exchangeable") {Rinv1 <- solve(exch_cormat(alpha0, t_i1))}
    if (corstr == "independence") {Rinv2 <- solve(diag(t_i2))} 
    if (corstr == "ar1") {Rinv2 <- solve(ar1_cormat(alpha0, t_i2))} 
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


# 4. Function to estimate Stage 2 tau
tau_stg2 <- function(tau0){
  tau <- as.numeric( 2 * tau0 / ( 1 + tau0 ^ 2 ) )
  return(tau)
}


# 5. Function to estimate Stage 1 alpha for AR1
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


# 6. Function to estimate Stage 2 alpha for AR1
alpha_stg2_ar1 <- function(alpha0){
  alpha <- as.numeric( 2 * alpha0 / ( 1 + alpha0^2 ) )
  return(alpha)
}

# 7. Function to estimate Stage 1 alpha for exchangeable 
estalpha1_exch <- function(mdat, Z, Qinv){
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
      
      g1 <- vector()
      for(t in 1:t_ij){
        matZ <- matrix(c(matZ_i1[t],matZ_i2[t]),nrow=2)
        g1[t] <- t(matZ) %*% Qinv %*% matZ
      } 
      G1 <- sum(g1)
      
      g2 <- vector()
      G2 <- 0
      if(t_ij > 1){
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
  alpha0 <- uniroot(alphafun, c(0,1), tol = 1e-10, extendInt = "yes")$root
  return(alpha0)
}


# 8. Function for estimating Stage 2 alpha for exchangeable 
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

# 9. Function to calculate the Sigma matrix 
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
      if (corstr == "ar1") {Ri <- ar1_cormat(alpha, ti1)}
      if (corstr == "exchangeable") {Ri <- exch_cormat(alpha, ti1)}
      Fi <- kronecker(Qi,Ri)
      Sigma_i <- Fi[1:ni, 1:ni]
      Sigma_list <- c(Sigma_list, list(Sigma_i))
    }
    else {
      if (corstr == "independence") {Ri <- diag(ti2)}
      if (corstr == "ar1") {Ri <- ar1_cormat(alpha, ti2)}
      if (corstr == "exchangeable") {Ri <- exch_cormat(alpha, ti2)}
      Fi <- kronecker(Qi,Ri)
      Sigma_i <- Fi[1:ni, 1:ni]
      Sigma_list <- c(Sigma_list, list(Sigma_i))
    }
  }
  Sigma_list
}

# 10. Function to calculate beta_hat 
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

# 11. Function for sandwich estimator 
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
    if (corstr == "ar1") {Ri <- ar1_cormat(alpha, ni)}
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


# 12. The main QLS function 
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