### generate matched sample
set.seed(2)
simnum <- 1000

# matched sample 
N <- 25
# total sample
M <- 2*N
# total number of visits
maxT <- 5
# vector of visits for all matches
visit <- rep(1:maxT, N)
# vector of matched ids
matchid <- rep(1:N, each = maxT)

# true baseline coefficients
b1_r <- 30
b2_r <- 0 #
b3_r <- -0.2
b4_r <- 0.5
beta_true <- cbind(b1_r, b3_r, b4_r)

# create lists for output
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


for (s in 1:simnum) {

  bi <- rnorm(N, mean = 0, sd = 5)
  b <- rep(bi, each = maxT)
  
  # generate outcome for exposed group (bav = 1)
  bav <- 1
  e <- rnorm(N*maxT, mean = 0, sd = 3)
  root1 <- b1_r + b2_r * bav + b3_r * visit + b4_r * bav * visit + b + e
  
  # generate outcome for exposed group (bav = 0)
  bav <- 0
  e <- rnorm(N*maxT, mean = 0, sd = 3)
  root0 <- b1_r + b2_r * bav + b3_r * visit + b4_r * bav * visit + b + e
  
  simdat <- as.data.frame(cbind(matchid, visit, root0, root1)) %>%
    pivot_longer(cols = starts_with("root"),
                 names_to = "bav",
                 names_prefix = "root",
                 values_to = "root") %>%
    group_by(matchid,visit) %>%
    mutate(patientid = row_number(),
           patientid = patientid+(matchid-1)*2) %>%
    ungroup()
  simdat$bav <- as.numeric(simdat$bav)
  simdat <- simdat %>% arrange(patientid)
  
  #1) independece
  gee_ind <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                    id = patientid, waves = visit, corstr = "independence")
  estbeta_ind[[s]] <- coef(gee_ind) # beta estimates
  se_ind[[s]] <- summary(gee_ind)$coefficient[,2] # standard errors
  ci_ind <- geeglm_ci(gee_ind) # 95% CI
  for (i in 1:length(beta_true)) {
    relbias_ind[[s]][i] <- (estbeta_ind[[s]][i] - beta_true[i]) / beta_true[i]
    mse_ind[[s]][i] <- (estbeta_ind[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_ind[i,1]&beta_true[i]<=ci_ind[i,2], cp_ind[[s]][i]<-1, cp_ind[[s]][i]<-0)
  }
  outvec_ind[[s]] <- c(s, N, M, estbeta_ind[[s]][1],estbeta_ind[[s]][2],estbeta_ind[[s]][3],
                       se_ind[[s]][1],se_ind[[s]][2],se_ind[[s]][3],relbias_ind[[s]][1],
                       relbias_ind[[s]][2],relbias_ind[[s]][3],mse_ind[[s]][1],mse_ind[[s]][2],
                       mse_ind[[s]][3],cp_ind[[s]][1],cp_ind[[s]][2],cp_ind[[s]][3])
  
  #2) exchangeable
  gee_exch <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                     id = patientid, waves = visit, corstr = "exchangeable")
  estbeta_exch[[s]] <- coef(gee_exch)
  se_exch[[s]] <- summary(gee_exch)$coefficient[,2]
  ci_exch <- geeglm_ci(gee_exch)
  for (i in 1:length(beta_true)) {
    relbias_exch[[s]][i] <- (estbeta_exch[[s]][i] - beta_true[i]) / beta_true[i]
    mse_exch[[s]][i] <- (estbeta_exch[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_exch[i,1]&beta_true[i]<=ci_exch[i,2], cp_exch[[s]][i]<-1, cp_exch[[s]][i]<-0)
  }
  outvec_exch[[s]] <- c(s, N, M, estbeta_exch[[s]][1],estbeta_exch[[s]][2],estbeta_exch[[s]][3],
                        se_exch[[s]][1],se_exch[[s]][2],se_exch[[s]][3],relbias_exch[[s]][1],
                        relbias_exch[[s]][2],relbias_exch[[s]][3],mse_exch[[s]][1],mse_exch[[s]][2],
                        mse_exch[[s]][3],cp_exch[[s]][1],cp_exch[[s]][2],cp_exch[[s]][3])
  
  #3) AR1
  gee_ar1 <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                    id = patientid, waves = visit, corstr = "ar1")
  estbeta_ar1[[s]] <- coef(gee_ar1)
  se_ar1[[s]] <- summary(gee_ar1)$coefficient[,2] 
  ci_ar1 <- geeglm_ci(gee_ar1)
  for (i in 1:length(beta_true)) {
    relbias_ar1[[s]][i] <- (estbeta_ar1[[s]][i] - beta_true[i]) / beta_true[i]
    mse_ar1[[s]][i] <- (estbeta_ar1[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_ar1[i,1]&beta_true[i]<=ci_ar1[i,2], cp_ar1[[s]][i]<-1, cp_ar1[[s]][i]<-0)
  }
  outvec_ar1[[s]] <- c(s, N, M,estbeta_ar1[[s]][1],estbeta_ar1[[s]][2],estbeta_ar1[[s]][3],
                       se_ar1[[s]][1],se_ar1[[s]][2],se_ar1[[s]][3],relbias_ar1[[s]][1],
                       relbias_ar1[[s]][2],relbias_ar1[[s]][3],mse_ar1[[s]][1],mse_ar1[[s]][2],
                       mse_ar1[[s]][3],cp_ar1[[s]][1],cp_ar1[[s]][2],cp_ar1[[s]][3])
  
  #4) unstructured
  gee_unstr <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                      id = patientid, waves = visit, corstr = "unstructured")
  estbeta_unstr[[s]] <- coef(gee_unstr)
  se_unstr[[s]] <- summary(gee_unstr)$coefficient[,2] 
  ci_unstr <- geeglm_ci(gee_unstr)
  for (i in 1:length(beta_true)) {
    relbias_unstr[[s]][i] <- (estbeta_unstr[[s]][i] - beta_true[i]) / beta_true[i]
    mse_unstr[[s]][i] <- (estbeta_unstr[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_unstr[i,1]&beta_true[i]<=ci_unstr[i,2], cp_unstr[[s]][i]<-1, cp_unstr[[s]][i]<-0)
  }
  outvec_unstr[[s]] <- c(s, N, M, estbeta_unstr[[s]][1],estbeta_unstr[[s]][2],estbeta_unstr[[s]][3],
                         se_unstr[[s]][1],se_unstr[[s]][2],se_unstr[[s]][3],relbias_unstr[[s]][1],
                         relbias_unstr[[s]][2],relbias_unstr[[s]][3],mse_unstr[[s]][1],mse_unstr[[s]][2],
                         mse_unstr[[s]][3],cp_unstr[[s]][1],cp_unstr[[s]][2],cp_unstr[[s]][3])
  
  #linear mixed effect model
  lme <- lmer(root ~  visit + bav:visit + (1|matchid), data = simdat)
  estbeta_lme[[s]] <- summary(lme)$coef[,1]
  se_lme[[s]] <- summary(lme)$coef[,2]
  ci_lme <- confint(lme)[3:5,]
  for (i in 1:length(beta_true)) {
    relbias_lme[[s]][i] <- (estbeta_lme[[s]][i] - beta_true[i]) / beta_true[i]
    mse_lme[[s]][i] <- (estbeta_lme[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_lme[i,1]&beta_true[i]<=ci_lme[i,2], cp_lme[[s]][i]<-1, cp_lme[[s]][i]<-0)
  }
  outvec_lme[[s]] <- c(s, N, M, estbeta_lme[[s]][1],estbeta_lme[[s]][2],estbeta_lme[[s]][3],
                       se_lme[[s]][1],se_lme[[s]][2],se_lme[[s]][3],relbias_lme[[s]][1],
                       relbias_lme[[s]][2],relbias_lme[[s]][3],mse_lme[[s]][1],mse_lme[[s]][2],
                       mse_lme[[s]][3],cp_lme[[s]][1],cp_lme[[s]][2],cp_lme[[s]][3])
  
  print(s)
}

names <- c("simnum","matched sample","total sample", "estbeta0", "estbeta1", "estbeta2",
           "SE(beta0)", "SE(beta1)", "SE(beta2)","RelBias(b0)","RelBias(b1)",
           "RelBias(b2)", "MSE(b0)","MSE(b1)","MSE(b2)","CovProb(b0)",
           "CovProb(b1)","CovProb(b2)")

outvec_ind <- do.call("rbind", outvec_ind)
colnames(outvec_ind) <- names
outvec_exch <- do.call("rbind",outvec_exch)
colnames(outvec_exch) <- names
outvec_ar1 <- do.call("rbind", outvec_ar1)
colnames(outvec_ar1) <- names
outvec_unstr <- do.call("rbind", outvec_unstr)
colnames(outvec_unstr) <- names
outvec_lme <- do.call("rbind", outvec_lme)
colnames(outvec_lme) <- names

outmean_ind <- c(simnum,colMeans(outvec_ind[,-1]), sd(outvec_ind[,"estbeta0"]),
                 sd(outvec_ind[,"estbeta1"]),sd(outvec_ind[,"estbeta2"]))
outmean_exch <- c(simnum,colMeans(outvec_exch[,-1]),sd(outvec_exch[,"estbeta0"]),
                  sd(outvec_exch[,"estbeta1"]),sd(outvec_exch[,"estbeta2"]))
outmean_ar1 <- c(simnum,colMeans(outvec_ar1[,-1]),sd(outvec_ar1[,"estbeta0"]),
                 sd(outvec_ar1[,"estbeta1"]),sd(outvec_ar1[,"estbeta2"]))
outmean_unstr <- c(simnum,colMeans(outvec_unstr[,-1]),sd(outvec_unstr[,"estbeta0"]),
                   sd(outvec_unstr[,"estbeta1"]),sd(outvec_unstr[,"estbeta2"]))
outmean_lme <- c(simnum, colMeans(outvec_lme[,-1]), sd(outvec_lme[,"estbeta0"]),
                 sd(outvec_lme[,"estbeta1"]),sd(outvec_lme[,"estbeta2"]))


Model <- c("Independence","","","Exchangeable","","", "AR(1)","","","Unstructured","","",
           "Linear mixed effect","","")
Parameters <- rep(c("$\\beta_0$", "$\\beta_1$", "$\\beta_2$"),5)
MeanEstimates <- c(outmean_ind[4:6], outmean_exch[4:6],outmean_ar1[4:6],outmean_unstr[4:6],
                   outmean_lme[4:6])
TrueValues <- rep(beta_true,5)
MeanSE <- c(outmean_ind[7:9], outmean_exch[7:9],outmean_ar1[7:9],outmean_unstr[7:9],
            outmean_lme[7:9])
SD <- c(outmean_ind[19:21], outmean_exch[19:21],outmean_ar1[19:21],outmean_unstr[19:21],
        outmean_lme[19:21])
MeanRelBias <- c(outmean_ind[10:12], outmean_exch[10:12],outmean_ar1[10:12],outmean_unstr[10:12],
                 outmean_lme[10:12])
MSE <- c(outmean_ind[13:15], outmean_exch[13:15],outmean_ar1[13:15],outmean_unstr[13:15],
         outmean_lme[13:15])
CovProb <- c(outmean_ind[16:18], outmean_exch[16:18],outmean_ar1[16:18],outmean_unstr[16:18],
             outmean_lme[16:18])
numout <- round(cbind(TrueValues, MeanEstimates, MeanSE, SD, MeanRelBias, MSE, CovProb),3)
MatchedSample <- rep(N, 15)
TotalSample <- rep(M,15)
t <- data.frame(Model, MatchedSample,TotalSample, Parameters, numout)

out_t1 <- kableExtra::kable(t, row.names=FALSE, escape = FALSE,
                            caption = "GEE models comparison for matched pairs with 1000 simulations",
                            col.names = c("Model","Num. pairs","Total sample","Parameters","True values","Mean estimates","Mean SE",
                                          "SD","Mean relative bias","MSE","Coverage prob")) %>% 
  kable_styling(full_width = F, position = "center") %>%
  row_spec(c(3,6,9,12,15),  background = "lightgrey") %>%
  column_spec(1:2,background = "transparent") 

out_t1



# FIT THE FOUR GEE MODELS 
# MEAN BETA
# SD of BETA
# MEAN SE(BETA)
# MEAN BIAS = mean(betahat - truebeta)
# MSE = mean((betahat - truebeta)^2)
# COVERAGE PROBABILITY 
# compute 95% confidence interval for each betahat
# see if truebeta lies in the 95% CI
# how many times does truebeta lies in the 95% CI out of the 1000 simulations

# for each simulation
# outputvec <- c(s, sample_sample, estbeta0ind, estbeta1ind, estbeta2ind, sebeta0ind, sebeta1ind, sebeta2ind, biasbeta0ind, sqerrbeta0ind, covbeta0ind,
# give names names(outputvec) <- c("simnum", "sample_size", "estbeta0id",...)
# repeat for ar1, exch, unstr)
# simresult <- outside of forloop create a list of length 1000
# simresult[[s]] <- outputvec
# output simresult
# outside of the forloop, do.call rbind
