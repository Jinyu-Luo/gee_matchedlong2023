library(mvtnorm)
library(nlme)
library(simsurv)
library(CorBin)
library(tidyverse)
library(optmatch)
library(geepack)

# sample size for full data
N <- 500
nvisit <- 6

age <- round(rnorm(N, mean = 50, sd = 10))
male <- rbinom(N, size = 1, prob = 0.65)
bsa <- rnorm(N, mean = 2, sd = 0.3)

### simulate BAV from logitBAV = g0 + g1*age + g2*male + g3*bsa

g0 <- 2 # range from -2 to 2
g1 <- -0.1
g2 <- 1.25
g3 <- 3

expg <- exp(g0 + g1*age + g2*male + g3*bsa)
prBAV <- expg / (1 + expg) 
BAV <- rbinom(N, size = 1, prob = prBAV)

# check
#summary(glm(BAV ~ age + male + bsa, family = binomial("logit")))


### simulate Y from logitY = beta0 + beta1*visit + beta2*BAV + beta3*visit*BAV

beta0 <- -2
beta1 <- 0
beta2 <- -1.5
beta3 <- 0.5
rho <- 0.6

visit <- rep(1:nvisit, N)
BAVlong <- rep(BAV, each = nvisit)
agelong <- rep(age, each = nvisit)
malelong <- rep(male, each = nvisit)
bsalong <- rep(bsa, each = nvisit)
id <- rep(1:N, each = nvisit)
expbeta <- exp(beta0 + beta1*visit + beta2*BAVlong + beta3*visit*BAVlong)
prY <- expbeta / (1 + expbeta)
Yi <- list()
for(i in 1:N){
  prYi <- prY[which(id == i)]
  Yi[[i]] <- t(cBern(1, p = prYi, rho = 0.6, type = "DCP"))
}
Y <- do.call("rbind", Yi)
fulldata <- data.frame(id = id, visit = visit, bav = BAVlong, root = Y, age = agelong, male = malelong, bsa = bsalong)

# check
summary(geeglm(root ~ visit + bav + visit * bav, data = fulldata,
               id = id, family = binomial("logit"),
               corstr = "ar1",
               scale.fix = TRUE))


# simulate drop out at visit j based on Y_{j-1}

visdata <- fulldata |> mutate(vdrop = 6)
for(j in 1:5){
  v <- visdata |> dplyr::filter(visit == j & vdrop == 6)
  simout <- simsurv(dist = "weibull", lambdas = 0.1, gammas = 1, betas = c(root = 1.5 + 0.05 * j),
               x = v, maxt = 1) 
  visdata <- right_join(simout, visdata, by = "id") |>
    dplyr::mutate(status = replace_na(status, 0),
                  vdrop = ifelse(status == 1, j+1, vdrop)) |>
    dplyr::select(-eventtime, -status) 
}

# delete observations if j < vdrop

visdata <- visdata |>
  dplyr::filter(visit < vdrop)

# check
summary(geeglm(root ~ visit + bav + visit * bav, data = visdata,
               id = id, family = binomial("logit"),
               corstr = "ar1",
               scale.fix = TRUE))

# match patients based on baseline data

basedata <- visdata |> dplyr::filter(visit == 1)

ps_model <- glm(bav ~ age + male + bsa, family = binomial, data = basedata)
pps_match <- pairmatch(ps_model, data = basedata)
matched_df <- data.frame(basedata, matched = pps_match, check.rows = TRUE) %>% 
  filter(!is.na(matched))
N_pat <- nrow(matched_df) # number of patients in the matched set 

# final matched data

matchid <- matched_df |> dplyr::select(id, matched)
finaldata <- visdata |>
  right_join(matchid, by = "id")

summary(geeglm(root ~ visit + bav + visit * bav, data = finaldata,
               id = id, family = binomial("logit"),
               corstr = "ar1",
               scale.fix = TRUE))
  
