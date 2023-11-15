# Reproduce Simulation

beta0 <- 1
beta1 <- 2
sigma <- 1
p <- 0.5
N <- 200
simnum <- 1000

estbeta <- vector("list", length = simnum)

for (s in 1:simnum){
  # generate predictor
  x <- rbinom(n = N, size = 1, prob = p)
  
  # generate error
  e <- rnorm(n = N, mean = 0, sd = sigma)
  
  # generate outcome
  y <- beta0 + beta1 * x + e
  
  simdat <- as.data.frame(cbind(y, x))
  
  outlm <- lm(y ~ x, data = simdat)
  estbeta[[s]] <- coef(outlm)
}

outbeta <- do.call("rbind", estbeta)
colMeans(outbeta)

