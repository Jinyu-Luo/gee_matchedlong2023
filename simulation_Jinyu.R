# title: Cardiovascular Data Simulation 
# author: Jinyu Luo
# version: 2024-06-26

# Required R packages ----------------------------------------------------------
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

# Global Coefficients ----------------------------------------------------------
load("data/realcoefs.RData")
true_coefs <- full_ar1$coefficients$Estimate
names(true_coefs) <- c("g0", "bav", "visit", "age", "male", "bsa", "bav_visit")
corr_alpha <- full_ar1$corr$Estimate

