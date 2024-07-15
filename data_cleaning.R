
# title: Data Cleaning
# author: Jinyu Luo
# version: "Final Version after meeting on 2024.06.26"

# Load Data -------------------------------------------------------------------
load(paste(path,"1 - data setup.RData", sep="/"))

# Required R packages ---------------------------------------------------------
library(tidyverse)
library(dplyr)
library(optmatch)
library(simsurv)
library(survival)
library(eha)
library(geepack)

# Select features from the large long table 
raw_long <- raw_long_d %>% 
  mutate(bav = ifelse(bav_confirmed == "1_BAV", "BAV", "TAV")) %>% 
  select(ptid, age, sex, bav, died, lka_d, yr2death, 
         bsa_baseline, bsa_echo, dateor, mdate, root, yr2vst, date_base) 

eda_dat <- raw_long %>% 
  group_by(ptid) %>% 
  mutate(count = row_number(), total = max(count)) %>% 
  slice(1) %>% 
  ungroup()

eda_dat %>% 
  group_by(total) %>% 
  summarise(n=n())
#    total     n
#    <int> <int>
# 1      1    54
# 2      2    76
# 3      3    42
# 4      4    31
# 5      5    20
# 6      6    15
# 7      7    13
# 8      8    12
# 9      9    10
# 10    10     5
# 11    11     6
# 12    12     2
# 13    13     5
# 14    14     1
# 15    15     1
# 16    16     3
# 17    18     2

eda_dat %>% 
  filter(total > 1) %>% 
  group_by(bav) %>% 
  summarise(n=n())
  

# Combine long table with baseline 
ar_data <- raw_long %>% 
  left_join(select(work_pt_d, ptid, bsa), by = "ptid") 

n_distinct(ar_data$ptid)

pid <- ar_data %>% select(ptid) %>% unique() %>% 
  mutate(pid = 1:298)

# Full Data 
full_data <- ar_data %>% 
  left_join(pid, by="ptid") %>% 
  relocate(pid, .before = ptid) %>% 
  select(-ptid) %>% 
  mutate(bav = ifelse(bav == "BAV", 1, 0), 
         sex = ifelse(sex == "1_male", 1, 0), 
         across(c(dateor, mdate, date_base), ~as.Date(., format="%m/%d/%Y")),
         lka_d = as.Date(lka_d, format="%d%b%Y")) %>% 
  filter(!is.na(root)) %>% 
  filter(!(lka_d < mdate)) %>%  # ensure that died patient don't have follow-up
  group_by(pid) %>% 
  mutate(visit = row_number(), total_visit = n(), 
         growth = root - lag(root, default = first(root)),
         y=case_when(root > 45 | growth > 5 ~ 1, TRUE ~ 0)) %>% 
  filter(min(yr2vst) == 0) %>% # The baseline measurement was made at the day of operation
  mutate(followup = visit - 1) %>% 
  ungroup() %>% 
  filter(total_visit != 1)

# EDA --------------------------------------------------------------------------
n_distinct(full_data$pid)
# [1] 101

# Modeling Dropouts ------------------------------------------------------------
full_data %>% 
  group_by(pid) %>% 
  mutate(event = ifelse(total_visit < 6, 1, 0)) %>% 
  slice(1) %>%
  ungroup() %>% 
  group_by(total_visit)%>% summarise(n =n()) %>% mutate(N = sum(n))
# A tibble: 12 × 3
#  total_visit     n     N
#         <int> <int> <int>
#            2    43   101
#            3    15   101
#            4    13   101
#            5     3   101
#            6     6   101
#            7     6   101
#            8     6   101
#            9     2   101
#           10     2   101
#           11     2   101
#           12     2   101
#           13     1   101

surv_dat <- full_data %>% 
  group_by(pid) %>% 
  slice(1) %>%
  mutate(one_fu = ifelse(total_visit == 2, 1, 0), 
         two_fu = ifelse(total_visit == 3, 1, 0), 
         three_fu = ifelse(total_visit == 4, 1, 0), 
         four_fu = ifelse(total_visit == 5, 1, 0), 
         five_fu = ifelse(total_visit == 6, 1, 0)) %>% 
  ungroup() %>% 
  select(pid, total_visit, y, age, sex, bsa, bav,
         one_fu, two_fu, three_fu, four_fu, five_fu) 

# (a) Total visit == 2 [ONE follow-up measurement]
surv_fit1 <- phreg(Surv(total_visit, one_fu) ~ y+bav+age+sex+bsa, 
                   data = surv_dat, dist="weibull")

# (b) Total visit == 3 [TWO follow-up measurement]
surv_fit2 <- phreg(Surv(total_visit, two_fu) ~ y+bav+age+sex+bsa, 
                   data = surv_dat, dist="weibull")

# (c) Total visit == 4 [THREE follow-up measurement]
surv_fit3 <- phreg(Surv(total_visit, three_fu) ~ y+bav+age+sex+bsa, 
                   data = surv_dat, dist="weibull")

# (d) Total visit == 5 [FOUR follow-up measurement]
surv_fit4 <- phreg(Surv(total_visit, four_fu) ~ y+bav+age+sex+bsa, 
                   data = surv_dat, dist="weibull")

# (e) Total visit == 6 [FIVE follow-up measurement]
surv_fit5 <- phreg(Surv(total_visit, five_fu) ~ y+bav+age+sex+bsa, 
                   data = surv_dat, dist="weibull")

surv_coefs <- data.frame(fit1 = coef(surv_fit1), 
                         fit2 = coef(surv_fit2), 
                         fit3 = coef(surv_fit3), 
                         fit4 = coef(surv_fit4), 
                         fit5 = coef(surv_fit5))

# Propensity Score Matching ----------------------------------------------------
baseline_data <- full_data %>% filter(visit == 1)
ps_model <- glm(bav ~ age + sex + bsa, family = binomial("logit"), data = baseline_data)
summary(ps_model)
ps_coefs <- coef(ps_model)
ps_match <- pairmatch(ps_model, data = baseline_data)
match_info <- summary(ps_match)
# Get matched pairs
matched_pairs <- data.frame(baseline_data, 
                            pairID = ps_match, 
                            check.rows = TRUE) %>% filter(!is.na(pairID))
# Transform the data to long format 
matched_long <- full_data %>% 
  filter(pid %in% matched_pairs$pid) %>% 
  left_join(matched_pairs %>% select(pid, pairID), by="pid")


# GEE --------------------------------------------------------------------------
## Full Covariate set ----------------------------------------------------------
### Exchangeable 
full_ex <- summary(geeglm(y ~ bav*visit+age+sex+bsa, 
                  id=pid, data = matched_long, 
                  wave = factor(visit),family = binomial('logit'), 
                  corstr = "exchangeable", scale.fix = TRUE))

### AR1 
full_ar1 <- summary(geeglm(y ~ bav*visit+age+sex+bsa, 
                           id=pid, data = matched_long, 
                           wave = factor(visit),family = binomial('logit'), 
                           corstr = "ar1", scale.fix = TRUE))

### Independent 
full_ind <- summary(geeglm(y ~ bav*visit+age+sex+bsa, 
                           id=pid, data = matched_long, 
                           wave = factor(visit),family = binomial('logit'), 
                           corstr = "independence", scale.fix = TRUE))

## Reduced Covariate set -------------------------------------------------------
### Exchangeable 
red_ex <- summary(geeglm(y ~ bav*visit, 
                          id=pid, data = matched_long, 
                          wave = factor(visit),family = binomial('logit'), 
                          corstr = "exchangeable", scale.fix = TRUE))
### AR1 
red_ar1 <- summary(geeglm(y ~ bav*visit, 
                           id=pid, data = matched_long, 
                           wave = factor(visit),family = binomial('logit'), 
                           corstr = "ar1", scale.fix = TRUE))

### Independent 
red_ind <- summary(geeglm(y ~ bav*visit, 
                           id=pid, data = matched_long, 
                           wave = factor(visit),family = binomial('logit'), 
                           corstr = "independence", scale.fix = TRUE))

# Save Coefficients 
# save(surv_coefs, ps_coefs, full_ex, full_ar1, full_ind,
#      red_ex, red_ar1, red_ind,
#      file = "data/realcoefs.RData")


# Save model fits 
# save(surv_fit1, surv_fit2, surv_fit3, surv_fit4, surv_fit5, ps_model,
#      match_info, file = "data/modelfits.RData")
# save(matched_long, file = "data/matched_samples.rds")

# QLS --------------------------------------------------------------------------
formula_full <- y ~ bav*visit+age+male+bsa
formula_red <- y ~ bav*visit
source("qls_functions.R")
qls_data <- matched_long[order(matched_long$pairID),] %>% 
  mutate(clusterID = as.integer(factor(pairID))) %>% 
  group_by(clusterID) %>% 
  mutate(cluster.var = ifelse(bav == 0, 1, 2), 
         order = row_number()) %>% 
  rename(id = pid, male = sex) %>% 
  relocate(c(clusterID, cluster.var, order), .after = id) %>% 
  filter(visit <= 6)


# Independence Correlation 
qls_ind_red <- qls(formula_red, data=qls_data, maxT = 6, 
                   time.var = qls_data$visit, corstr = "independence")
qls_ind_full <- qls(formula_full, data=qls_data, maxT = 6, 
                   time.var = qls_data$visit, corstr = "independence")

# AR1 
qls_ar1_red <- qls(formula_red, data=qls_data, time.var = qls_data$visit, maxT = 6, corstr = "ar1")
qls_ar1_full <- qls(formula_full, data=qls_data, time.var = qls_data$visit, maxT = 6, corstr = "ar1")

# Exchangeable
qls_exch_red <- qls(formula_red, data=qls_data, time.var = qls_data$visit, maxT = 6, corstr = "exchangeable")
qls_exch_full <- qls(formula_full, data=qls_data, time.var = qls_data$visit, maxT = 6, corstr = "exchangeable")

N <- n_distinct(matched_long$pid)
p <- 4 
DF <- N/(N-p)

df_se <- data.frame(
  method = c(rep("GEE", 12), rep("QLS", 12)), 
  mdl = c(rep("ind",4),  rep("ar1", 4), rep("exch", 4), 
          rep("ind",4),  rep("ar1", 4), rep("exch", 4)), 
  df_se = c(sqrt(diag(DF*red_ind$cov.scaled)), 
            sqrt(diag(DF*red_ar1$cov.scaled)), 
            sqrt(diag(DF*red_ex$cov.scaled)), 
            sqrt(diag(DF*qls_ind_red$vcov)), 
            sqrt(diag(DF*qls_ar1_red$vcov)), 
            sqrt(diag(DF*qls_exch_red$vcov)))
)

save(df_se, file = "data/se_df.RData")
save(qls_ind_red, qls_ind_full, qls_ar1_red, qls_ar1_full, qls_exch_red, qls_exch_full, file = "data/qlsfits.RData")
