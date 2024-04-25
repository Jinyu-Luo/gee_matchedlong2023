# ============= Data Cleaning =============== # 

# Load Data 
load("1 - data setup.RData")

library(tidyverse)
library(dplyr)
library(optmatch)
library(simsurv)
library(survival)
library(eha)
# Step 1: Merge raw_long_d and work_pt_d
new_long <- raw_long_d %>%
  left_join(select(work_pt_d, ptid, diabetes, hyper, chlstrl, bsa), by = "ptid") %>% 
  select(ptid, age, sex, bav_confirmed, bsa, died, root,
         mdate, lka_d, dateor, yr2vst, day2vst, yr2death, date_base)

# Step 2: Data cleaning 
new_data <- new_long %>% 
  mutate(bav = ifelse(bav_confirmed == "0_TAV", 0, 1),
         bav = factor(bav, labels = c("TAV", "BAV")), 
         death = factor(died, levels = c(0, 1), labels = c("Survived", "Died")), 
         sex = ifelse(sex == "1_male", 1, 0), 
         sex = factor(sex, levels = c(0, 1), labels = c("female", "male")), 
         across(c(dateor, mdate, date_base), ~as.Date(., format = "%m/%d/%Y"))) %>% 
  filter(!is.na(root)) %>% # remove records without root size measurement 
  mutate(vst_year = year(mdate)) %>%  # extract the year of the current visit
  group_by(ptid) %>% 
  mutate(visit = row_number(), total_visit = n(), 
         yr2vst_diff = c(0, diff(yr2vst)),
         root_growth = root - lag(root), 
         vst_1yr = first(vst_year), 
         yr_diff = vst_year - vst_1yr, 
         age_adj = age+yr_diff)%>%   
  select(-bav_confirmed, -yr_diff) %>% 
  filter(min(yr2vst)==0 & total_visit <= 8) %>% 
  filter(total_visit >= 2) %>% 
  ungroup() 

new_data <- new_data %>% group_by(ptid) %>% 
  mutate(y = case_when(root > 45 ~ 1, 
                       root < 45 & is.na(root_growth) ~ 0, 
                       root < 45 & root_growth > 5 ~ 1, 
                       TRUE ~ 0)) 

baseline_df <- new_data %>% group_by(ptid) %>% slice(1) %>% 
  mutate(bav = ifelse(bav == "BAV", 1, 0), 
         sex = ifelse(sex == "male", 1, 0))

ps_model <- glm(bav ~ age + sex + bsa, family = binomial("logit"), data = baseline_df)
summary(ps_model)


# Perform propensity score matching 
pps_match <- pairmatch(ps_model, data = baseline_df)
summary(pps_match) # 64 matched pairs and 95 unmatched patients


# Extract matched pairs 
matched_pairs <- data.frame(baseline_df,  pairID = pps_match, check.rows = TRUE) %>% filter(!is.na(pairID))

library(cobalt)
# Check Balance 
covs <- subset(baseline_df, select = c(age, sex, bsa))
pps <- ps_model$fitted.values 
bal.tab(pps_match, covs = covs, distance = ~ pps)

# Get the matched data in long format 
matched_long <- new_data %>% filter(ptid %in% matched_pairs$ptid) %>% 
  filter(visit <= 9) %>% 
  mutate(id = substr(ptid, 5, nchar(ptid))) %>% 
  left_join(matched_pairs %>% select(ptid, pairID), by = "ptid")

time_to_event <- matched_long %>% 
  mutate(Event = ifelse(total_visit < 8, 1, 0)) %>% 
  distinct(ptid, .keep_all = TRUE) %>%
  group_by(ptid) %>%
  summarize(Time = first(total_visit),
            y = first(y),
            Event = max(Event))

surv_data <- time_to_event %>% 
  left_join(matched_long %>% select(-y), by = "ptid") %>% 
  select(ptid, visit, total_visit, y, age, sex, bav, bsa, Time, Event) 

glimpse(surv_data)

cox_model <- coxph(Surv(Time, Event) ~ age + bav+ y+sex+bsa, data = surv_data)
cox_model

# Call:
#   coxph(formula = Surv(Time, Event) ~ age + bav + y + sex + bsa, 
#         data = surv_data)
# 
# n= 137, number of events= 121 
# 
#             coef exp(coef) se(coef)      z Pr(>|z|)    
# age     -0.03810   0.96262  0.01086 -3.507 0.000454 ***
# bavBAV   0.16574   1.18027  0.19779  0.838 0.402052    
# y        0.93280   2.54161  0.63862  1.461 0.144110    
# sexmale -0.67235   0.51051  0.20912 -3.215 0.001304 ** 
# bsa     -0.52135   0.59372  0.68091 -0.766 0.443876    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
#         exp(coef) exp(-coef) lower .95 upper .95
# age        0.9626     1.0388    0.9423    0.9833
# bavBAV     1.1803     0.8473    0.8010    1.7392
# y          2.5416     0.3935    0.7270    8.8859
# sexmale    0.5105     1.9588    0.3388    0.7691
# bsa        0.5937     1.6843    0.1563    2.2551
# 
# Concordance= 0.665  (se = 0.034 )
# Likelihood ratio test= 33.79  on 5 df,   p=3e-06
# Wald test            = 37.5  on 5 df,   p=5e-07
# Score (logrank) test = 41.14  on 5 df,   p=9e-08

weibull_mod <- phreg(Surv(Time, Event) ~ age + bav+ y+sex+bsa, data = surv_data, 
                     dist = "weibull", shape = 1)
weibull_mod
# Call:
#   phreg(formula = Surv(Time, Event) ~ age + bav + y + sex + bsa, 
#         data = surv_data, dist = "weibull", shape = 1)
# 
# Covariate          W.mean      Coef Exp(Coef)  se(Coef)    Wald p
# age                72.491    -0.019     0.981     0.010     0.055 
# bav 
#              TAV    0.567     0         1           (reference)
#              BAV    0.433     0.159     1.172     0.188     0.398 
# y                   0.026     0.654     1.924     0.602     0.277 
# sex 
#           female    0.424     0         1           (reference)
#             male    0.576    -0.462     0.630     0.203     0.022 
# bsa                 1.818    -0.456     0.634     0.585     0.436 
# 
# log(scale)                   -0.695               1.127     0.538 
# 
# Shape is fixed at  1 
# 
# Events                    121 
# Total time at risk           609 
# Max. log. likelihood      -308.97 
# LR test statistic         15.13 
# Degrees of freedom        5 
# Overall p-value           0.0098103

basehaz(cox_model, centered = FALSE)

coxy <- coxph(Surv(Time, Event) ~ y, data = surv_data)
summary(coxy)














