# title: Propensity score matching
# author: Jinyu Luo
# version: 2024-07-18
rm(list = ls())

library(optmatch)
library(tidyverse)
library(lme4)

load("~/Documents/Projects/Practicum_GEE/Outputs/sim_fulldata.RData")

ps_match <- function(dat){
  base_info <- dat %>% filter(visit == 1) %>% select(-total_visit)
  ps_model <- glm(bav ~ age + male + bsa, family = binomial, data = base_info)
  pps_match <- pairmatch(ps_model, data = base_info)
  matched_df <- data.frame(base_info, matched = pps_match, check.rows = TRUE) %>% 
    filter(!is.na(matched))
  matchid <- matched_df %>% select(id, matched)
  finaldata <- dat %>% right_join(matchid, by = "id")
  finaldata
}

matched_df <- sim_df %>% 
  mutate(matched_full =  map(full_data, ~ps_match(.x))) %>% 
  mutate(matched_ddat = map(dropout_data, ~ps_match(.x))) %>% 
  matched_df %>% select(matched_full, matched_ddat)

save(matched_df, file = "Outputs/matched_data.RData")
