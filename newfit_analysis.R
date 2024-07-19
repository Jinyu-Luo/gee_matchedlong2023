
load("~/Documents/Projects/Practicum_GEE/Outputs/simfit.RData")
source("~/Documents/Projects/Practicum_GEE/fits_analysis.R")

# Convergence Check ------------------------------------------------------------
gee_convergence <- calculate_convergence(gee_fits_df)
gee_sim_results <- gee_convergence$data

qls_convergence <- calculate_convergence(qls_fits_df)
qls_sim_results <- qls_convergence$data

divergence <- list(gee = gee_convergence$diverged, qls = qls_convergence$diverged)


# Extreme SE check -------------------------------------------------------------
gee_extremSE <- extr_se_search(gee_sim_results, se_max = 5)
qls_extremSE <- extr_se_search(qls_sim_results, se_max = 5)
cat("GEE\n", gee_extremSE$messages)
cat("QLS\n", gee_extremSE$messages)
extreme_SEs <- list(gee = gee_extremSE$extreme_id, qls = qls_extremSE$extreme_id)
gee_sim_results <- gee_extremSE$filtered_results
qls_sim_results <- qls_extremSE$filtered_results


# Coverage Probability ---------------------------------------------------------
# GEE 
# (1) Coverage probabilities for models without adjustment by age, male, and BSA
coverage_gee_unadj <- calculate_coverage(gee_sim_results$ind_mdl_red, true_set) %>%
  rbind(calculate_coverage(gee_sim_results$ar1_mdl_red, true_set)) %>%
  rbind(calculate_coverage(gee_sim_results$exch_mdl_red, true_set)) %>%
  na.omit()

# (2) Coverage probabilities for models adjusted by age, male, and BSA
coverage_gee_adj <- calculate_coverage(gee_sim_results$ind_mdl_full, true_set) %>%
  rbind(calculate_coverage(gee_sim_results$ar1_mdl_full, true_set)) %>%
  rbind(calculate_coverage(gee_sim_results$exch_mdl_full, true_set)) %>%
  na.omit()

# QLS
# (1) Coverage probabilities for models without adjustment by age, male, and BSA
coverage_qls_unadj <- calculate_coverage(qls_sim_results$ind_mdl_red, true_set) %>%
  rbind(calculate_coverage(qls_sim_results$ar1_mdl_red, true_set)) %>%
  rbind(calculate_coverage(qls_sim_results$exch_mdl_red, true_set)) %>%
  na.omit()

# (2) Coverage probabilities for models adjusted by age, male, and BSA
coverage_qls_adj <- calculate_coverage(qls_sim_results$ind_mdl_full, true_set) %>%
  rbind(calculate_coverage(qls_sim_results$ar1_mdl_full, true_set)) %>%
  rbind(calculate_coverage(qls_sim_results$exch_mdl_full, true_set)) %>%
  na.omit()

coverage <- list(gee = cbind(coverage_gee_unadj, coverage_gee_adj) %>% rbind(c(rep("NCFV", 3), rep("CFV", 3))), 
                 qls = cbind(coverage_qls_unadj, coverage_qls_adj) %>% rbind(c(rep("NCFV", 3), rep("CFV", 3))))


# Extract Data For Plots -------------------------------------------------------
bav_red <- extract_term_data(gee_sim_results, term = "bav", mod = "_red") %>% mutate(type = "GEE") %>%
  rbind(extract_term_data(qls_sim_results, term = "bav", mod = "_red") %>% mutate(type = "QLS")) 

visit_red <- extract_term_data(gee_sim_results, term = "visit", mod = "_red") %>% mutate(type = "GEE") %>%
  rbind(extract_term_data(qls_sim_results, term = "visit", mod = "_red") %>% mutate(type = "QLS")) 

bav_visit_red <- extract_term_data(gee_sim_results, term = "bav:visit", mod = "_red") %>% mutate(type = "GEE") %>%
  rbind(extract_term_data(qls_sim_results, term = "bav:visit", mod = "_red") %>% mutate(type = "QLS")) 

red_mdl_df <- rbind(bav_red, visit_red, bav_visit_red) %>% 
  relocate(type, .after = sim_id)

red_summary <- red_mdl_df %>%
  mutate(term = case_when(term == "bav" ~ "BAV",
                          term == "visit" ~ "Visit",
                          TRUE ~ "BAV:Visit")) %>%
  group_by(term, model, type) %>%
  summarise(mean_estimate = mean(estimate),
            mean_lower = mean(lower),
            mean_upper = mean(upper),
            mean_adj_lower = mean(adj_lower),
            mean_adj_upper = mean(adj_upper),
            .groups = 'drop') %>%
  pivot_longer(cols = c(mean_lower, mean_upper, mean_adj_lower, mean_adj_upper),
               names_to = "ci_type", values_to = "ci_value") %>%
  mutate(ci_type = case_when(
    ci_type == "mean_lower" ~ "Unadjusted Lower",
    ci_type == "mean_upper" ~ "Unadjusted Upper",
    ci_type == "mean_adj_lower" ~ "Adjusted Lower",
    ci_type == "mean_adj_upper" ~ "Adjusted Upper")) %>%
  pivot_wider(names_from = ci_type, values_from = ci_value) %>%
  pivot_longer(cols = c("Unadjusted Lower", "Unadjusted Upper", "Adjusted Lower", "Adjusted Upper"),
               names_to = "ci_type", values_to = "ci_value") %>%
  separate(ci_type, into = c("adjusted", "boundary"), sep = " ") %>%
  pivot_wider(names_from = boundary, values_from = ci_value)

bav_full <- extract_term_data(gee_sim_results, term = "bav", mod = "_full") %>% mutate(type = "GEE") %>%
  rbind(extract_term_data(qls_sim_results, term = "bav", mod = "_full") %>% mutate(type = "QLS")) %>%
  relocate(type, .after = sim_id)

visit_full <- extract_term_data(gee_sim_results, term = "visit", mod = "_full") %>% mutate(type = "GEE") %>%
  rbind(extract_term_data(qls_sim_results, term = "visit", mod = "_full") %>% mutate(type = "QLS")) %>%
  relocate(type, .after = sim_id)

bav_visit_full <- extract_term_data(gee_sim_results, term = "bav:visit", mod = "_full") %>% mutate(type = "GEE") %>%
  rbind(extract_term_data(qls_sim_results, term = "bav:visit", mod = "_full") %>% mutate(type = "QLS")) %>%
  relocate(type, .after = sim_id)

age_df <- extract_term_data(gee_sim_results, term = "age", mod = "_full") %>% mutate(type = "GEE") %>%
  rbind(extract_term_data(qls_sim_results, term = "age", mod = "_full") %>% mutate(type = "QLS")) %>%
  relocate(type, .after = sim_id)

male_df <- extract_term_data(gee_sim_results, term = "male", mod = "_full") %>% mutate(type = "GEE") %>%
  rbind(extract_term_data(qls_sim_results, term = "male", mod = "_full") %>% mutate(type = "QLS")) %>%
  relocate(type, .after = sim_id)

bsa_df <- extract_term_data(gee_sim_results, term = "bsa", mod = "_full") %>% mutate(type = "GEE") %>%
  rbind(extract_term_data(qls_sim_results, term = "bsa", mod = "_full") %>% mutate(type = "QLS")) %>%
  relocate(type, .after = sim_id)

full_mdl_df <- rbind(bav_full, visit_full, bav_visit_full, age_df, male_df, bsa_df)

full_summary <- full_mdl_df %>%
  mutate(term = case_when(term == "bav" ~ "BAV",
                          term == "visit" ~ "Visit",
                          term == "age" ~ "Age",
                          term == "male" ~ "Male",
                          term == "bsa" ~ "BSA",
                          TRUE ~ "BAV:Visit")) %>%
  group_by(term, model, type) %>%
  summarise(mean_estimate = mean(estimate),
            mean_lower = mean(lower),
            mean_upper = mean(upper),
            mean_adj_lower = mean(adj_lower),
            mean_adj_upper = mean(adj_upper),
            .groups = 'drop') %>%
  pivot_longer(cols = c(mean_lower, mean_upper, mean_adj_lower, mean_adj_upper),
               names_to = "ci_type", values_to = "ci_value") %>%
  mutate(ci_type = case_when(
    ci_type == "mean_lower" ~ "Unadjusted Lower",
    ci_type == "mean_upper" ~ "Unadjusted Upper",
    ci_type == "mean_adj_lower" ~ "Adjusted Lower",
    ci_type == "mean_adj_upper" ~ "Adjusted Upper")) %>%
  pivot_wider(names_from = ci_type, values_from = ci_value) %>%
  pivot_longer(cols = c("Unadjusted Lower", "Unadjusted Upper", "Adjusted Lower", "Adjusted Upper"),
               names_to = "ci_type", values_to = "ci_value") %>%
  separate(ci_type, into = c("adjusted", "boundary"), sep = " ") %>%
  pivot_wider(names_from = boundary, values_from = ci_value)

# Relative Bias ----------------------------------------------------------------
red_relbias_df <- cal_rel_bias(gee_sim_results$ind_mdl_red, true_set) %>% 
  cbind(method = "GEE") %>% cbind(corstr = "Independence") %>% 
  rbind(cal_rel_bias(gee_sim_results$ar1_mdl_red, true_set) %>% 
          cbind(method = "GEE") %>% cbind(corstr = "AR1")) %>% 
  rbind(cal_rel_bias(gee_sim_results$exch_mdl_red, true_set) %>% 
          cbind(method = "GEE") %>% cbind(corstr = "Exchangeable")) %>% 
  rbind(cal_rel_bias(qls_sim_results$ind_mdl_red, true_set) %>% 
          cbind(method = "QLS") %>% cbind(corstr = "Independence")) %>% 
  rbind(cal_rel_bias(qls_sim_results$ar1_mdl_red, true_set) %>% 
          cbind(method = "QLS") %>% cbind(corstr = "AR1")) %>% 
  rbind(cal_rel_bias(qls_sim_results$exch_mdl_red, true_set) %>% 
          cbind(method = "QLS") %>% cbind(corstr = "Exchangeable")) %>% 
  cbind(cov_set = "No")

full_relbias_df <- cal_rel_bias(gee_sim_results$ind_mdl_full, true_set) %>% 
  cbind(method = "GEE") %>% cbind(corstr = "Independence") %>% 
  rbind(cal_rel_bias(gee_sim_results$ar1_mdl_full, true_set) %>% 
          cbind(method = "GEE") %>% cbind(corstr = "AR1")) %>% 
  rbind(cal_rel_bias(gee_sim_results$exch_mdl_full, true_set) %>% 
          cbind(method = "GEE") %>% cbind(corstr = "Exchangeable")) %>% 
  rbind(cal_rel_bias(qls_sim_results$ind_mdl_full, true_set) %>% 
          cbind(method = "QLS") %>% cbind(corstr = "Independence")) %>% 
  rbind(cal_rel_bias(qls_sim_results$ar1_mdl_full, true_set) %>% 
          cbind(method = "QLS") %>% cbind(corstr = "AR1")) %>% 
  rbind(cal_rel_bias(qls_sim_results$exch_mdl_full, true_set) %>% 
          cbind(method = "QLS") %>% cbind(corstr = "Exchangeable")) %>% 
  cbind(cov_set = "Yes")

relbias_df <- rbind(red_relbias_df, full_relbias_df) %>% 
  relocate(c(method, corstr, cov_set), .before = b0)

# Mean Rho ---------------------------------------------------------------------
rho_df <- calculate_mean_rho(gee_sim_results$ar1_mdl_red) %>% 
  rbind(calculate_mean_rho(gee_sim_results$ar1_mdl_full)) %>% 
  rbind(calculate_mean_rho(gee_sim_results$exch_mdl_red)) %>% 
  rbind(calculate_mean_rho(gee_sim_results$exch_mdl_full)) %>% 
  select(mean_rho, bias) %>% 
  cbind(corstr = c(rep("AR1", 2), rep("Exchangeable", 2))) %>% 
  cbind(mod = rep(c("No", "Yes"), 2)) %>% 
  cbind(rho=0.3) %>% 
  relocate(c(corstr, mod, rho), .before = mean_rho) %>% 
  cbind(calculate_mean_rho(qls_sim_results$ar1_mdl_red) %>% 
          rbind(calculate_mean_rho(qls_sim_results$ar1_mdl_full)) %>%
          rbind(calculate_mean_rho(qls_sim_results$exch_mdl_red)) %>% 
          rbind(calculate_mean_rho(qls_sim_results$ar1_mdl_full)) %>% 
          select(mean_rho, bias))

# Tau --------------------------------------------------------------------------
red_tau <- data.frame(
  ind = pull_tau(qls_sim_results$ind_mdl_red),
  ar1 = pull_tau(qls_sim_results$ar1_mdl_red),
  exch = pull_tau(qls_sim_results$exch_mdl_red)) %>% 
  cbind(covset = "No") 

full_tau <- data.frame(
  ind = pull_tau(qls_sim_results$ind_mdl_full),
  ar1 = pull_tau(qls_sim_results$ar1_mdl_full),
  exch = pull_tau(qls_sim_results$exch_mdl_full)) %>% 
  cbind(covset = "Yes")

tau_df <- red_tau %>% rbind(full_tau) %>% 
  pivot_longer(!covset, names_to = "corstr", values_to = "tau") 


save(gee_sim_results, qls_sim_results, divergence, coverage, 
     red_summary, full_summary, relbias_df, rho_df, tau_df,
     file = "Outputs/newfits_analysis.RData")
