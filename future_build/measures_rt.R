# measures_rt.R

library(tidyverse)

results_path <- "results.csv"

results <- read.csv(results_path, stringsAsFactors = FALSE) %>%
  tibble::as_tibble()

# ---------------------------
# overall metric performance
# ---------------------------
metric_summary <- results %>%
  summarise(
    mnlfa_metric_type1 = mean(mnlfa_metric_lrt_reject[true_metric_noninvariance == FALSE], na.rm = TRUE),
    mnlfa_metric_power = mean(mnlfa_metric_lrt_reject[true_metric_noninvariance == TRUE], na.rm = TRUE),
    
    tree_metric_type1  = mean(tree_metric_reject[true_metric_noninvariance == FALSE], na.rm = TRUE),
    tree_metric_power  = mean(tree_metric_reject[true_metric_noninvariance == TRUE], na.rm = TRUE)
  )

# ---------------------------
# overall scalar performance
# ---------------------------
scalar_summary <- results %>%
  summarise(
    mnlfa_scalar_type1 = mean(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == FALSE], na.rm = TRUE),
    mnlfa_scalar_power = mean(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == TRUE], na.rm = TRUE),
    
    tree_scalar_type1  = mean(tree_scalar_reject[true_scalar_noninvariance == FALSE], na.rm = TRUE),
    tree_scalar_power  = mean(tree_scalar_reject[true_scalar_noninvariance == TRUE], na.rm = TRUE)
  )

# ---------------------------
# stratified performance
# ---------------------------
condition_summary <- results %>%
  group_by(N, reliability, moderator) %>%
  summarise(
    n_replications = n(),
    
    mnlfa_metric_type1 = mean(mnlfa_metric_lrt_reject[true_metric_noninvariance == FALSE], na.rm = TRUE),
    mnlfa_metric_power = mean(mnlfa_metric_lrt_reject[true_metric_noninvariance == TRUE], na.rm = TRUE),
    
    tree_metric_type1  = mean(tree_metric_reject[true_metric_noninvariance == FALSE], na.rm = TRUE),
    tree_metric_power  = mean(tree_metric_reject[true_metric_noninvariance == TRUE], na.rm = TRUE),
    
    mnlfa_scalar_type1 = mean(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == FALSE], na.rm = TRUE),
    mnlfa_scalar_power = mean(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == TRUE], na.rm = TRUE),
    
    tree_scalar_type1  = mean(tree_scalar_reject[true_scalar_noninvariance == FALSE], na.rm = TRUE),
    tree_scalar_power  = mean(tree_scalar_reject[true_scalar_noninvariance == TRUE], na.rm = TRUE),
    
    tree_scalar_available = mean(!is.na(tree_scalar_reject)),
    
    .groups = "drop"
  )

print(metric_summary)
print(scalar_summary)
print(condition_summary)
