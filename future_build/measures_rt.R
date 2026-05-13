# measures_rt.R

library(dplyr)
library(tidyr)
library(tibble)

results_path <- "results.csv"

results <- read.csv(results_path, stringsAsFactors = FALSE) %>%
  tibble::as_tibble()

# ---------------------------
# helper functions
# ---------------------------
rate_by_truth <- function(decision, truth) {
  tibble(
    type1 = mean(decision[truth == FALSE], na.rm = TRUE),
    power = mean(decision[truth == TRUE],  na.rm = TRUE),
    n_type1 = sum(!is.na(decision[truth == FALSE])),
    n_power = sum(!is.na(decision[truth == TRUE]))
  )
}

rate_tree_by_truth <- function(reject, correct_split, truth) {
  detection <- reject == TRUE
  correct_detection <- reject == TRUE & correct_split == TRUE
  
  tibble(
    type1_detection = mean(detection[truth == FALSE], na.rm = TRUE),
    power_detection = mean(detection[truth == TRUE], na.rm = TRUE),
    
    type1_wrong_split = mean(detection[truth == TRUE & correct_split == FALSE], na.rm = TRUE),
    power_correct_split = mean(correct_detection[truth == TRUE], na.rm = TRUE),
    
    n_type1 = sum(!is.na(reject[truth == FALSE])),
    n_power = sum(!is.na(reject[truth == TRUE]))
  )
}
# ---------------------------
# overall Type I error and power
# ---------------------------
overall_performance <- bind_rows(
  rate_by_truth(results$mnlfa_metric_lrt_reject, results$true_metric_noninvariance) %>%
    mutate(method = "MNLFA", level = "metric"),
  
  rate_tree_by_truth(results$tree_metric_reject, results$tree_metric_correct_split, results$true_metric_noninvariance) %>%
    transmute(
      type1 = type1_detection,
      power = power_correct_split,
      n_type1,
      n_power,
      method = "SEMTREE",
      level = "metric"
    ),
  
  rate_by_truth(results$mnlfa_scalar_lrt_reject, results$true_scalar_noninvariance) %>%
    mutate(method = "MNLFA", level = "scalar"),
  
  rate_tree_by_truth(results$tree_scalar_reject, results$tree_scalar_correct_split, results$true_scalar_noninvariance) %>%
    transmute(
      type1 = type1_detection,
      power = power_correct_split,
      n_type1,
      n_power,
      method = "SEMTREE",
      level = "scalar"
    )
) %>%
  select(level, method, type1, power, n_type1, n_power)

print(overall_performance)

# ---------------------------
# condition-specific Type I error and power
# ---------------------------
condition_performance <- results %>%
  group_by(
    N, reliability, popmodel, moderator, analysis_form,
    delta_lambda, delta_nu
  ) %>%
  summarise(
    n_rows = n(),
    
    mnlfa_metric_type1 = mean(mnlfa_metric_lrt_reject[true_metric_noninvariance == FALSE], na.rm = TRUE),
    mnlfa_metric_power = mean(mnlfa_metric_lrt_reject[true_metric_noninvariance == TRUE],  na.rm = TRUE),
    
    tree_metric_type1 = mean(tree_metric_reject[true_metric_noninvariance == FALSE], na.rm = TRUE),
    tree_metric_power = mean(
      (tree_metric_reject == TRUE & tree_metric_correct_split == TRUE)[true_metric_noninvariance == TRUE],
      na.rm = TRUE
    ),
    
    mnlfa_scalar_type1 = mean(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == FALSE], na.rm = TRUE),
    mnlfa_scalar_power = mean(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == TRUE],  na.rm = TRUE),
    
    tree_scalar_type1 = mean(tree_scalar_reject[true_scalar_noninvariance == FALSE], na.rm = TRUE),
    tree_scalar_power = mean(
      (tree_scalar_reject == TRUE & tree_scalar_correct_split == TRUE)[true_scalar_noninvariance == TRUE],
      na.rm = TRUE
    ),
    
    n_metric_type1 = sum(!is.na(mnlfa_metric_lrt_reject[true_metric_noninvariance == FALSE])),
    n_metric_power = sum(!is.na(mnlfa_metric_lrt_reject[true_metric_noninvariance == TRUE])),
    n_scalar_type1 = sum(!is.na(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == FALSE])),
    n_scalar_power = sum(!is.na(mnlfa_scalar_lrt_reject[true_scalar_noninvariance == TRUE])),
    
    .groups = "drop"
  )

print(condition_performance)

condition_performance_long <- condition_performance %>%
  pivot_longer(
    cols = c(
      mnlfa_metric_type1,
      mnlfa_metric_power,
      tree_metric_type1,
      tree_metric_power,
      mnlfa_scalar_type1,
      mnlfa_scalar_power,
      tree_scalar_type1,
      tree_scalar_power
    ),
    names_to = "measure",
    values_to = "rate"
  ) %>%
  mutate(
    method = ifelse(grepl("^mnlfa", measure), "MNLFA", "SEMTREE"),
    level = ifelse(grepl("metric", measure), "metric", "scalar"),
    estimand = ifelse(grepl("type1", measure), "type1", "power")
  ) %>%
  select(
    N, reliability, popmodel, moderator, analysis_form,
    delta_lambda, delta_nu,
    method, level, estimand, rate
  )

print(condition_performance_long)