#Tables
library(dplyr)
library(tidyr)
library(knitr)

# Truth label for each population model
pop_truth <- eval_confirmatory %>%
  filter(moderator != "noise") %>%
  group_by(popmodel) %>%
  summarise(
    true_metric = first(true_metric_noninvariance),
    true_scalar = first(true_scalar_noninvariance),
    .groups = "drop"
  ) %>%
  mutate(
    truth = case_when(
      !true_metric & !true_scalar ~ "No MNI",
      true_metric & !true_scalar  ~ "Metric only",
      !true_metric & true_scalar  ~ "Scalar only",
      true_metric & true_scalar   ~ "Metric + scalar",
      TRUE ~ NA_character_
    )
  )
#-------------------------------------------------------------------------------

table2_cell <- perf_cell_all %>%
  filter(moderator != "noise") %>%
  mutate(
    method = recode(
      method,
      MNLFA = "Linear MNLFA",
      MNLFAQ = "Quadratic MNLFA",
      SEMTREE = "SEM tree"
    ),
    result = if_else(
      n_eval > 0,
      sprintf("%.3f (%.3f)", rate, mcse),
      "NA"
    )
  ) %>%
  select(
    popmodel, N, moderator,
    delta_lambda, delta_nu,
    method, level, result
  ) %>%
  pivot_wider(
    names_from = c(level, method),
    values_from = result,
    names_glue = "{level}: {method}"
  ) %>%
  arrange(
    popmodel, N, moderator,
    delta_lambda, delta_nu
  )

# did the pivot produced exactly one row
# per distinct data-generating condition.
n_expected <- perf_cell_all %>%
  filter(moderator != "noise") %>%
  distinct(
    popmodel, N, moderator,
    delta_lambda, delta_nu
  ) %>%
  nrow()

stopifnot(nrow(table2_cell) == n_expected)

saveRDS(
  table2_cell,
  "results/table2_condition_specific.rds"
)

#-------------------------------------------------------------------------------
table3_cell <- perf_cell_all %>%
  filter(
    moderator %in% c("linear", "quadratic", "sigmoid"),
    true_noninvariance == TRUE,
    estimand == "Power"
  ) %>%
  mutate(
    method = recode(
      method,
      MNLFA = "Linear MNLFA",
      MNLFAQ = "Quadratic MNLFA",
      SEMTREE = "SEM tree"
    ),
    result = if_else(
      n_eval > 0,
      sprintf("%.3f (%.3f)", rate, mcse),
      "NA"
    )
  ) %>%
  select(
    popmodel, N, moderator,
    delta_lambda, delta_nu,
    level, method, result
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = result
  ) %>%
  arrange(
    popmodel, level, N,
    moderator, delta_lambda, delta_nu
  )

table3_cell

saveRDS(
  table3_cell,
  "results/table3_condition_specific.rds"
)
#-------------------------------------------------------------------------------
# FULL-DESIGN SEQUENTIAL CLASSIFICATION

eval_sequential_all <- eval_main %>%
  filter(moderator != "noise") %>%
  mutate(
    truth_class = case_when(
      !true_metric_noninvariance & !true_scalar_noninvariance ~ "No MNI",
      true_metric_noninvariance & !true_scalar_noninvariance  ~ "Metric only",
      !true_metric_noninvariance & true_scalar_noninvariance  ~ "Scalar only",
      true_metric_noninvariance & true_scalar_noninvariance   ~ "Metric + scalar",
      TRUE ~ NA_character_
    ),
    
    # Sequential procedure:
    # metric test first; scalar test only if metric is retained
    seq_decision = case_when(
      is.na(metric_reject) ~ NA_character_,
      metric_reject == TRUE ~ "Metric noninvariance",
      is.na(scalar_reject) ~ NA_character_,
      scalar_reject == TRUE ~ "Scalar noninvariance",
      TRUE ~ "Invariance retained"
    ),
    
    true_seq_decision = case_when(
      true_metric_noninvariance ~ "Metric noninvariance",
      true_scalar_noninvariance ~ "Scalar noninvariance",
      !true_metric_noninvariance & !true_scalar_noninvariance ~
        "Invariance retained",
      TRUE ~ NA_character_
    ),
    
    seq_correct = if_else(
      !is.na(seq_decision) & !is.na(true_seq_decision),
      seq_decision == true_seq_decision,
      NA
    )
  )

perf_sequential_all <- eval_sequential_all %>%
  group_by(
    method,
    popmodel,
    N,
    moderator,
    delta_lambda,
    delta_nu,
    truth_class
  ) %>%
  summarise(
    n_total = n(),
    
    n_seq_eval = sum(!is.na(seq_decision)),
    
    p_metric_decision = if_else(
      n_seq_eval > 0,
      sum(seq_decision == "Metric noninvariance", na.rm = TRUE) /
        n_seq_eval,
      NA_real_
    ),
    
    p_scalar_decision = if_else(
      n_seq_eval > 0,
      sum(seq_decision == "Scalar noninvariance", na.rm = TRUE) /
        n_seq_eval,
      NA_real_
    ),
    
    p_invariance_retained = if_else(
      n_seq_eval > 0,
      sum(seq_decision == "Invariance retained", na.rm = TRUE) /
        n_seq_eval,
      NA_real_
    ),
    
    sequential_accuracy = if_else(
      n_seq_eval > 0,
      sum(seq_correct, na.rm = TRUE) / n_seq_eval,
      NA_real_
    ),
    
    .groups = "drop"
  )


table4_cell <- perf_sequential_all %>%
  mutate(
    method = recode(
      method,
      MNLFA = "Linear MNLFA",
      MNLFAQ = "Quadratic MNLFA",
      SEMTREE = "SEM tree"
    )
  ) %>%
  transmute(
    popmodel,
    truth_class,
    N,
    moderator,
    delta_lambda,
    delta_nu,
    method,
    n_seq_eval,
    `Metric decision` = sprintf("%.3f", p_metric_decision),
    `Scalar decision` = sprintf("%.3f", p_scalar_decision),
    `Invariance retained` = sprintf("%.3f", p_invariance_retained),
    `Correct classification` = sprintf("%.3f", sequential_accuracy)
  ) %>%
  arrange(
    popmodel, N, moderator,
    delta_lambda, delta_nu, method
  )

table4_cell

saveRDS(
  table4_cell,
  "results/table4_condition_specific.rds"
)

#-------------------------------------------------------------------------------
save(
  table2_cell,
  table3_cell,
  table4_cell,
  file = "manuscript_tables.RData"
)
