#Tables
# ============================================================
# TABLE 2
# Level-specific performance by population model and method
# ============================================================

# Truth label for each population model
pop_truth <- eval_main %>%
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


# Aggregate cell-level rejection rates
table2_long <- perf_cell %>%
  filter(moderator != "noise") %>%
  group_by(
    popmodel,
    method,
    level,
    estimand
  ) %>%
  summarise(
    n_cells = n(),
    
    mean_rate = mean(
      rate,
      na.rm = TRUE
    ),
    
    # Monte Carlo SE of the equally weighted mean
    mcse_mean =
      sqrt(
        sum(mcse^2, na.rm = TRUE)
      ) / n_cells,
    
    min_rate = min(
      rate,
      na.rm = TRUE
    ),
    
    max_rate = max(
      rate,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  left_join(
    pop_truth %>%
      select(popmodel, truth),
    by = "popmodel"
  )

table2_long
#---------------------------------------------------------------------------
table2 <- table2_long %>%
  select(
    popmodel,
    truth,
    method,
    level,
    mean_rate,
    mcse_mean
  ) %>%
  pivot_wider(
    names_from = level,
    values_from = c(
      mean_rate,
      mcse_mean
    )
  ) %>%
  rename(
    metric_rejection = mean_rate_Metric,
    scalar_rejection = mean_rate_Scalar,
    metric_mcse = mcse_mean_Metric,
    scalar_mcse = mcse_mean_Scalar
  ) %>%
  arrange(
    factor(
      popmodel,
      levels = c(
        "0", "1.1", "1.11", "1.12",
        "1.2", "1.21", "1.22",
        "1.3", "1.32"
      )
    ),
    method
  )

table2
#------------
table2_print <- table2 %>%
  mutate(
    Metric = sprintf(
      "%.3f (%.3f)",
      metric_rejection,
      metric_mcse
    ),
    
    Scalar = sprintf(
      "%.3f (%.3f)",
      scalar_rejection,
      scalar_mcse
    )
  ) %>%
  select(
    Population_Model = popmodel,
    Truth = truth,
    Method = method,
    Metric,
    Scalar
  )

table2_print
#----------------------------------------------------------------------------
# ============================================================
# TABLE 3
# Performance by moderator functional form
# ============================================================

table3_long <- perf_cell %>%
  filter(
    moderator != "noise"
  ) %>%
  group_by(
    level,
    estimand,
    moderator,
    method
  ) %>%
  summarise(
    n_cells = n(),
    
    mean_rate = mean(
      rate,
      na.rm = TRUE
    ),
    
    mcse_mean =
      sqrt(
        sum(mcse^2, na.rm = TRUE)
      ) / n_cells,
    
    min_rate = min(
      rate,
      na.rm = TRUE
    ),
    
    max_rate = max(
      rate,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

table3_long
#--------------
table3 <- table3_long %>%
  mutate(
    result = sprintf(
      "%.3f (%.3f)",
      mean_rate,
      mcse_mean
    )
  ) %>%
  select(
    level,
    estimand,
    moderator,
    method,
    result
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = result
  ) %>%
  arrange(
    factor(
      level,
      levels = c("Metric", "Scalar")
    ),
    factor(
      estimand,
      levels = c(
        "Power",
        "Rejection u. invariance"
      )
    ),
    factor(
      moderator,
      levels = c(
        "linear",
        "quadratic",
        "sigmoid"
      )
    )
  )

table3
#------------------------
table3_noise <- perf_cell %>%
  filter(
    moderator == "noise"
  ) %>%
  group_by(
    level,
    method
  ) %>%
  summarise(
    n_cells = n(),
    
    mean_rejection = mean(
      rate,
      na.rm = TRUE
    ),
    
    mcse_mean =
      sqrt(
        sum(mcse^2, na.rm = TRUE)
      ) / n_cells,
    
    .groups = "drop"
  )

table3_noise

# ============================================================
# TABLE 4
# Sequential performance by true MNI pattern
# ============================================================

table4 <- perf_sequential_cell %>%
  filter(
    moderator != "noise"
  ) %>%
  group_by(
    truth_class,
    method
  ) %>%
  summarise(
    n_cells = n(),
    
    metric_decision =
      mean(
        p_metric_decision,
        na.rm = TRUE
      ),
    
    scalar_stage_reached =
      mean(
        scalar_stage_reached,
        na.rm = TRUE
      ),
    
    scalar_conditional_rejection =
      mean(
        scalar_conditional_rejection,
        na.rm = TRUE
      ),
    
    scalar_decision =
      mean(
        p_scalar_decision,
        na.rm = TRUE
      ),
    
    invariance_retained =
      mean(
        p_invariance_retained,
        na.rm = TRUE
      ),
    
    sequential_accuracy =
      mean(
        sequential_accuracy,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  arrange(
    factor(
      truth_class,
      levels = c(
        "No MNI",
        "Metric only",
        "Scalar only",
        "Metric + scalar"
      )
    ),
    method
  )

table4
#-------------------
table4_print <- table4 %>%
  mutate(
    across(
      c(
        metric_decision,
        scalar_stage_reached,
        scalar_decision,
        invariance_retained,
        sequential_accuracy
      ),
      ~ round(.x, 3)
    )
  ) %>%
  select(
    Truth = truth_class,
    Method = method,
    `Metric decision` = metric_decision,
    `Scalar stage reached` = scalar_stage_reached,
    `Scalar decision` = scalar_decision,
    `Invariance retained` = invariance_retained,
    `Correct classification` = sequential_accuracy
  )

table4_print

#------------------------------------------------------------------------------
print(n=27,table2_print)
table3
table4_print
#------------------------------------------------------------------------------