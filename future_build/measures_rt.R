library(dplyr)
library(tidyr)

# ============================================================
# 1. Create common method-independent result variables
# ============================================================

eval_dat <- results %>%
  mutate(
    
    metric_reject = case_when(
      method %in% c("MNLFA", "MNLFAQ") ~
        mnlfa_metric_lrt_reject,
      method == "SEMTREE" ~
        tree_metric_reject,
      TRUE ~ NA
    ),
    
    scalar_reject = case_when(
      method %in% c("MNLFA", "MNLFAQ") ~
        mnlfa_scalar_lrt_reject,
      method == "SEMTREE" ~
        tree_scalar_reject,
      TRUE ~ NA
    )
  )


# ============================================================
# 2. Separate main simulation and noisy-predictor sensitivity
# ============================================================

eval_main <- eval_dat %>%
  filter(num_noisy_predictors == 0)

eval_noise_sensitivity <- eval_dat %>%
  filter(num_noisy_predictors > 0)


# ============================================================
# 3. Define effective effect magnitudes
# ============================================================

eval_main <- eval_main %>%
  mutate(
    
    metric_magnitude = if_else(
      true_metric_noninvariance,
      delta_lambda,
      NA_real_
    ),
    
    scalar_magnitude = if_else(
      true_scalar_noninvariance,
      delta_nu,
      NA_real_
    )
  )


# ============================================================
# 4. Separate metric and scalar testing
#    PRIMARY = nominal design cells
# ============================================================

performance_long <- eval_main %>%
  select(
    method,
    popmodel,
    N,
    moderator,
    delta_lambda,
    delta_nu,
    metric_magnitude,
    scalar_magnitude,
    true_metric_noninvariance,
    true_scalar_noninvariance,
    metric_reject,
    scalar_reject
  ) %>%
  pivot_longer(
    cols = c(
      metric_reject,
      scalar_reject
    ),
    names_to = "level",
    values_to = "reject"
  ) %>%
  mutate(
    
    level = recode(
      level,
      metric_reject = "Metric",
      scalar_reject = "Scalar"
    ),
    
    true_noninvariance = case_when(
      level == "Metric" ~
        true_metric_noninvariance,
      level == "Scalar" ~
        true_scalar_noninvariance,
      TRUE ~ NA
    ),
    
    estimand = case_when(
      true_noninvariance == TRUE ~
        "Power",
      true_noninvariance == FALSE ~
        "Rejection u. invariance",
      TRUE ~ NA_character_
    )
  )

# ============================================================
# 5. Cell-specific power / rejection rate + MCSE
# ============================================================

perf_cell <- performance_long %>%
  group_by(
    method,
    popmodel,
    N,
    moderator,
    delta_lambda,
    delta_nu,
    level,
    true_noninvariance,
    estimand
  ) %>%
  summarise(
    
    n_total = n(),
    
    n_eval = sum(
      !is.na(reject)
    ),
    
    n_reject = sum(
      reject == TRUE,
      na.rm = TRUE
    ),
    
    rate = if_else(
      n_eval > 0,
      n_reject / n_eval,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    
    mcse = if_else(
      n_eval > 0,
      sqrt(
        rate * (1 - rate) /
          n_eval
      ),
      NA_real_
    ),
    
    mc_lower = pmax(
      0,
      rate - 1.96 * mcse
    ),
    
    mc_upper = pmin(
      1,
      rate + 1.96 * mcse
    ),
    
    missing_rate = if_else(
      n_total > 0,
      1 - n_eval / n_total,
      NA_real_
    )
  )

perf_cell

# ============================================================
# 6. Counterfactual sequential testing
# ============================================================

eval_main <- eval_main %>%
  mutate(
    
    # Would the scalar stage have been reached?
    seq_scalar_reached = case_when(
      is.na(metric_reject) ~ NA,
      metric_reject == TRUE ~ FALSE,
      metric_reject == FALSE ~ TRUE
    ),
    
    # Final decision under classical sequential testing
    seq_decision = case_when(
      
      metric_reject == TRUE ~
        "Metric noninvariance",
      
      metric_reject == FALSE &
        scalar_reject == TRUE ~
        "Scalar noninvariance",
      
      metric_reject == FALSE &
        scalar_reject == FALSE ~
        "Invariance retained",
      
      TRUE ~ NA_character_
    ),
    
    # True first level at which invariance fails
    true_seq_decision = case_when(
      
      true_metric_noninvariance == TRUE ~
        "Metric noninvariance",
      
      true_metric_noninvariance == FALSE &
        true_scalar_noninvariance == TRUE ~
        "Scalar noninvariance",
      
      true_metric_noninvariance == FALSE &
        true_scalar_noninvariance == FALSE ~
        "Invariance retained",
      
      TRUE ~ NA_character_
    ),
    
    truth_class = case_when(
      
      !true_metric_noninvariance &
        !true_scalar_noninvariance ~
        "No MNI",
      
      true_metric_noninvariance &
        !true_scalar_noninvariance ~
        "Metric only",
      
      !true_metric_noninvariance &
        true_scalar_noninvariance ~
        "Scalar only",
      
      true_metric_noninvariance &
        true_scalar_noninvariance ~
        "Metric + scalar",
      
      TRUE ~ NA_character_
    ),
    
    seq_correct =
      seq_decision == true_seq_decision
  )

# ============================================================
# 7. Cell-specific sequential performance
# ============================================================

perf_sequential_cell <- eval_main %>%
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
    
    # --------------------------------------------------------
    # Metric stage
    # --------------------------------------------------------
    
    n_metric_eval =
      sum(!is.na(metric_reject)),
    
    metric_rejection_rate =
      if_else(
        n_metric_eval > 0,
        sum(metric_reject == TRUE,
            na.rm = TRUE) /
          n_metric_eval,
        NA_real_
      ),
    
    # --------------------------------------------------------
    # Scalar stage
    # --------------------------------------------------------
    
    n_scalar_reached =
      sum(
        seq_scalar_reached == TRUE,
        na.rm = TRUE
      ),
    
    scalar_stage_reached =
      if_else(
        n_metric_eval > 0,
        n_scalar_reached /
          n_metric_eval,
        NA_real_
      ),
    
    # Rejection conditional on actually reaching scalar stage
    scalar_conditional_rejection =
      if_else(
        n_scalar_reached > 0,
        mean(
          scalar_reject[
            seq_scalar_reached == TRUE
          ],
          na.rm = TRUE
        ),
        NA_real_
      ),
    
    # Probability from the start that final decision is scalar MNI
    scalar_sequential_detection =
      mean(
        seq_decision ==
          "Scalar noninvariance",
        na.rm = TRUE
      ),
    
    # --------------------------------------------------------
    # Final sequential decisions
    # --------------------------------------------------------
    
    p_metric_decision =
      mean(
        seq_decision ==
          "Metric noninvariance",
        na.rm = TRUE
      ),
    
    p_scalar_decision =
      mean(
        seq_decision ==
          "Scalar noninvariance",
        na.rm = TRUE
      ),
    
    p_invariance_retained =
      mean(
        seq_decision ==
          "Invariance retained",
        na.rm = TRUE
      ),
    
    n_seq_eval =
      sum(!is.na(seq_correct)),
    
    sequential_accuracy =
      if_else(
        n_seq_eval > 0,
        mean(
          seq_correct,
          na.rm = TRUE
        ),
        NA_real_
      ),
    
    .groups = "drop"
  )
#------------------------------------
perf_sequential_cell <- perf_sequential_cell %>%
  mutate(
    
    mcse_metric =
      if_else(
        n_metric_eval > 0,
        sqrt(
          metric_rejection_rate *
            (1 - metric_rejection_rate) /
            n_metric_eval
        ),
        NA_real_
      ),
    
    mcse_scalar_stage =
      if_else(
        n_metric_eval > 0,
        sqrt(
          scalar_stage_reached *
            (1 - scalar_stage_reached) /
            n_metric_eval
        ),
        NA_real_
      ),
    
    mcse_scalar_conditional =
      if_else(
        n_scalar_reached > 0,
        sqrt(
          scalar_conditional_rejection *
            (1 - scalar_conditional_rejection) /
            n_scalar_reached
        ),
        NA_real_
      ),
    
    mcse_scalar_sequential =
      if_else(
        n_seq_eval > 0,
        sqrt(
          scalar_sequential_detection *
            (1 - scalar_sequential_detection) /
            n_seq_eval
        ),
        NA_real_
      ),
    
    mcse_sequential_accuracy =
      if_else(
        n_seq_eval > 0,
        sqrt(
          sequential_accuracy *
            (1 - sequential_accuracy) /
            n_seq_eval
        ),
        NA_real_
      )
  )

perf_sequential_cell

# ============================================================
# 8. Technical performance
# ============================================================

technical_cell <- eval_main %>%
  group_by(
    method,
    popmodel,
    N,
    moderator,
    delta_lambda,
    delta_nu
  ) %>%
  summarise(
    
    replications = n(),
    
    metric_available =
      sum(!is.na(metric_reject)),
    
    scalar_available =
      sum(!is.na(scalar_reject)),
    
    metric_missing =
      sum(is.na(metric_reject)),
    
    scalar_missing =
      sum(is.na(scalar_reject)),
    
    metric_missing_rate =
      mean(is.na(metric_reject)),
    
    scalar_missing_rate =
      mean(is.na(scalar_reject)),
    
    mean_runtime =
      mean(runtime_sec, na.rm = TRUE),
    
    median_runtime =
      median(runtime_sec, na.rm = TRUE),
    
    sd_runtime =
      sd(runtime_sec, na.rm = TRUE),
    
    .groups = "drop"
  )
#-----------------------
eval_dat
eval_main
eval_noise_sensitivity

performance_long
perf_cell

perf_sequential_cell

technical_cell
#-------------------------
method_failure_summary
# A tibble: 3 × 6
#method  attempted_runs failed_runs successful_runs failure_rate failure_percent
#<chr>            <int>       <int>           <int>        <dbl>           <dbl>
#  1 MNLFA            57600          11           57589     0.000191          0.0191
#2 MNLFAQ           57600          31           57569     0.000538          0.0538
#3 SEMTREE          57600           0           57600     0                 0     


overall_failure_summary <- eval_main %>%
  summarise(
    attempted_runs = n(),
    
    failed_runs = sum(
      is.na(metric_reject) |
        is.na(scalar_reject)
    ),
    
    failure_rate =
      failed_runs / attempted_runs,
    
    failure_percent =
      100 * failure_rate
  )

overall_failure_summary

#attempted_runs failed_runs failure_rate failure_percent
#1         172800          42 0.0002430556      0.02430556
#----------------------------

