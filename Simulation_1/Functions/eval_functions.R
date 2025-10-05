# ---- eval_functions.R ----
# Purpose: post-simulation aggregation and full Monte Carlo evaluation helpers

suppressPackageStartupMessages({
  library(tidyverse)
})

# -------------------- Basic aggregations already in your pipeline --------------------

# Combine per-condition summary tibbles (rows) and compute ΔAIC/ΔBIC.
aggregate_results <- function(results_list) {
  df <- dplyr::bind_rows(results_list)
  df %>%
    mutate(
      dAIC = AIC_alt - AIC_null,
      dBIC = BIC_alt - BIC_null
    )
}

# Summaries of detection by family and overall (retained for continuity).
summarise_detection <- function(agg_df) {
  agg_df %>%
    group_by(n, mod_items_load, mod_items_int, mod_items_res,
             slope_load, slope_int, slope_logres, moderator_type, .drop = FALSE) %>%
    summarise(
      mean_det_load = mean(det_rate_load),
      mean_det_int  = mean(det_rate_int),
      mean_det_res  = mean(det_rate_res),
      mean_det_any  = mean(det_rate_any),
      mean_dAIC     = mean(dAIC),
      mean_dBIC     = mean(dBIC),
      .groups = "drop"
    )
}

# -------------------- Full single-parameter evaluation metrics --------------------

# Utility: MCSE for a percentage (0–100) computed from Bernoulli(p) with n trials.
.mcse_pct <- function(pct, n) {
  # pct in [0,100]; return MCSE in percentage points.
  if (is.na(pct) || is.na(n) || n <= 0) return(NA_real_)
  p <- pct / 100
  sqrt(p * (1 - p) / n) * 100
}

# Wrapper: evaluate parameter by condition. Input is a data frame
# with columns that at least include: rep, estimate, se, p, converged, true_value..

evaluate_by_condition <- function(reps_df, condition_keys) {
  df <- reps_df %>%
    filter(is.finite(estimate), is.finite(se)) %>%
    mutate(sig = (p < 0.05))
  
  group_keys <- c(condition_keys, "method", "target_param_name", "true_value")
  
  df %>%
    group_by(across(all_of(group_keys))) %>%
    summarise(
      total_runs = n(),
      converged_runs = sum(isTRUE(converged) | method == "semtree", na.rm = TRUE),
      reject_rate = mean(sig, na.rm = TRUE) * 100,
      mcse_reject_rate = sqrt(reject_rate * (100 - reject_rate) / pmax(1, converged_runs)),
      
      mean_estimate = mean(estimate, na.rm = TRUE),
      bias = mean(estimate, na.rm = TRUE) - unique(true_value)[1],
      abs_bias = abs(bias),
      emp_se = sd(estimate, na.rm = TRUE),
      mod_se = mean(se, na.rm = TRUE),
      rel_error_se = ((mod_se / emp_se) - 1) * 100,
      
      mse = mean((estimate - unique(true_value)[1])^2, na.rm = TRUE),
      rmse = sqrt(mse),
      
      mcse_bias = sd(estimate, na.rm = TRUE) / sqrt(pmax(1, converged_runs)),
      mcse_emp_se = emp_se / sqrt(2 * pmax(1, converged_runs - 1)),
      
      coverage = mean(
        (estimate - qnorm(.975) * se) <= unique(true_value)[1] &
          (estimate + qnorm(.975) * se) >= unique(true_value)[1],
        na.rm = TRUE
      ) * 100,
      mean_ci_width = mean((estimate + qnorm(.975) * se) - (estimate - qnorm(.975) * se), na.rm = TRUE),
      
      .groups = "drop"
    )
}



# -------------------- Runtime summaries --------------------

# Summarise runtimes per method and condition.
# Expects columns: method, runtime_seconds (numeric), plus condition keys.
summarise_runtimes <- function(df, condition_keys) {
  df %>%
    dplyr::group_by(dplyr::across(all_of(c(condition_keys, "method")))) %>%
    dplyr::summarise(
      n_reps   = dplyr::n(),
      n_valid  = sum(is.finite(runtime_seconds)),
      Mean     = ifelse(n_valid > 0, mean(runtime_seconds, na.rm = TRUE), NA_real_),
      Min      = ifelse(n_valid > 0, min(runtime_seconds,  na.rm = TRUE), NA_real_),
      Max      = ifelse(n_valid > 0, max(runtime_seconds,  na.rm = TRUE), NA_real_),
      Total    = ifelse(n_valid > 0, sum(runtime_seconds,  na.rm = TRUE), NA_real_),
      .groups  = "drop"
    )
}

