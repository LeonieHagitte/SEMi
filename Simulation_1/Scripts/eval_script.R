# ---- eval_script.R ----
# Purpose: aggregate and summarise simulation outcomes

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

# source functions
source(here("Simulation_1", "Functions", "eval_functions.R"))

results_dir <- here("Simulation_1", "Results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# --------- 1) Load aggregate per-condition results ----------
results_list <- readRDS(file.path(results_dir, "sim_raw_results.rds"))

# Aggregate and compute ΔAIC/ΔBIC 
agg <- aggregate_results(results_list)
summary_det <- summarise_detection(agg)

# Persist detection summaries
readr::write_csv(agg,         file.path(results_dir, "evaluation_aggregated.csv"))
readr::write_csv(summary_det, file.path(results_dir, "evaluation_summary_detection.csv"))
saveRDS(agg,         file.path(results_dir, "evaluation_aggregated.rds"))
saveRDS(summary_det, file.path(results_dir, "evaluation_summary_detection.rds"))

# --------- 2) Full metrics only if per-replications are available ----------
reps_path <- file.path(results_dir, "sim_rep_records.rds")
have_reps <- file.exists(reps_path)

if (!have_reps) {
  warning(
    "Per-replication file 'sim_rep_records.rds' not found. ",
    "Full parameter metrics and runtime summaries will be skipped. ",
    "To enable them, write per-rep records during simulation (see notes)."
  )
  message("Saved detection-only summaries: ",
          "evaluation_aggregated.(csv|rds), evaluation_summary_detection.(csv|rds)")
} else {
  
  reps_df <- readRDS(reps_path)
  
  # Coerce possible list-cols for item sets to compact strings
  coerce_items <- function(x) {
    if (is.list(x)) vapply(x, function(el) paste0(el, collapse = ","), character(1)) else as.character(x)
  }
  if ("mod_items_load" %in% names(reps_df)) reps_df$mod_items_load <- coerce_items(reps_df$mod_items_load)
  if ("mod_items_int"  %in% names(reps_df)) reps_df$mod_items_int  <- coerce_items(reps_df$mod_items_int)
  if ("mod_items_res"  %in% names(reps_df)) reps_df$mod_items_res  <- coerce_items(reps_df$mod_items_res)
  
  # Condition keys aligned with simulation design
  condition_keys <- c("n", "moderator_type",
                      "mod_items_load", "mod_items_int", "mod_items_res",
                      "slope_load", "slope_int", "slope_logres")
  
  # ---- 2a) Full single-target metrics by condition x method (includes semtree) ----
  # evaluate_by_condition() should already group by method and compute:
  # convergence rate, rejection rate, bias/rel. bias, MSE/RMSE, emp SE, model SE,
  # SE rel. error, 95% coverage, CI width, MCSEs, etc.
  eval_full <- evaluate_by_condition(reps_df, condition_keys = condition_keys)
  
  # ---- 2b) Runtimes by method x condition (includes semtree) ----
  runtime_cols <- c("method", "runtime_seconds")
  if (all(runtime_cols %in% names(reps_df))) {
    runtimes_summary <- summarise_runtimes(reps_df[, c(condition_keys, runtime_cols)], condition_keys)
  } else {
    runtimes_summary <- tibble()
    warning("Runtime columns not found; skipping runtime summaries.")
  }
  
  # ---- 2c) Optional: SEM-tree specific detection summaries (new) ----
  # These use columns you recorded only for method == 'semtree':
  #   any_split (logical), split_on_moderator (logical), n_splits (integer)
  if ("method" %in% names(reps_df) &&
      any(reps_df$method == "semtree")) {
    
    semtree_df <- reps_df %>% filter(method == "semtree")
    
    # helper that returns NA if no rows (rather than 0)
    mean_or_na <- function(x) if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE)
    
    semtree_summary <- semtree_df %>%
      group_by(across(all_of(condition_keys))) %>%
      summarise(
        det_rate_tree_any    = mean_or_na(any_split),
        det_rate_tree_on_mod = mean_or_na(split_on_moderator),
        mean_tree_splits     = mean_or_na(n_splits),
        .groups = "drop"
      )
    
    readr::write_csv(semtree_summary, file.path(results_dir, "evaluation_semtree_detection.csv"))
    saveRDS(semtree_summary,          file.path(results_dir, "evaluation_semtree_detection.rds"))
  }
  
  # Save main outputs
  readr::write_csv(eval_full,        file.path(results_dir, "evaluation_metrics_full.csv"))
  readr::write_csv(runtimes_summary, file.path(results_dir, "evaluation_runtimes.csv"))
  saveRDS(eval_full,        file.path(results_dir, "evaluation_metrics_full.rds"))
  saveRDS(runtimes_summary, file.path(results_dir, "evaluation_runtimes.rds"))
  
  message("Saved: evaluation_metrics_full.(csv|rds), evaluation_runtimes.(csv|rds), ",
          "plus detection summaries: evaluation_aggregated.(csv|rds), evaluation_summary_detection.(csv|rds)",
          if (exists("semtree_summary")) ", evaluation_semtree_detection.(csv|rds)" else "")
}
