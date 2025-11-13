# results_script.R
# =============================================================================
# Purpose
#   Present (as tibbles) and persist simulation results:
#     - raw per-replication results (rows from the driver)
#     - grouped summaries for MNLFA (δ metrics) and SEM-Tree (split metrics)
#
# Inputs
#   - `results`: tibble with columns produced by the driver, e.g.
#       model_type, N, reliability, lambda, moderator_1_type, rep,
#       delta_bias, delta_rmse, delta_coverage,
#       tree_split_on_M, tree_n_splits_M
#   - `out_dir`: directory where files are written (created if absent)
#   - `prefix`: filename prefix (e.g., "study_onefactor")
#
# Outputs (all written into out_dir with a timestamp suffix)
#   - CSV:   raw results, MNLFA summary, SEM-Tree summary
#   - RDS:   same three objects
#   - RData: list(raw=..., mnlfa=..., tree=...)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

# Optional: MCSE helpers if you want MCSEs in summaries
calculate_mcse_bias <- function(x) {
  x <- na.omit(x); if (!length(x)) return(NA_real_)
  sqrt(var(x) / length(x))
}
calculate_mcse_rmse <- function(x) {
  x <- na.omit(x); if (!length(x)) return(NA_real_)
  K <- length(x); m <- mean(x)
  sqrt(sum((x - m)^2) / (K * (K - 1)))
}

# --------------------------- summary builders --------------------------------

summarise_mnlfa <- function(results) {
  results %>%
    dplyr::group_by(model_type, N, moderator_1_type,
                    assumed_mnlfa_model, assumed_mnlfa_link) %>%
    summarise(
      mean_bias_delta  = mean(delta_bias,  na.rm = TRUE),
      mean_rmse_delta  = mean(delta_rmse,  na.rm = TRUE),
      mean_cov_delta   = mean(delta_coverage, na.rm = TRUE),
      mcse_bias        = calculate_mcse_bias(delta_bias),
      mcse_rmse        = calculate_mcse_rmse(delta_rmse),
      n_reps = dplyr::n(),
      n_reps           = sum(!is.na(delta_bias) | !is.na(delta_rmse) | !is.na(delta_coverage)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(model_type, N, moderator_1_type, assumed_mnlfa_model)
}

summarise_tree <- function(results) {
  results %>%
    group_by(model_type, N, moderator_1_type) %>%
    summarise(
      split_rate_M   = mean(tree_split_on_M,  na.rm = TRUE),   # proportion of reps with at least one split on M
      mean_splits_M  = mean(tree_n_splits_M, na.rm = TRUE),    # average number of M-splits
      n_reps         = sum(!is.na(tree_split_on_M)),
      .groups = "drop"
    ) %>%
    arrange(model_type, N, moderator_1_type)
}
# -----------------------------------------------------------------------------

# --------------------------- saving utility ----------------------------------
.ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

.timestamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

# ---- detection summary helper -----------------------------------------
# Power (TPR) and Type-I error (FPR) for SEM-Tree (M and M2) and MNLFA.
summarise_detection <- function(results) {
  results %>%
    dplyr::mutate(
      # "Truth": M moderation present if model_type is 1.1/1.2 and moderator_1_type ≠ "noise"
      truth_M_present = (model_type %in% c("1.1","1.2")) & (moderator_1_type != "noise"),
      # Tree detections (booleans already from run_one())
      tree_detect_M   = as.logical(tree_split_on_M),
      tree_detect_M2  = as.logical(tree_split_on_M2),     # always FPR target
      # MNLFA detection flag (TRUE/FALSE/NA), produced by run_one()
      mnlfa_detect    = as.logical(mnlfa_detect)
    ) %>%
    dplyr::group_by(model_type, N, moderator_1_type) %>%
    dplyr::summarise(
      # SEM-Tree: M
      TPR_tree_M = mean(tree_detect_M[ truth_M_present], na.rm = TRUE),
      FPR_tree_M = mean(tree_detect_M[!truth_M_present], na.rm = TRUE),
      # SEM-Tree: M2 (always noise ⇒ only FPR meaningful)
      FPR_tree_M2 = mean(tree_detect_M2, na.rm = TRUE),
      # MNLFA: δ significance
      TPR_mnlfa   = mean(mnlfa_detect[ truth_M_present], na.rm = TRUE),
      FPR_mnlfa   = mean(mnlfa_detect[!truth_M_present], na.rm = TRUE),
      n_reps      = dplyr::n(),
      .groups = "drop"
    )
}


# Main convenience wrapper: present and persist
present_and_save_results <- function(results, out_dir = "results", prefix = "study_onefactor") {
  stopifnot(is.data.frame(results))
  .ensure_dir(out_dir)
  ts <- .timestamp()
  
  # Build summaries
  summary_mnlfa  <- summarise_mnlfa(results)
  summary_tree   <- summarise_tree(results)
  summary_detect <- summarise_detection(results)  
  
  # Present (console)
  cat("\n=== RAW results (head) ===\n")
  print(tibble::as_tibble(results) %>% head(10))
  cat("\n=== MNLFA summary ===\n")
  print(summary_mnlfa)
  cat("\n=== SEM-Tree summary ===\n")
  print(summary_tree)
  cat("\n=== Detection summary ===\n")
  print(summary_detect) 
  
  # Filenames
  f_raw_csv    <- file.path(out_dir, sprintf("%s_raw_%s.csv",      prefix, ts))
  f_mnlfa_csv  <- file.path(out_dir, sprintf("%s_mnlfa_%s.csv",    prefix, ts))
  f_tree_csv   <- file.path(out_dir, sprintf("%s_tree_%s.csv",     prefix, ts))
  f_detect_csv <- file.path(out_dir, sprintf("%s_detect_%s.csv",   prefix, ts))  
  
  f_raw_rds    <- file.path(out_dir, sprintf("%s_raw_%s.rds",      prefix, ts))
  f_mnlfa_rds  <- file.path(out_dir, sprintf("%s_mnlfa_%s.rds",    prefix, ts))
  f_tree_rds   <- file.path(out_dir, sprintf("%s_tree_%s.rds",     prefix, ts))
  f_detect_rds <- file.path(out_dir, sprintf("%s_detect_%s.rds",   prefix, ts))  
  
  f_bundle_RData <- file.path(out_dir, sprintf("%s_bundle_%s.RData", prefix, ts))
  
  # Persist as CSV
  readr::write_csv(results,        f_raw_csv)
  readr::write_csv(summary_mnlfa,  f_mnlfa_csv)
  readr::write_csv(summary_tree,   f_tree_csv)
  readr::write_csv(summary_detect, f_detect_csv)   
  
  # Persist as RDS
  saveRDS(results,        f_raw_rds)
  saveRDS(summary_mnlfa,  f_mnlfa_rds)
  saveRDS(summary_tree,   f_tree_rds)
  saveRDS(summary_detect, f_detect_rds)            
  
  # Bundle as .RData
  raw <- results; mnlfa <- summary_mnlfa; tree <- summary_tree 
  detect <- summary_detect
  save(raw, mnlfa, tree, detect, file = f_bundle_RData)                                   
  
  invisible(list(
    csv   = list(raw = f_raw_csv, mnlfa = f_mnlfa_csv, tree = f_tree_csv, detect = f_detect_csv),  
    rds   = list(raw = f_raw_rds, mnlfa = f_mnlfa_rds, tree = f_tree_rds, detect = f_detect_rds),  
    rdata = f_bundle_RData
  ))
}

