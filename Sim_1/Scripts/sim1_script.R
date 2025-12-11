# minimal driver 
# ----------------------------------------------------------------------

library(here)
library(dplyr)
library(purrr)
library(OpenMx)

# dependencies
source(here("Sim_1", "Functions", "dataprep_functions2.R"))
source(here("Sim_1", "Functions", "analysis.R"))
source(here("Sim_1", "Functions", "dataprep_wrapper.R"))
source(here("Sim_1", "Functions", "eval_functions.R"))
source(here("Sim_1", "Scripts", "results_script.R"))

results_dir <- here("Sim_1", "Results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# --------------------- DESIGN (edit here) ------------------------------------
MOD_TYPES <- c("linear","sigmoid","quadratic","noise")

DESIGN <- tidyr::expand_grid(
  model_type        = c("NULL","1.1","1.2"),   # DGP pattern of true deltas
  N                 = c(250, 500),
  reliability       = 0.80,
  lambda            = 0.70,
  delta_lambda_full = 0.20,
  delta_theta_full  = 0.20,
  delta_lambda_12   = 0.20,
  moderator_1_type  = MOD_TYPES,               
  rep               = 1#:10
)

set.seed(42)

# ---------------
# helper to print item sets as compact strings
.items_str <- function(idx) if (length(idx)) paste0("y", idx, collapse = ",") else "none"


# --------------------- RUN ONE REPLICATION -----------------------------------
run_one <- function(row) {
  # 1) Simulate per-case data
  sim <- simulate_moderated_onefactor(
    model_type        = row$model_type,
    N                 = row$N,
    reliability       = row$reliability,
    lambda            = row$lambda,
    moderator_1_type  = row$moderator_1_type,
    delta_lambda_full = row$delta_lambda_full,
    delta_theta_full  = row$delta_theta_full,
    delta_lambda_12   = row$delta_lambda_12
  )
  dat    <- sim$data
  truths <- sim$true_deltas
  
  # 2) Analysis syntaxes 
  mnlfa_variants <- list(
    MNLFA_linear_full    = build_mnlfa_linear_full(),
    MNLFA_linear_partial = build_mnlfa_linear_partial(),
    MNLFA_none           = build_cfa_baseline()   # include unmoderated MNLFA
  )
  cfa_base <- build_cfa_baseline()
  delta_names <- paste0("delta_lambda_y", 1:3)
  
  # Iterate over *all* variants
  rows <- lapply(names(mnlfa_variants), function(vname) {
    mnlfa_model <- mnlfa_variants[[vname]]
    
    # --- CI targets present in THIS syntax ---
    assumed_M <- assumptions_from_syntax(mnlfa_model, method = "MNLFA")
    lambda_targets <- if (length(assumed_M$load_items))  paste0("delta_lambda_y", assumed_M$load_items) else character(0)
    theta_targets  <- if (length(assumed_M$theta_items)) paste0("delta_theta_y",  assumed_M$theta_items) else character(0)
    delta_names    <- c(lambda_targets, theta_targets)
    need_intervals <- length(delta_names) > 0
# --------------------------------------    
    # 3) Fit methods (timed)
    mx_t <- system.time({
      fit_mnlfa <- run_analysis(
        dat, mnlfa_model,
        method       = "MNLFA",
        mx_intervals = need_intervals,   # only ask for CIs if targets exist
        delta_names  = delta_names       # may be character(0)
      )
    })[["elapsed"]]
    
    tr_t <- system.time({
      fit_tree <- run_analysis(dat, cfa_base, method = "SEMTREE")
    })[["elapsed"]]
    
    # --- δ metrics with error handling ---
    if (inherits(fit_mnlfa, c("MxModel","MxRAMModel"))) {
      dm <- delta_metrics_any(fit_mnlfa, truths, label_prefix = "delta_")
    } else {
      if (inherits(fit_mnlfa, "simpleError")) print(fit_mnlfa)  # optional logging
      dm <- list(bias = NA_real_, rmse_abs = NA_real_, coverage = NA_real_)
    }
    delta_bias     <- if (is.null(dm$bias))      NA_real_ else dm$bias
    delta_rmse     <- if (is.null(dm$rmse_abs))  NA_real_ else dm$rmse_abs
    delta_coverage <- if (is.null(dm$coverage))  NA_real_ else dm$coverage
    
    # --- MNLFA detection (single robust block, removed duplication) ---
    if (inherits(fit_mnlfa, c("MxModel","MxRAMModel"))) {
      det_m <- tryCatch(
        detect_moderation_mxsem(fit_mnlfa, label_prefix = "delta_", alpha = 0.05),
        error = function(e) list(any = NA, n_sig = NA_integer_, per_param = NULL)
      )
    } else {
      det_m <- list(any = NA, n_sig = NA_integer_, per_param = NULL)
    }
    mnlfa_detect <- isTRUE(det_m$any)
    mnlfa_n_sig  <- if (is.null(det_m$n_sig)) NA_integer_ else as.integer(det_m$n_sig)
    
    # --- SEM-Tree split detection ---
    predictors <- getPredictorsFromTree(fit_tree)
    k_M  <- if (is.null(predictors)) 0L else sum(predictors == "M")
    k_M2 <- if (is.null(predictors)) 0L else sum(predictors == "M2")
    tree_split_on_M  <- (k_M  > 0)
    tree_split_on_M2 <- (k_M2 > 0)
    
    # --- truth & assumptions ---
    truth     <- truth_from_sim(sim, link_label = row$moderator_1_type)
    assumed_M <- assumptions_from_syntax(mnlfa_model, method = "MNLFA")
    assumed_T <- assumptions_from_syntax(cfa_base,    method = "SEMTREE")
    
    # 4) one tidy row per variant
    tibble::tibble(
      model_type        = row$model_type,
      N                 = row$N,
      reliability       = row$reliability,
      lambda            = row$lambda,
      moderator_1_type  = row$moderator_1_type,
      variant_label     = vname,  #identify which MNLFA model
      
      mnlfa_detect      = mnlfa_detect,
      mnlfa_n_sig       = mnlfa_n_sig,
      tree_split_on_M   = tree_split_on_M,
      tree_n_splits_M   = as.integer(k_M),
      tree_split_on_M2  = tree_split_on_M2,
      tree_n_splits_M2  = as.integer(k_M2),
      
      delta_lambda_full = row$delta_lambda_full,
      delta_theta_full  = row$delta_theta_full,
      delta_lambda_12   = row$delta_lambda_12,
      rep               = row$rep,
      
      delta_bias        = delta_bias,
      delta_rmse        = delta_rmse,
      delta_coverage    = delta_coverage,
      
      runtime_mxsem     = as.numeric(mx_t),
      runtime_semtree   = as.numeric(tr_t),
      
      # --- TRUE (DGP)
      true_link           = truth$true_link,
      true_load_items     = .items_str(truth$true_load_items),
      true_theta_items    = .items_str(truth$true_theta_items),
      true_ref_moderated  = truth$true_ref_moderated,
      
      # --- ASSUMED (MNLFA)
      assumed_mnlfa_model         = assumed_M$model_label,
      assumed_mnlfa_link          = assumed_M$link,
      assumed_mnlfa_load_items    = .items_str(assumed_M$load_items),
      assumed_mnlfa_theta_items   = .items_str(assumed_M$theta_items),
      assumed_mnlfa_ref_moderated = assumed_M$ref_moderated,
      
      # --- ASSUMED (SEMTREE)
      assumed_semtree_model      = assumed_T$model_label,
      assumed_semtree_parametric = assumed_T$parametric,
      assumed_semtree_link       = assumed_T$link
    )
  })
  
  dplyr::bind_rows(rows)  # return all rows together
}


# --------------------- EXECUTE -----------------------------------------------
results <- DESIGN %>%
  split(.$rep) %>%
  purrr::map_dfr(~ purrr::pmap_dfr(.x, function(...) run_one(list(...))))


# --------------------- SUMMARIES ---------------------------------------------
summary_mnlfa <- results %>%
  dplyr::group_by(model_type, N) %>%
  dplyr::summarise(
    mean_cov_delta  = mean(delta_coverage, na.rm = TRUE),
    mean_bias_delta = mean(delta_bias,     na.rm = TRUE),
    mean_rmse_delta = mean(delta_rmse,     na.rm = TRUE),
    n               = dplyr::n(),
    .groups = "drop"
  )


summary_tree <- results %>%
  group_by(model_type, N) %>%
  summarise(
    split_rate_M    = mean(tree_split_on_M, na.rm = TRUE),
    mean_splits_M   = mean(tree_n_splits_M, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_mnlfa)
print(summary_tree)

# -----------------------------------------------------------------------------
evaluation_runtimes <- results %>%
  dplyr::select(N, moderator_1_type, variant_label, runtime_mxsem, runtime_semtree) %>%
  tidyr::pivot_longer(
    c(runtime_mxsem, runtime_semtree),
    names_to  = "method",
    values_to = "runtime_seconds"
  ) %>%
  dplyr::mutate(
    method = dplyr::recode(method,
                           runtime_mxsem   = ifelse(variant_label == "MNLFA_none", "mxsem_null", "mxsem"),
                           runtime_semtree = "semtree"
    )
  )

# ---- Persist results (bottom of sim1_script.R) ------------------------------
# 'results' must be a tibble/data.frame with columns produced by run_one():
#   model_type, N, reliability, lambda, moderator_1_type, rep,
#   delta_bias, delta_rmse, delta_coverage, tree_split_on_M, tree_n_splits_M, ...

stopifnot(exists("results"), is.data.frame(results))

# quick peek
cat("\nFinished. First rows:\n"); print(utils::head(results, 10))

# Save raw + summaries to disk (CSV, RDS, and a bundled .RData)
out_files <- present_and_save_results(
  results,
  out_dir = "results",              # folder is created if missing
  prefix  = "sim1_onefactor"        # tag this run (appears in filenames)
)

# show file paths
str(out_files)

