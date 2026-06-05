library(furrr)
library(dplyr)
library(tibble)
library(future)


source("dataprep_rt.R")
source("analysis_rt.R")

# -----------------------------------
set.seed(42)
# ---------------------------------
n_rep <- 100 # needs to be larger, or as large as SLURM_ARRAY

# ----------------------------------
MOD_TYPES <- c("linear","sigmoid","quadratic","noise")

DESIGN <- tidyr::expand_grid(
  popmodel     = c("0","1.1", "1.11", "1.12","1.2","1.21","1.22","1.3","1.32"),
  N            = c(300, 500, 700, 1000),
  reliability  = c(0.60, 0.70, 0.80, 0.95),
  lambda       = 0.70,
  intercepts   = 1,
  # latentmean  = 0,
  # delta_eta   = c(-1, -0.5, 0.5, 1),
  delta_lambda = c(-0.3, -0.2, 0.2, 0.3),
  delta_nu     = c(-1, -0.5, 0.5, 1),
  moderator    = MOD_TYPES,
  analysis_form = c("linear", "quadratic"),
  rep_id = 1:n_rep
)%>%
  dplyr::arrange(
    popmodel, N, reliability, lambda, intercepts,
    delta_lambda, delta_nu, moderator,analysis_form, rep_id
  ) %>%
  dplyr::mutate(
    job_id = dplyr::row_number()
  ) 

DESIGN <- DESIGN %>%
  dplyr::mutate(
    # randomly draw one seed per each row
    seed = round(runif(nrow(DESIGN),0,.Machine$integer.max))
    
  )
# ---------- randomly permute lines for even distribution of
#     run times across jobs ----------------
DESIGN <- DESIGN %>%
  slice_sample(prop = 1) # this is a random permutation

# ---------- Split on rep_id, for similarly long runtimes ------
chunk_id <- NULL
n_chunks <- NULL

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  chunk_id <- as.integer(args[1])
  n_chunks <- as.integer(args[2])
  
  rep_ids <- sort(unique(DESIGN$rep_id))
  
  if (n_chunks > length(rep_ids)) {
    stop("n_chunks cannot be larger than number of replications.")
  }
  
  rep_chunks <- split(
    rep_ids,
    cut(seq_along(rep_ids), n_chunks, labels = FALSE)
  )
  
  reps_this_chunk <- rep_chunks[[chunk_id]]
  
  DESIGN <- DESIGN %>%
    dplyr::filter(rep_id %in% reps_this_chunk)
}
# --------------------------------------------------------------
mnlfa_moderation_estimate_names <- function(p = 4) {
  base_names <- c(
    "mnlfa_est_dnu_am1",
    "mnlfa_est_dnu_am2",
    "mnlfa_est_dnu_am12",
    "mnlfa_est_dlambda_am1",
    "mnlfa_est_dlambda_am2",
    "mnlfa_est_dlambda_am12"
  )
  as.vector(outer(base_names, paste0("x", seq_len(p)), paste, sep = "_"))
}

empty_mnlfa_moderation_estimates <- function(p = 4) {
  cols <- mnlfa_moderation_estimate_names(p = p)
  tibble::as_tibble(as.list(stats::setNames(rep(NA_real_, length(cols)), cols)))
}

flatten_mnlfa_moderation_estimates <- function(mnlfa_result, p = 4) {
  out <- empty_mnlfa_moderation_estimates(p = p)
  
  if (inherits(mnlfa_result, "error") ||
      is.null(mnlfa_result$configural_moderation_estimates)) {
    return(out)
  }
  
  est <- mnlfa_result$configural_moderation_estimates
  
  est_wide <- est %>%
    tidyr::pivot_wider(
      names_from = item,
      values_from = c(
        mnlfa_est_dnu_am1,
        mnlfa_est_dnu_am2,
        mnlfa_est_dnu_am12,
        mnlfa_est_dlambda_am1,
        mnlfa_est_dlambda_am2,
        mnlfa_est_dlambda_am12
      )
    )
  
  for (nm in intersect(names(out), names(est_wide))) {
    out[[nm]] <- est_wide[[nm]]
  }
  
  out
}


##############################################################################
run_one <- function(row) { #run_one <- function(seed, N, popmodel, moderator) 

  set.seed(row$seed)

  popmodel_use <- row$popmodel
  
  # ---------------------------
  params <- gen_paramsC(
    popmodel      = popmodel_use,
    lambda        = row$lambda,
    nu            = row$intercepts,
    reliability   = row$reliability,
    moderator     = row$moderator,
    delta_lambda  = row$delta_lambda,
    delta_nu      = row$delta_nu
  )
  
  sim <- gen_dataC(
    N = row$N,
    params = params
  )
  
  df <- sim$data
  
  df <- add_analysis_form(
    data = df,
    analysis_form = row$analysis_form,
    k = params$k
  )
  # ensure required columns exist for analyses
  if (!"m0" %in% names(df)) df$m0 <- 0
  
  # ---------------------------
  res <- run_analysis(
    data = df,
    methods = c("MNLFA", "SEMTREE"),
    nfactors = 1,
    alpha = 0.05,
    predictors = c("am1", "am2", "m0")
  )
  # ---------------------------
  mnlfa_error_msg <- NA_character_
  semtree_error_msg <- NA_character_
  
  if (inherits(res$mnlfa, "error")) {
    mnlfa_error_msg <- conditionMessage(res$mnlfa)
  }
  
  if (inherits(res$semtree, "error")) {
    semtree_error_msg <- conditionMessage(res$semtree)
  }
  
  
  mnlfa_mod_est <- flatten_mnlfa_moderation_estimates(
    mnlfa_result = res$mnlfa,
    p = 4
  )
  
  # ---------------------------
  
  mnlfa_metric_lrt_chisq <- NA_real_
  mnlfa_metric_lrt_df <- NA_real_
  mnlfa_metric_lrt_p <- NA_real_
  mnlfa_metric_lrt_reject <- NA
  
  mnlfa_scalar_lrt_chisq <- NA_real_
  mnlfa_scalar_lrt_df <- NA_real_
  mnlfa_scalar_lrt_p <- NA_real_
  mnlfa_scalar_lrt_reject <- NA
  
  mnlfa_omnibus_lrt_chisq <- NA_real_
  mnlfa_omnibus_lrt_df <- NA_real_
  mnlfa_omnibus_lrt_p <- NA_real_
  mnlfa_omnibus_lrt_reject <- NA
  
  if (!inherits(res$mnlfa, "error")) {
    if (!is.null(res$mnlfa$metric_lrt)) {
      mnlfa_metric_lrt_chisq <- res$mnlfa$metric_lrt$chisq_diff
      mnlfa_metric_lrt_df <- res$mnlfa$metric_lrt$df_diff
      mnlfa_metric_lrt_p <- res$mnlfa$metric_lrt$p_value
      mnlfa_metric_lrt_reject <- res$mnlfa$metric_lrt$reject_h0
    }
    
    if (!is.null(res$mnlfa$scalar_lrt)) {
      mnlfa_scalar_lrt_chisq <- res$mnlfa$scalar_lrt$chisq_diff
      mnlfa_scalar_lrt_df <- res$mnlfa$scalar_lrt$df_diff
      mnlfa_scalar_lrt_p <- res$mnlfa$scalar_lrt$p_value
      mnlfa_scalar_lrt_reject <- res$mnlfa$scalar_lrt$reject_h0
    }
    
    if (!is.null(res$mnlfa$omnibus_lrt)) {
      mnlfa_omnibus_lrt_chisq <- res$mnlfa$omnibus_lrt$chisq_diff
      mnlfa_omnibus_lrt_df <- res$mnlfa$omnibus_lrt$df_diff
      mnlfa_omnibus_lrt_p <- res$mnlfa$omnibus_lrt$p_value
      mnlfa_omnibus_lrt_reject <- res$mnlfa$omnibus_lrt$reject_h0
    }
  }
  
  mnlfa_model <- NA_character_
  if (!inherits(res$mnlfa, "error")) {
    if (!is.null(res$mnlfa$fitScalar)) {
      mnlfa_model <- "scalar"
    } else if (!is.null(res$mnlfa$fitMetric)) {
      mnlfa_model <- "metric"
    } else if (!is.null(res$mnlfa$fitConfig)) {
      mnlfa_model <- "configural"
    }
  }
  
  mnlfa_det <- NA
  
  if (identical(mnlfa_metric_lrt_reject, TRUE) ||
      identical(mnlfa_scalar_lrt_reject, TRUE)) {
    mnlfa_det <- TRUE
  } else if (identical(mnlfa_metric_lrt_reject, FALSE) &&
             (identical(mnlfa_scalar_lrt_reject, FALSE) ||
              is.na(mnlfa_scalar_lrt_reject))) {
    mnlfa_det <- FALSE
  }
  
  mnlfa_final_decision <- NA_character_
  
  if (identical(mnlfa_metric_lrt_reject, TRUE)) {
    mnlfa_final_decision <- "noninvariance_at_metric_lrt"
  } else if (identical(mnlfa_metric_lrt_reject, FALSE) &&
             identical(mnlfa_scalar_lrt_reject, TRUE)) {
    mnlfa_final_decision <- "noninvariance_at_scalar_lrt"
  } else if (identical(mnlfa_metric_lrt_reject, FALSE) &&
             (identical(mnlfa_scalar_lrt_reject, FALSE) ||
              is.na(mnlfa_scalar_lrt_reject))) {
    mnlfa_final_decision <- "scalar_invariance_retained_lrt"
  }
  # truth indicators
  # ---------------------------
  has_metric <- row$popmodel %in% c("1.1", "1.12", "1.2", "1.22", "1.3", "1.32")
  has_scalar <- row$popmodel %in% c("1.11", "1.12", "1.21", "1.22", "1.32")
  
  true_structured_moderator <- row$moderator != "noise"
  
  true_metric_noninvariance <- has_metric &&
    row$delta_lambda != 0 &&
    true_structured_moderator
  
  true_scalar_noninvariance <- has_scalar &&
    row$delta_nu != 0 &&
    true_structured_moderator
  
  true_any_noninvariance <- true_metric_noninvariance ||
    true_scalar_noninvariance
  # ---------------------------
  
  error_msg <- NA_character_
  
  tree_metric_split <- NA
  tree_scalar_split <- NA
  
  tree_metric_split_on_am1 <- NA
  tree_metric_split_on_am2 <- NA
  tree_metric_n_splits_am1 <- NA_integer_
  tree_metric_n_splits_am2 <- NA_integer_
  
  tree_scalar_split_on_am1 <- NA
  tree_scalar_split_on_am2 <- NA
  tree_scalar_n_splits_am1 <- NA_integer_
  tree_scalar_n_splits_am2 <- NA_integer_
  
  tree_metric_split_on_m0 <- NA
  tree_metric_n_splits_m0 <- NA_integer_
  
  tree_scalar_split_on_m0 <- NA
  tree_scalar_n_splits_m0 <- NA_integer_
  
  tree_metric_p <- NA_real_
  tree_metric_p_uncorrected <- NA_real_
  tree_metric_reject <- NA
  
  tree_scalar_p <- NA_real_
  tree_scalar_p_uncorrected <- NA_real_
  tree_scalar_reject <- NA
  
  tree_metric_correct_split <- NA
  tree_scalar_correct_split <- NA
  
  if (!inherits(res$semtree, "error")) {
    
    tree_metric_split <- res$semtree$metric_split
    tree_scalar_split <- res$semtree$scalar_split
    
    if (!is.null(res$semtree$metric_test)) {
      tree_metric_p <- res$semtree$metric_test$p_value
      tree_metric_p_uncorrected <- res$semtree$metric_test$p_uncorrected
      tree_metric_reject <- res$semtree$metric_test$reject_h0
    }
    
    if (!is.null(res$semtree$scalar_test)) {
      tree_scalar_p <- res$semtree$scalar_test$p_value
      tree_scalar_p_uncorrected <- res$semtree$scalar_test$p_uncorrected
      tree_scalar_reject <- res$semtree$scalar_test$reject_h0
    }
    
    metric_info <- semtree_detects_moderation(
      res$semtree$metric_tree,
      moderators = c("am1", "am2", "m0")
    )
    
    scalar_info <- semtree_detects_moderation(
      res$semtree$scalar_tree,
      moderators = c("am1", "am2", "m0")
    )
    
    tree_metric_split_on_m0 <- metric_info$tree_split_on_m0
    tree_metric_n_splits_m0 <- metric_info$tree_n_splits_m0
    
    tree_scalar_split_on_m0 <- scalar_info$tree_split_on_m0
    tree_scalar_n_splits_m0 <- scalar_info$tree_n_splits_m0
    
    tree_metric_split_on_am1 <- metric_info$tree_split_on_am1
    tree_metric_split_on_am2 <- metric_info$tree_split_on_am2
    tree_metric_n_splits_am1 <- metric_info$tree_n_splits_am1
    tree_metric_n_splits_am2 <- metric_info$tree_n_splits_am2
    
    tree_scalar_split_on_am1 <- scalar_info$tree_split_on_am1
    tree_scalar_split_on_am2 <- scalar_info$tree_split_on_am2
    tree_scalar_n_splits_am1 <- scalar_info$tree_n_splits_am1
    tree_scalar_n_splits_am2 <- scalar_info$tree_n_splits_am2
    
    true_metric_moderators <- character(0)
    true_scalar_moderators <- character(0)
    
    if (true_structured_moderator) {
      if (row$popmodel %in% c("1.1", "1.12", "1.2", "1.22")) {
        true_metric_moderators <- c(true_metric_moderators, "am1")
      }
      
      if (row$popmodel == "1.3") {
        true_metric_moderators <- c(true_metric_moderators, "am1", "am2")
      }
      
      if (row$popmodel == "1.32") {
        true_metric_moderators <- c(true_metric_moderators, "am1")
        true_scalar_moderators <- c(true_scalar_moderators, "am2")
      }
      
      if (row$popmodel %in% c("1.11", "1.12", "1.21", "1.22")) {
        true_scalar_moderators <- c(true_scalar_moderators, "am1")
      }
    }
    
    tree_metric_selected <- c(
      if (isTRUE(tree_metric_split_on_am1)) "am1",
      if (isTRUE(tree_metric_split_on_am2)) "am2",
      if (isTRUE(tree_metric_split_on_m0))  "m0"
    )
    
    tree_scalar_selected <- c(
      if (isTRUE(tree_scalar_split_on_am1)) "am1",
      if (isTRUE(tree_scalar_split_on_am2)) "am2",
      if (isTRUE(tree_scalar_split_on_m0))  "m0"
    )
    
    tree_metric_correct_split <- if (length(true_metric_moderators) == 0) {
      length(tree_metric_selected) == 0
    } else {
      any(tree_metric_selected %in% true_metric_moderators)
    }
    
    tree_scalar_correct_split <- if (length(true_scalar_moderators) == 0) {
      length(tree_scalar_selected) == 0
    } else {
      any(tree_scalar_selected %in% true_scalar_moderators)
    }
  }

  
  # ---------------------------
  tibble(
    job_id         = as.integer(row$job_id), # row indentifier for simulation condition x replication
    popmodel       = as.character(row$popmodel), # population-generating model
    N              = as.integer(row$N), # sample size in the used generated data set
    reliability    = as.numeric(row$reliability), # target indicator reliability used during data generation
    lambda         = as.numeric(row$lambda), # baseline factor loading value
    intercepts     = as.numeric(row$intercepts), # baseline intercept value
    delta_lambda   = as.numeric(row$delta_lambda), # magnitude of loading moderation
    delta_nu       = as.numeric(row$delta_nu), # magnitude of intercept moderation
    moderator      = as.character(row$moderator), # functional form of moderation in the population
    analysis_form = as.character(row$analysis_form), # functional form assumed in analaysis
    rep_id = as.integer(row$rep_id), # replication number within a design condition
    
    # those show what is true in the dgm
    true_any_noninvariance    = as.logical(true_any_noninvariance), # true if either metric or scalar non-invariance is present
    true_metric_noninvariance = as.logical(true_metric_noninvariance), # true if either metric non-invariance is present
    true_scalar_noninvariance = as.logical(true_scalar_noninvariance), # true if either scalar non-invariance is present
    true_structured_moderator = as.logical(true_structured_moderator), #if the moderator contains systematic structure (!= "noise").
    
    # mnlfa statistics
    mnlfa_model    = as.character(mnlfa_model), # most restrictive successfully estimated MNLFA model
    mnlfa_det      = as.logical(mnlfa_det), # TRUE = some non-invariance detected; FALSE = invariance retained; NA = failed/inconclusive
    mnlfa_final_decision     = as.character(mnlfa_final_decision), # where noninvariance was detected
    
    mnlfa_metric_lrt_chisq   = as.numeric(mnlfa_metric_lrt_chisq),
    mnlfa_metric_lrt_df      = as.numeric(mnlfa_metric_lrt_df),
    mnlfa_metric_lrt_p       = as.numeric(mnlfa_metric_lrt_p),
    mnlfa_metric_lrt_reject  = as.logical(mnlfa_metric_lrt_reject),
    
    mnlfa_scalar_lrt_chisq   = as.numeric(mnlfa_scalar_lrt_chisq),
    mnlfa_scalar_lrt_df      = as.numeric(mnlfa_scalar_lrt_df),
    mnlfa_scalar_lrt_p       = as.numeric(mnlfa_scalar_lrt_p),
    mnlfa_scalar_lrt_reject  = as.logical(mnlfa_scalar_lrt_reject),
    
    # joint test of all moderation effects simultaneously
    mnlfa_omnibus_lrt_chisq  = as.numeric(mnlfa_omnibus_lrt_chisq),
    mnlfa_omnibus_lrt_df     = as.numeric(mnlfa_omnibus_lrt_df),
    mnlfa_omnibus_lrt_p      = as.numeric(mnlfa_omnibus_lrt_p),
    mnlfa_omnibus_lrt_reject = as.logical(mnlfa_omnibus_lrt_reject),
    
    # tree statistics (structural tree outcomes, not hypothesis-test decisions)
    tree_metric_split        = as.logical(tree_metric_split), # true if the metric tree produced at least one split
    tree_scalar_split        = as.logical(tree_scalar_split), # true if the scalar tree produced at least one split
    
    tree_metric_p            = as.numeric(tree_metric_p), # multiplicity-corrected global split-test p-value
    tree_metric_p_uncorrected = as.numeric(tree_metric_p_uncorrected),
    tree_metric_reject       = as.logical(tree_metric_reject), # true if corrected p-value ≤ alpha
    
    tree_scalar_p            = as.numeric(tree_scalar_p),
    tree_scalar_p_uncorrected = as.numeric(tree_scalar_p_uncorrected),
    tree_scalar_reject       = as.logical(tree_scalar_reject),
    
    # moderators used (these come from 'semtree_detects_moderation()')
    tree_metric_split_on_am1 = as.logical(tree_metric_split_on_am1), # true if any metric-tree split used 'am1'
    tree_metric_split_on_am2 = as.logical(tree_metric_split_on_am2), # true if any metric-tree split used 'am2'
    tree_metric_split_on_m0 = as.logical(tree_metric_split_on_m0), # true if any metric-tree split used 'm0'
    tree_metric_n_splits_am1 = as.integer(tree_metric_n_splits_am1), # number of metric-tree splits using 'am1'
    tree_metric_n_splits_am2 = as.integer(tree_metric_n_splits_am2), # number of metric-tree splits using 'am2'
    tree_metric_n_splits_m0 = as.integer(tree_metric_n_splits_m0), # number of metric-tree splits using 'm0'
    
    tree_scalar_split_on_am1 = as.logical(tree_scalar_split_on_am1),
    tree_scalar_split_on_am2 = as.logical(tree_scalar_split_on_am2),
    tree_scalar_split_on_m0 = as.logical(tree_scalar_split_on_m0),
    tree_scalar_n_splits_am1 = as.integer(tree_scalar_n_splits_am1),
    tree_scalar_n_splits_am2 = as.integer(tree_scalar_n_splits_am2),
    tree_scalar_n_splits_m0 = as.integer(tree_scalar_n_splits_m0),
    
    tree_metric_correct_split = as.logical(tree_metric_correct_split), # true if the metric tree selected at least one correct moderator 
    tree_scalar_correct_split = as.logical(tree_scalar_correct_split), # true if the scalar tree selected at least one correct moderator 

    error_msg = as.character(error_msg),
    mnlfa_error_msg = as.character(mnlfa_error_msg),
    semtree_error_msg = as.character(semtree_error_msg)
  ) %>%
  dplyr::bind_cols(mnlfa_mod_est)

  }



safe_run_one <- function(row) {
  tryCatch(
    run_one(row),
    error = function(e) {
      message("Error in job ", row$job_id, ": ", e$message)
      
      has_metric <- row$popmodel %in% c("1.1", "1.12", "1.2", "1.22", "1.3","1.32")
      has_scalar <- row$popmodel %in% c("1.11", "1.12", "1.21", "1.22","1.32")
      
      true_structured_moderator <- row$moderator != "noise"
      
      true_metric_noninvariance <- has_metric &&
        row$delta_lambda != 0 &&
        true_structured_moderator
      
      true_scalar_noninvariance <- has_scalar &&
        row$delta_nu != 0 &&
        true_structured_moderator
      
      true_any_noninvariance <- true_metric_noninvariance ||
        true_scalar_noninvariance
      
      tibble(
        job_id         = as.integer(row$job_id),
        popmodel       = as.character(row$popmodel),
        N              = as.integer(row$N),
        reliability    = as.numeric(row$reliability),
        lambda         = as.numeric(row$lambda),
        intercepts     = as.numeric(row$intercepts),
        delta_lambda   = as.numeric(row$delta_lambda),
        delta_nu       = as.numeric(row$delta_nu),
        moderator      = as.character(row$moderator),
        analysis_form = as.character(row$analysis_form),
        rep_id = as.integer(row$rep_id),
        
        true_any_noninvariance    = as.logical(true_any_noninvariance),
        true_metric_noninvariance = as.logical(true_metric_noninvariance),
        true_scalar_noninvariance = as.logical(true_scalar_noninvariance),
        true_structured_moderator = as.logical(true_structured_moderator),
        
        mnlfa_model = NA_character_,
        mnlfa_det = NA,
        
        mnlfa_final_decision = NA_character_,
        
        mnlfa_metric_lrt_chisq = NA_real_,
        mnlfa_metric_lrt_df = NA_real_,
        mnlfa_metric_lrt_p = NA_real_,
        mnlfa_metric_lrt_reject = NA,
        
        mnlfa_scalar_lrt_chisq = NA_real_,
        mnlfa_scalar_lrt_df = NA_real_,
        mnlfa_scalar_lrt_p = NA_real_,
        mnlfa_scalar_lrt_reject = NA,
        
        mnlfa_omnibus_lrt_chisq = NA_real_,
        mnlfa_omnibus_lrt_df = NA_real_,
        mnlfa_omnibus_lrt_p = NA_real_,
        mnlfa_omnibus_lrt_reject = NA,
        
        tree_metric_split = NA,
        tree_scalar_split = NA,
        
        tree_metric_p = NA_real_,
        tree_metric_p_uncorrected = NA_real_,
        tree_metric_reject = NA,
        
        tree_scalar_p = NA_real_,
        tree_scalar_p_uncorrected = NA_real_,
        tree_scalar_reject = NA,
        
        tree_metric_split_on_am1 = NA,
        tree_metric_split_on_am2 = NA,
        tree_metric_split_on_m0 = NA,
        tree_metric_n_splits_am1 = NA_integer_,
        tree_metric_n_splits_am2 = NA_integer_,
        tree_metric_n_splits_m0 = NA_integer_,
        
        tree_scalar_split_on_am1 = NA,
        tree_scalar_split_on_am2 = NA,
        tree_scalar_split_on_m0 = NA,
        tree_scalar_n_splits_am1 = NA_integer_,
        tree_scalar_n_splits_am2 = NA_integer_,
        tree_scalar_n_splits_m0 = NA_integer_,
        
        tree_metric_correct_split = NA,
        tree_scalar_correct_split = NA,
        
        mnlfa_error_msg = NA_character_,
        semtree_error_msg = NA_character_,
        error_msg = conditionMessage(e)
      )%>%
    dplyr::bind_cols(empty_mnlfa_moderation_estimates(p = 4))
    }
  )
}

#############################################

n_workers <- max(1, parallelly::availableCores() - 1)

plan(multisession, workers = n_workers)


## ----------- get splitter from command line argument ----------
## order the simulation to only round the i-th of j many chunks
## with i and j being the first two arguments

# get command line arguments
#args <- commandArgs(trailingOnly = TRUE)

#if (length(args)>0) {
  
#  if (length(args)==1) args[2]=100
  
#  chunk_id <- as.integer(args[1])
#  n_chunks <- as.integer(args[2]) 

#  all_indices <- seq_len(nrow(DESIGN))

#  chunks <- split(all_indices, cut(seq_along(all_indices),
#                                  n_chunks, labels = FALSE))
#  my_indices <- chunks[[chunk_id]]
  
#  DESIGN <- DESIGN[my_indices, ]

#}

# only for debugging purposes, run only first two rows:
# DESIGN <- DESIGN[1:4, ]


#
# -- Start Simulation --

t1 <- Sys.time()

# run across all rows (use future package's parallelization)
results <- future.apply::future_sapply(seq_len(nrow(DESIGN)), function(i) {
  safe_run_one(DESIGN[i, , drop = FALSE])
},simplify = TRUE)

results <- t(results)

t2 <- Sys.time()

elapsed_total_min <- as.numeric(difftime(t2, t1, units = "mins"))
elapsed_total_min

if (is.null(chunk_id)) {
  saveRDS(results, "results_parallel.rds")
} else {
  saveRDS(results, paste0("results_parallel_",chunk_id,"_of_",n_chunks,".rds"))
}

# this should be done later in a 
# collection script
#append_results(results, results_path)
########################################
