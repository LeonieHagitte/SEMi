library(tidyverse)
library(furrr)
library(dplyr)
library(tibble)

source("dataprep_rt.R")
source("analysis_rt.R")


manifest_path <- "progress_manifest.csv"
results_path  <- "results.csv"
lock_dir <- "locks"


MOD_TYPES <- c("linear","sigmoid","quadratic","noise")

DESIGN <- tidyr::expand_grid(
  popmodel     = c("0","1.1", "1.11", "1.12","1.2","1.21","1.22","1.3"),
  N            = c(300, 500, 700, 1000),
  reliability  = c(0.60, 0.70, 0.80, 0.95),
  lambda       = 0.70,
  intercepts   = 1,
  # latentmean  = 0,
  # delta_eta   = c(-1, -0.5, 0.5, 1),
  delta_lambda = c(-0.3, -0.2, 0.2, 0.3),
  delta_nu     = c(-1, -0.5, 0.5, 1),
  moderator    = MOD_TYPES,
  seed          = 1:1
)%>%
  dplyr::arrange(
    popmodel, N, reliability, lambda, intercepts,
    delta_lambda, delta_nu, moderator, seed
  ) %>%
  dplyr::mutate(
    job_id = dplyr::row_number()
  )

manifest <- transform(
  DESIGN,
  status = "PENDING",           # PENDING | RUNNING | DONE | ERROR
  started_at = NA_character_,
  finished_at = NA_character_,
  worker = NA_character_,
  error_msg = NA_character_
)

if (!file.exists(manifest_path)) {
  write.csv(manifest, manifest_path, row.names = FALSE)
}
if (!dir.exists(lock_dir)) dir.create(lock_dir, recursive = TRUE)

if (!file.exists(results_path)) {
  write.csv(
    data.frame(
      job_id = integer(),
      popmodel = character(),
      N = integer(),
      reliability = numeric(),
      lambda = numeric(),
      intercepts = numeric(),
      delta_lambda = numeric(),
      delta_nu = numeric(),
      moderator = character(),
      seed = integer(),
      
      mnlfa_model = character(),
      mnlfa_det = logical(),
      
      mnlfa_metric_delta_cfi = numeric(),
      mnlfa_metric_delta_rmsea = numeric(),
      mnlfa_metric_retain = logical(),
      
      mnlfa_scalar_delta_cfi = numeric(),
      mnlfa_scalar_delta_rmsea = numeric(),
      mnlfa_scalar_retain = logical(),
      
      mnlfa_final_decision = character(),
      
      mnlfa_metric_lrt_chisq = numeric(),
      mnlfa_metric_lrt_df = numeric(),
      mnlfa_metric_lrt_p = numeric(),
      mnlfa_metric_lrt_reject = logical(),
      
      mnlfa_scalar_lrt_chisq = numeric(),
      mnlfa_scalar_lrt_df = numeric(),
      mnlfa_scalar_lrt_p = numeric(),
      mnlfa_scalar_lrt_reject = logical(),
      
      mnlfa_omnibus_lrt_chisq = numeric(),
      mnlfa_omnibus_lrt_df = numeric(),
      mnlfa_omnibus_lrt_p = numeric(),
      mnlfa_omnibus_lrt_reject = logical(),
      
      tree_metric_split = logical(),
      tree_scalar_split = logical(),
      
      tree_metric_p = numeric(),
      tree_metric_p_uncorrected = numeric(),
      tree_metric_reject = logical(),
      
      tree_scalar_p = numeric(),
      tree_scalar_p_uncorrected = numeric(),
      tree_scalar_reject = logical(),
      
      tree_metric_split_on_m1 = logical(),
      tree_metric_split_on_m2 = logical(),
      tree_metric_n_splits_m1 = integer(),
      tree_metric_n_splits_m2 = integer(),
      
      tree_scalar_split_on_m1 = logical(),
      tree_scalar_split_on_m2 = logical(),
      tree_scalar_n_splits_m1 = integer(),
      tree_scalar_n_splits_m2 = integer()
    ),
    results_path,
    row.names = FALSE
  )
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
  
  # ensure required columns exist for analyses
  if (!"m0" %in% names(df)) df$m0 <- 0
  
  # ---------------------------
  res <- run_analysis(
    data = df,
    methods = c("MNLFA", "SEMTREE"),
    nfactors = 1,
    alpha = 0.05,
    predictors = c("m1", "m2", "m0")
  )
  
  # ---------------------------
  mnlfa_metric_delta_cfi <- NA_real_
  mnlfa_metric_delta_rmsea <- NA_real_
  mnlfa_metric_retain <- NA
  
  mnlfa_scalar_delta_cfi <- NA_real_
  mnlfa_scalar_delta_rmsea <- NA_real_
  mnlfa_scalar_retain <- NA
  
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
    if (!is.null(res$mnlfa$metric_fit)) {
      mnlfa_metric_delta_cfi <- res$mnlfa$metric_fit$delta_cfi
      mnlfa_metric_delta_rmsea <- res$mnlfa$metric_fit$delta_rmsea
      mnlfa_metric_retain <- res$mnlfa$metric_fit$retain
    }
    if (!is.null(res$mnlfa$scalar_fit)) {
      mnlfa_scalar_delta_cfi <- res$mnlfa$scalar_fit$delta_cfi
      mnlfa_scalar_delta_rmsea <- res$mnlfa$scalar_fit$delta_rmsea
      mnlfa_scalar_retain <- res$mnlfa$scalar_fit$retain
    }
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
  
  if (identical(mnlfa_metric_retain, FALSE) || identical(mnlfa_scalar_retain, FALSE)) {
    mnlfa_det <- TRUE
  } else if (identical(mnlfa_metric_retain, TRUE) &&
             (identical(mnlfa_scalar_retain, TRUE) || is.na(mnlfa_scalar_retain))) {
    mnlfa_det <- FALSE
  }
  
  mnlfa_final_decision <- NA_character_
  
  if (identical(mnlfa_metric_retain, FALSE)) {
    mnlfa_final_decision <- "noninvariance_at_metric"
  } else if (identical(mnlfa_metric_retain, TRUE) &&
             identical(mnlfa_scalar_retain, FALSE)) {
    mnlfa_final_decision <- "noninvariance_at_scalar"
  } else if (identical(mnlfa_metric_retain, TRUE) &&
             (is.null(res$mnlfa$fitScalar) || identical(mnlfa_scalar_retain, TRUE))) {
    mnlfa_final_decision <- "scalar_invariance_retained"
  }
  
  # ---------------------------
  tree_metric_split <- NA
  tree_scalar_split <- NA
  
  tree_metric_split_on_m1 <- NA
  tree_metric_split_on_m2 <- NA
  tree_metric_n_splits_m1 <- NA_integer_
  tree_metric_n_splits_m2 <- NA_integer_
  
  tree_scalar_split_on_m1 <- NA
  tree_scalar_split_on_m2 <- NA
  tree_scalar_n_splits_m1 <- NA_integer_
  tree_scalar_n_splits_m2 <- NA_integer_
  
  tree_metric_p <- NA_real_
  tree_metric_p_uncorrected <- NA_real_
  tree_metric_reject <- NA
  
  tree_scalar_p <- NA_real_
  tree_scalar_p_uncorrected <- NA_real_
  tree_scalar_reject <- NA
  
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
      moderators = c("m1", "m2", "m0")
    )
    
    scalar_info <- semtree_detects_moderation(
      res$semtree$scalar_tree,
      moderators = c("m1", "m2", "m0")
    )
    
    tree_metric_split_on_m1 <- metric_info$tree_split_on_m1
    tree_metric_split_on_m2 <- metric_info$tree_split_on_m2
    tree_metric_n_splits_m1 <- metric_info$tree_n_splits_m1
    tree_metric_n_splits_m2 <- metric_info$tree_n_splits_m2
    
    tree_scalar_split_on_m1 <- scalar_info$tree_split_on_m1
    tree_scalar_split_on_m2 <- scalar_info$tree_split_on_m2
    tree_scalar_n_splits_m1 <- scalar_info$tree_n_splits_m1
    tree_scalar_n_splits_m2 <- scalar_info$tree_n_splits_m2
  }
  
  # ---------------------------
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
    seed           = as.integer(row$seed),
    
    mnlfa_model    = as.character(mnlfa_model),
    mnlfa_det      = as.logical(mnlfa_det),
    
    mnlfa_metric_delta_cfi   = as.numeric(mnlfa_metric_delta_cfi),
    mnlfa_metric_delta_rmsea = as.numeric(mnlfa_metric_delta_rmsea),
    mnlfa_metric_retain      = as.logical(mnlfa_metric_retain),
    
    mnlfa_scalar_delta_cfi   = as.numeric(mnlfa_scalar_delta_cfi),
    mnlfa_scalar_delta_rmsea = as.numeric(mnlfa_scalar_delta_rmsea),
    mnlfa_scalar_retain      = as.logical(mnlfa_scalar_retain),
    
    mnlfa_final_decision     = as.character(mnlfa_final_decision),
    
    mnlfa_metric_lrt_chisq   = as.numeric(mnlfa_metric_lrt_chisq),
    mnlfa_metric_lrt_df      = as.numeric(mnlfa_metric_lrt_df),
    mnlfa_metric_lrt_p       = as.numeric(mnlfa_metric_lrt_p),
    mnlfa_metric_lrt_reject  = as.logical(mnlfa_metric_lrt_reject),
    
    mnlfa_scalar_lrt_chisq   = as.numeric(mnlfa_scalar_lrt_chisq),
    mnlfa_scalar_lrt_df      = as.numeric(mnlfa_scalar_lrt_df),
    mnlfa_scalar_lrt_p       = as.numeric(mnlfa_scalar_lrt_p),
    mnlfa_scalar_lrt_reject  = as.logical(mnlfa_scalar_lrt_reject),
    
    mnlfa_omnibus_lrt_chisq  = as.numeric(mnlfa_omnibus_lrt_chisq),
    mnlfa_omnibus_lrt_df     = as.numeric(mnlfa_omnibus_lrt_df),
    mnlfa_omnibus_lrt_p      = as.numeric(mnlfa_omnibus_lrt_p),
    mnlfa_omnibus_lrt_reject = as.logical(mnlfa_omnibus_lrt_reject),
    
    tree_metric_split        = as.logical(tree_metric_split),
    tree_scalar_split        = as.logical(tree_scalar_split),
    
    tree_metric_p            = as.numeric(tree_metric_p),
    tree_metric_p_uncorrected = as.numeric(tree_metric_p_uncorrected),
    tree_metric_reject       = as.logical(tree_metric_reject),
    
    tree_scalar_p            = as.numeric(tree_scalar_p),
    tree_scalar_p_uncorrected = as.numeric(tree_scalar_p_uncorrected),
    tree_scalar_reject       = as.logical(tree_scalar_reject),
    
    tree_metric_split_on_m1 = as.logical(tree_metric_split_on_m1),
    tree_metric_split_on_m2 = as.logical(tree_metric_split_on_m2),
    tree_metric_n_splits_m1 = as.integer(tree_metric_n_splits_m1),
    tree_metric_n_splits_m2 = as.integer(tree_metric_n_splits_m2),
    
    tree_scalar_split_on_m1 = as.logical(tree_scalar_split_on_m1),
    tree_scalar_split_on_m2 = as.logical(tree_scalar_split_on_m2),
    tree_scalar_n_splits_m1 = as.integer(tree_scalar_n_splits_m1),
    tree_scalar_n_splits_m2 = as.integer(tree_scalar_n_splits_m2)
  )
  }

# ------- TEST --------------------

TEST_DESIGN <- tidyr::expand_grid(
  popmodel     = c("0", "1.22"),
  N            = c(300, 1000),
  reliability  = c(0.60, 0.95),
  lambda       = 0.70,
  intercepts   = 1,
  delta_lambda = c(-0.3, 0.3),
  delta_nu     = c(-1, 1),
  moderator    = c("linear", "noise"),
  seed         = 1
) %>%
  dplyr::arrange(
    popmodel, N, reliability, lambda, intercepts,
    delta_lambda, delta_nu, moderator, seed
  ) %>%
  dplyr::mutate(job_id = dplyr::row_number())

DESIGN <- TEST_DESIGN

safe_run_one <- function(row) {
  tryCatch(
    run_one(row),
    error = function(e) {
      message("Error in job ", row$job_id, ": ", e$message)
      return(NULL)
    }
  )
}

t1 <- Sys.time()
#out <- run_one(DESIGN[1, ])
results_list <- vector("list", nrow(DESIGN))

for (i in seq_len(nrow(DESIGN))) {
  cat("Running job", i, "of", nrow(DESIGN), "\n")
  
  res <- safe_run_one(DESIGN[i, ])
  
  if (!is.null(res)) {
    results_list[[i]] <- res
  }
}
t2 <- Sys.time()

elapsed_one <- as.numeric(difftime(t2, t1, units = "secs"))
elapsed_one

#out
#str(out)
append_results(dplyr::bind_rows(results_list), results_path)

# ------- FULL RUN ----------------
#t1 <- Sys.time()

#all_results <- purrr::map_dfr(
#  seq_len(nrow(DESIGN)),
#  ~ run_one(DESIGN[.x, ])
#)

#t2 <- Sys.time()

#elapsed_total_min <- as.numeric(difftime(t2, t1, units = "mins"))
#elapsed_total_sec <- as.numeric(difftime(t2, t1, units = "secs"))

#elapsed_total_min
#elapsed_total_sec

#append_results(all_results, results_path)




# -------------------------------
# TEST: single controlled run
# -------------------------------



#set.seed(123)

# ---- Define ONE condition ----
#test_row <- tibble(
#  job_id = 1,
#  popmodel = "1.22",   # contains moderation
#  N = 500,
#  reliability = 0.80,
#  lambda = 0.70,
#  intercepts = 1,
#  delta_lambda = 0.3,
#  delta_nu = 1,
#  moderator = "linear",
#  seed = 123
#)

# ---- Generate data ----
#params <- gen_paramsC(
#  popmodel      = test_row$popmodel,
#  lambda        = test_row$lambda,
#  nu            = test_row$intercepts,
#  reliability   = test_row$reliability,
#  moderator     = test_row$moderator,
#  delta_lambda  = test_row$delta_lambda,
#  delta_nu      = test_row$delta_nu
#)

#sim <- gen_dataC(
#  N = test_row$N,
#  params = params
#)

#df <- sim$data

# ensure required columns
#if (!"m0" %in% names(df)) df$m0 <- 0

# ---- Run analysis ----
#res <- run_analysis(
#  data = df,
#  methods = c("MNLFA", "SEMTREE"),
#  nfactors = 1,
#  alpha = 0.05,
#  predictors = c("m1", "m2", "m0")
#)