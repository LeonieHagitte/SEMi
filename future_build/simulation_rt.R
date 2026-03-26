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
  seed          = 1:2
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
      
      tree_metric_split = logical(),
      tree_scalar_split = logical(),
      
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
    predictors = c("m1", "m2","m0") 
  )
  
  # ---------------------------
  mnlfa_metric_delta_cfi <- NA_real_
  mnlfa_metric_delta_rmsea <- NA_real_
  mnlfa_metric_retain <- NA
  
  mnlfa_scalar_delta_cfi <- NA_real_
  mnlfa_scalar_delta_rmsea <- NA_real_
  mnlfa_scalar_retain <- NA
  
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
  
  mnlfa_det <- FALSE
  if (identical(mnlfa_metric_retain, FALSE)) mnlfa_det <- TRUE
  if (identical(mnlfa_scalar_retain, FALSE)) mnlfa_det <- TRUE
  
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
  
  if (!inherits(res$semtree, "error")) {
    tree_metric_split <- res$semtree$metric_split
    tree_scalar_split <- res$semtree$scalar_split
    
    metric_info <- semtree_detects_moderation(
      res$semtree$metric_tree,
      moderators = c("m1", "m2","m0") 
    )
    
    scalar_info <- semtree_detects_moderation(
      res$semtree$scalar_tree,
      moderators = c("m1", "m2","m0") 
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
    job_id         = row$job_id,
    popmodel       = row$popmodel,
    N              = row$N,
    reliability    = row$reliability,
    lambda         = row$lambda,
    intercepts     = row$intercepts,
    delta_lambda   = row$delta_lambda,
    delta_nu       = row$delta_nu,
    moderator      = row$moderator,
    seed           = row$seed,
    
    mnlfa_model    = mnlfa_model,
    mnlfa_det      = mnlfa_det,
    
    tree_metric_split_on_m1 = tree_metric_split_on_m1,
    tree_metric_split_on_m2 = tree_metric_split_on_m2,
    tree_metric_n_splits_m1 = tree_metric_n_splits_m1,
    tree_metric_n_splits_m2 = tree_metric_n_splits_m2,
    
    tree_scalar_split_on_m1 = tree_scalar_split_on_m1,
    tree_scalar_split_on_m2 = tree_scalar_split_on_m2,
    tree_scalar_n_splits_m1 = tree_scalar_n_splits_m1,
    tree_scalar_n_splits_m2 = tree_scalar_n_splits_m2,
    
    mnlfa_metric_delta_cfi   = mnlfa_metric_delta_cfi,
    mnlfa_metric_delta_rmsea = mnlfa_metric_delta_rmsea,
    mnlfa_metric_retain      = mnlfa_metric_retain,
    
    mnlfa_scalar_delta_cfi   = mnlfa_scalar_delta_cfi,
    mnlfa_scalar_delta_rmsea = mnlfa_scalar_delta_rmsea,
    mnlfa_scalar_retain      = mnlfa_scalar_retain,
    
    mnlfa_final_decision     = mnlfa_final_decision,
    
    tree_metric_split     = tree_metric_split,
    tree_scalar_split     = tree_scalar_split
  )
  }

# ------- TEST --------------------
t1 <- Sys.time()
out <- run_one(DESIGN[1, ])
t2 <- Sys.time()

elapsed_one <- as.numeric(difftime(t2, t1, units = "secs"))
elapsed_one
out
str(out)

# ------- FULL RUN ----------------
t1 <- Sys.time()

all_results <- purrr::map_dfr(
  seq_len(nrow(DESIGN)),
  ~ run_one(DESIGN[.x, ])
)

t2 <- Sys.time()

elapsed_total_min <- as.numeric(difftime(t2, t1, units = "mins"))
elapsed_total_sec <- as.numeric(difftime(t2, t1, units = "secs"))

elapsed_total_min
elapsed_total_sec

append_results(all_results, results_path)
