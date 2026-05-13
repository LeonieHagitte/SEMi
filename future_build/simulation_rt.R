library(furrr)
library(dplyr)
library(tibble)
library(future)


source("dataprep_rt.R")
source("analysis_rt.R")


manifest_path <- "progress_manifest.csv"
results_path  <- "results.csv"
lock_dir <- "locks"
# -----------------------------------
set.seed(42)
# ---------------------------------
n_rep <- 10

# ----------------------------------
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

manifest <- DESIGN %>%
  select(-seed) %>%
  mutate(
    status = "PENDING",
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
      analysis_form = character(),
      rep_id = integer(),
      
      true_any_noninvariance = logical(),
      true_metric_noninvariance = logical(),
      true_scalar_noninvariance = logical(),
      true_structured_moderator = logical(),
      
      mnlfa_model = character(),
      mnlfa_det = logical(),
      
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
      
      tree_metric_split_on_am1 = logical(),
      tree_metric_split_on_am2 = logical(),
      tree_metric_split_on_m0 = logical(),
      tree_metric_n_splits_am1 = integer(),
      tree_metric_n_splits_am2 = integer(),
      tree_metric_n_splits_m0 = integer(),
      
      tree_scalar_split_on_am1 = logical(),
      tree_scalar_split_on_am2 = logical(),
      tree_scalar_split_on_m0 = logical(),
      tree_scalar_n_splits_am1 = integer(),
      tree_scalar_n_splits_am2 = integer(),
      tree_scalar_n_splits_m0 = integer(),
      
      error_msg = character(),
      mnlfa_error_msg = character(),
      semtree_error_msg = character()
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
    analysis_form = as.character(row$analysis_form),
    rep_id = as.integer(row$rep_id),
    
    true_any_noninvariance    = as.logical(true_any_noninvariance),
    true_metric_noninvariance = as.logical(true_metric_noninvariance),
    true_scalar_noninvariance = as.logical(true_scalar_noninvariance),
    true_structured_moderator = as.logical(true_structured_moderator),
    
    mnlfa_model    = as.character(mnlfa_model),
    mnlfa_det      = as.logical(mnlfa_det),
    
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
    
    tree_metric_split_on_am1 = as.logical(tree_metric_split_on_am1),
    tree_metric_split_on_am2 = as.logical(tree_metric_split_on_am2),
    tree_metric_split_on_m0 = as.logical(tree_metric_split_on_m0),
    tree_metric_n_splits_am1 = as.integer(tree_metric_n_splits_am1),
    tree_metric_n_splits_am2 = as.integer(tree_metric_n_splits_am2),
    tree_metric_n_splits_m0 = as.integer(tree_metric_n_splits_m0),
    
    tree_scalar_split_on_am1 = as.logical(tree_scalar_split_on_am1),
    tree_scalar_split_on_am2 = as.logical(tree_scalar_split_on_am2),
    tree_scalar_split_on_m0 = as.logical(tree_scalar_split_on_m0),
    tree_scalar_n_splits_am1 = as.integer(tree_scalar_n_splits_am1),
    tree_scalar_n_splits_am2 = as.integer(tree_scalar_n_splits_am2),
    tree_scalar_n_splits_m0 = as.integer(tree_scalar_n_splits_m0),

    error_msg = as.character(error_msg),
    mnlfa_error_msg = as.character(mnlfa_error_msg),
    semtree_error_msg = as.character(semtree_error_msg)
  )
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
        
        mnlfa_error_msg = NA_character_,
        semtree_error_msg = NA_character_,
        error_msg = conditionMessage(e)
      )
    }
  )
}
# -------------------------------------------
done_jobs <- read.csv(results_path)$job_id

DESIGN <- DESIGN %>%
  filter(!job_id %in% done_jobs)
#############################################

plan(list(
  tweak(multisession, workers = length(parallelly::availableWorkers()) - 1),
  sequential
))



## ----------- get splitter from command line argument ----------
## order the simulation to only round the i-th of j many chunks
## with i and j being the first two arguments

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args)>0) {
  
  if (length(args)==1) args[2]=100
  
  chunk_id <- as.integer(args[1])
  n_chunks <- as.integer(args[2]) 

  all_indices <- seq_len(nrow(DESIGN))

  chunks <- split(all_indices, cut(seq_along(all_indices),
                                   n_chunks, labels = FALSE))
  my_indices <- chunks[[chunk_id]]
  
  DESIGN <- DESIGN[my_indices, ]

}

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
