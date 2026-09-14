# ============================================================
# Sensitivity analysis:
# lower effect sizes + SD-matched moderator transformations
# ============================================================

library(dplyr)
library(tidyr)
library(future)
library(future.apply)
library(parallelly)

source("dataprep_rt.R")
source("analysis_rt.R")

set.seed(123)

# ------------------------------------------------------------
# Calibrated moderator transformation
# ------------------------------------------------------------

target_sd <- 1 / sqrt(3)

mod_h_cal <- function(M,
                      type = c("linear", "sigmoid", "quadratic", "noise"),
                      k = 10) {
  
  type <- match.arg(type)
  
  switch(
    type,
    
    linear = M,
    
    quadratic = {
      raw <- 2 * M^2 - 1
      
      mean_raw <- -1 / 3
      sd_raw   <- sqrt(16 / 45)
      
      (raw - mean_raw) / sd_raw * target_sd
    },
    
    sigmoid = {
      raw <- -1 + 2 / (1 + exp(-k * M))
      
      mean_raw <- 0
      sd_raw <- sqrt(1 - tanh(k / 2) / (k / 2))
      
      (raw - mean_raw) / sd_raw * target_sd
    },
    
    noise = rep(0, length(M))
  )
}
#------------------------------------------------------------------------------
gen_dataC_sensitivity <- function(N, params, return_latent = TRUE) {
  
  p <- 4
  
  lambda0   <- params$lambda
  nu0       <- params$nu
  rel       <- params$reliability
  moderator <- params$moderator
  k         <- params$k
  dlam      <- params$delta_lambda
  dnu       <- params$delta_nu
  deta      <- params$delta_eta
  psi       <- params$psi
  mu0       <- params$mu_eta
  popmodel  <- params$popmodel
  
  # person-level moderators
  eps <- .Machine$double.eps
  
  m1 <- runif(N, -1 + eps, 1 - eps)
  m2 <- runif(N, -1 + eps, 1 - eps)
  m0 <- runif(N, -1 + eps, 1 - eps)
  
  # IMPORTANT CHANGE:
  # use calibrated population transformations
  hm1  <- mod_h_cal(m1, type = moderator, k = k)
  hm2  <- mod_h_cal(m2, type = moderator, k = k)
  hm12 <- hm1 * hm2
  hm0  <- rep(0, N)
  
  # item moderation pattern
  dlam1_vec  <- rep(0, p)
  dnu1_vec   <- rep(0, p)
  dlam2_vec  <- rep(0, p)
  dnu2_vec   <- rep(0, p)
  dlam12_vec <- rep(0, p)
  dnu12_vec  <- rep(0, p)
  
  if (popmodel == "0" || popmodel == "NULL") {
    
    # no moderation
    
  } else if (popmodel == "1.1") {
    
    dlam1_vec[] <- dlam
    
  } else if (popmodel == "1.11") {
    
    dnu1_vec[] <- dnu
    
  } else if (popmodel == "1.12") {
    
    dlam1_vec[] <- dlam
    dnu1_vec[]  <- dnu
    
  } else if (popmodel == "1.2") {
    
    dlam1_vec[1:2] <- dlam
    
  } else if (popmodel == "1.21") {
    
    dnu1_vec[1:2] <- dnu
    
  } else if (popmodel == "1.22") {
    
    dlam1_vec[1:2] <- dlam
    dnu1_vec[1:2]  <- dnu
    
  } else if (popmodel == "1.3") {
    
    dlam1_vec[]  <- dlam
    dlam2_vec[]  <- dlam
    dlam12_vec[] <- dlam^2
    
  } else if (popmodel == "1.32") {
    
    dlam1_vec[] <- dlam
    dnu2_vec[]  <- dnu
    
  } else {
    
    stop("Unknown popmodel: ", popmodel)
  }
  
  # person-specific loadings/intercepts
  Lambda_Np <- matrix(lambda0, N, p) +
    hm1  %o% dlam1_vec +
    hm2  %o% dlam2_vec +
    hm12 %o% dlam12_vec
  
  Nu_Np <- matrix(nu0, N, p) +
    hm1  %o% dnu1_vec +
    hm2  %o% dnu2_vec +
    hm12 %o% dnu12_vec
  
  # latent variable
  MU_eta_i <- rep(mu0, N)
  eta <- rnorm(N, mean = MU_eta_i, sd = sqrt(psi))
  
  # residual variance
  theta0 <- (lambda0^2 * psi * (1 - rel)) / rel
  
  ThetaVar_Np <- matrix(
    theta0,
    nrow = N,
    ncol = p,
    byrow = TRUE
  )
  
  E <- matrix(rnorm(N * p), N, p) * sqrt(ThetaVar_Np)
  
  X <- Nu_Np + Lambda_Np * eta + E
  
  colnames(X) <- paste0("x", 1:p)
  
  data <- as.data.frame(X)
  
  data$m1   <- m1
  data$hm1  <- hm1
  data$m2   <- m2
  data$hm2  <- hm2
  data$hm12 <- hm12
  data$m0   <- m0
  data$hm0  <- hm0
  
  out <- list(
    data = data,
    params = params
  )
  
  if (return_latent) {
    out$eta <- eta
  }
  
  out
}
#-------------------------------------------------------------------------------

n_rep <- 1000

SENSITIVITY_DESIGN <- tidyr::expand_grid(
  popmodel     = c("1.1", "1.11", "1.12"),
  N            = c(500, 1000),
  reliability  = 0.75,
  lambda       = 0.70,
  intercepts   = 1,
  delta_lambda = 0.10,
  delta_nu     = 0.25,
  moderator    = c("linear", "quadratic", "sigmoid"),
  method       = c("SEMTREE", "MNLFA", "MNLFAQ"),
  rep_id       = seq_len(n_rep),
  num_noisy_predictors = 0
) %>%
  mutate(
    sensitivity_case = "low_effect_sd_matched",
    job_id = row_number(),
    seed = sample.int(.Machine$integer.max, n(), replace = TRUE)
  )


#-------------------------------------------------------------------------------
mnlfa_moderation_estimate_names <- function(p = 4) {
  base_names <- c(
    "mnlfa_est_dnu_am1",
    "mnlfa_est_dnu_am2",
    "mnlfa_est_dnu_am12",
    "mnlfa_est_dlambda_am1",
    "mnlfa_est_dlambda_am2",
    "mnlfa_est_dlambda_am12"
  )
  
  as.vector(
    outer(base_names, paste0("x", seq_len(p)), paste, sep = "_")
  )
}


empty_mnlfa_moderation_estimates <- function(p = 4) {
  
  cols <- mnlfa_moderation_estimate_names(p = p)
  
  tibble::as_tibble(
    as.list(
      stats::setNames(
        rep(NA_real_, length(cols)),
        cols
      )
    )
  )
}


flatten_mnlfa_moderation_estimates <- function(mnlfa_result, p = 4) {
  
  out <- empty_mnlfa_moderation_estimates(p = p)
  
  if (inherits(mnlfa_result, "error") ||
      is.null(mnlfa_result$configural_moderation_estimates)) {
    return(out)
  }
  
  est <- mnlfa_result$configural_moderation_estimates
  
  est_wide <- est |>
    tidyr::pivot_wider(
      names_from = item,
      values_from = c(
        mnlfa_est_dnu_am1,
        mnlfa_est_dnu_am2,
        mnlfa_est_dnu_am12,
        mnlfa_est_dlambda_am1,
        mnlfa_est_dlambda_am2,
        mnlfa_est_dlambda_am12
      ),
      names_glue = "{.value}_{item}"
    )
  
  out[names(out)] <- est_wide[names(out)]
  
  out
}


get_tree_predictors <- function(
    popmodel,
    noisy_predictor_names = character(0)) {
  
  base_predictors <- switch(
    as.character(popmodel),
    
    "0"    = "am1",
    "1.1"  = "am1",
    "1.11" = "am1",
    "1.12" = "am1",
    "1.2"  = "am1",
    "1.21" = "am1",
    "1.22" = "am1",
    "1.3"  = c("am1", "am2"),
    "1.32" = c("am1", "am2"),
    
    stop("Unknown popmodel: ", popmodel)
  )
  
  unique(c(
    base_predictors,
    noisy_predictor_names
  ))
}
#-------------------------------------------------------------------------------
run_one_sens <- function(row) { #run_one_sens <- function(seed, N, popmodel, moderator) 
  
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
  
  sim <- gen_dataC_sensitivity(
    N = row$N,
    params = params
  )
  
  df <- sim$data
  
  analysis_form <- if (row$method == "MNLFAQ") {
    "quadratic"
  } else if (row$method %in% c("MNLFA", "SEMTREE")) {
    "linear"
  } else {
    stop("Unknown method: ", row$method)
  }
  
  df <- add_analysis_form(
    data = df,
    analysis_form = analysis_form,
    k = params$k
  )
  
  temp_colnames <- colnames(df)
  
  df <- add_noisy_predictors(
    data = df,
    num_noisy_predictors = row$num_noisy_predictors
  )
  
  noisy_predictor_names <- setdiff(colnames(df), temp_colnames)
  
  tree_predictors <- get_tree_predictors(
    popmodel = row$popmodel,
    noisy_predictor_names = if (row$method == "SEMTREE") noisy_predictor_names else character(0)
  )
  
  tree_moderators_to_check <- c("am1", "am2", "m0", noisy_predictor_names)
  
  runtime_start <- Sys.time()
  
  res <- run_analysis(
    data = df,
    methods = row$method,
    nfactors = 1,
    alpha = 0.05,
    tree_predictors = tree_predictors
  )
  
  runtime_end <- Sys.time()
  
  runtime_sec <- as.numeric(
    difftime(runtime_end, runtime_start, units = "secs")
  )
  
  mnlfa_result <- if (row$method == "MNLFAQ") {
    res$mnlfaq
  } else if (row$method == "MNLFA") {
    res$mnlfa
  } else {
    NULL
  }
  
  # truth indicators
  # ---------------------------
  has_metric <- row$popmodel %in% c("1.1", "1.12", "1.2", "1.22", "1.3", "1.32")
  has_scalar <- row$popmodel %in% c("1.11", "1.12", "1.21", "1.22", "1.32")
  
  # ---------------------------
  mnlfa_error_msg <- NA_character_
  semtree_error_msg <- NA_character_
  
  if (inherits(mnlfa_result, "error")) {
    mnlfa_error_msg <- conditionMessage(mnlfa_result)
  }
  
  if (inherits(res$semtree, "error")) {
    semtree_error_msg <- conditionMessage(res$semtree)
  }
  
  
  mnlfa_mod_est <- flatten_mnlfa_moderation_estimates(
    mnlfa_result = mnlfa_result,
    p = 4
  )
  
  # ---------------------------
  mnlfa_kl_configural <- NA_real_
  mnlfa_kl_metric <- NA_real_
  mnlfa_kl_scalar <- NA_real_
  
  #if (!inherits(mnlfa_result, "error")) {                          # commented out for runtime improvement
  #  mnlfa_kl_configural <- average_kl_mnlfa(df, params, mnlfa_result$fitConfig)
  #  mnlfa_kl_metric     <- average_kl_mnlfa(df, params, mnlfa_result$fitMetric)
  #  mnlfa_kl_scalar     <- average_kl_mnlfa(df, params, mnlfa_result$fitScalar)
  #}
  
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
  
  if (!is.null(mnlfa_result) && !inherits(mnlfa_result, "error")) {
    if (!is.null(mnlfa_result$metric_lrt)) {
      mnlfa_metric_lrt_chisq <- mnlfa_result$metric_lrt$chisq_diff
      mnlfa_metric_lrt_df <- mnlfa_result$metric_lrt$df_diff
      mnlfa_metric_lrt_p <- mnlfa_result$metric_lrt$p_value
      mnlfa_metric_lrt_reject <- mnlfa_result$metric_lrt$reject_h0
    }
    
    if (!is.null(mnlfa_result$scalar_lrt)) {
      mnlfa_scalar_lrt_chisq <- mnlfa_result$scalar_lrt$chisq_diff
      mnlfa_scalar_lrt_df <- mnlfa_result$scalar_lrt$df_diff
      mnlfa_scalar_lrt_p <- mnlfa_result$scalar_lrt$p_value
      mnlfa_scalar_lrt_reject <- mnlfa_result$scalar_lrt$reject_h0
    }
    
    if (!is.null(mnlfa_result$omnibus_lrt)) {
      mnlfa_omnibus_lrt_chisq <- mnlfa_result$omnibus_lrt$chisq_diff
      mnlfa_omnibus_lrt_df <- mnlfa_result$omnibus_lrt$df_diff
      mnlfa_omnibus_lrt_p <- mnlfa_result$omnibus_lrt$p_value
      mnlfa_omnibus_lrt_reject <- mnlfa_result$omnibus_lrt$reject_h0
    }
  }
  
  mnlfa_model <- NA_character_
  if (!is.null(mnlfa_result) && !inherits(mnlfa_result, "error")) {
    if (!is.null(mnlfa_result$fitScalar)) {
      mnlfa_model <- "scalar"
    } else if (!is.null(mnlfa_result$fitMetric)) {
      mnlfa_model <- "metric"
    } else if (!is.null(mnlfa_result$fitConfig)) {
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
  
  tree_scalar_tested_anyway <- FALSE
  tree_scalar_power_anyway  <- NA
  tree_sequential_decision <- NA_character_
  tree_scalar_p <- NA_real_
  tree_scalar_p_uncorrected <- NA_real_
  tree_scalar_reject <- NA
  
  tree_metric_correct_split <- NA
  tree_scalar_correct_split <- NA
  tree_metric_split_on_noisy <- NA
  tree_scalar_split_on_noisy <- NA
  
  if (!is.null(res$semtree) && !inherits(res$semtree, "error")) {
    
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
    
    tree_sequential_decision <- res$semtree$tree_sequential_decision
    
    tree_scalar_tested_anyway <- !is.na(tree_scalar_p)
    tree_scalar_power_anyway  <- tree_scalar_reject
    
    metric_info <- semtree_detects_moderation(
      res$semtree$metric_tree,
      moderators = tree_moderators_to_check
    )
    
    scalar_info <- semtree_detects_moderation(
      res$semtree$scalar_tree,
      moderators = tree_moderators_to_check
    )
    
    tree_metric_split_on_noisy <- if (length(noisy_predictor_names) == 0) {
      FALSE
    } else {
      any(unlist(metric_info[paste0("tree_split_on_", noisy_predictor_names)]), na.rm = TRUE)
    }
    
    tree_scalar_split_on_noisy <- if (length(noisy_predictor_names) == 0) {
      FALSE
    } else {
      any(unlist(scalar_info[paste0("tree_split_on_", noisy_predictor_names)]), na.rm = TRUE)
    }
    
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
    analysis_form = as.character(analysis_form), # functional form assumed in analaysis
    method = as.character(row$method),
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
    
    # KL divergence between true DGM-implied and fitted MNLFA-implied moments
    mnlfa_kl_configural = as.numeric(mnlfa_kl_configural),
    mnlfa_kl_metric     = as.numeric(mnlfa_kl_metric),
    mnlfa_kl_scalar     = as.numeric(mnlfa_kl_scalar),
    
    # tree statistics (structural tree outcomes, not hypothesis-test decisions)
    tree_metric_split        = as.logical(tree_metric_split), # true if the metric tree produced at least one split
    tree_scalar_split        = as.logical(tree_scalar_split), # true if the scalar tree produced at least one split
    
    tree_metric_p            = as.numeric(tree_metric_p), # multiplicity-corrected global split-test p-value
    tree_metric_p_uncorrected = as.numeric(tree_metric_p_uncorrected),
    tree_metric_reject       = as.logical(tree_metric_reject), # true if corrected p-value ≤ alpha
    
    tree_scalar_p            = as.numeric(tree_scalar_p),
    tree_scalar_p_uncorrected = as.numeric(tree_scalar_p_uncorrected),
    tree_scalar_reject       = as.logical(tree_scalar_reject),
    
    tree_sequential_decision = as.character(tree_sequential_decision),
    tree_scalar_tested_anyway = as.logical(tree_scalar_tested_anyway),
    tree_scalar_power_anyway = as.logical(tree_scalar_power_anyway),
    
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
    tree_metric_split_on_noisy = as.logical(tree_metric_split_on_noisy),
    tree_scalar_split_on_noisy = as.logical(tree_scalar_split_on_noisy),
    num_noisy_predictors = as.integer(row$num_noisy_predictors),
    runtime_sec = as.numeric(runtime_sec),
    
    error_msg = as.character(error_msg),
    mnlfa_error_msg = as.character(mnlfa_error_msg),
    semtree_error_msg = as.character(semtree_error_msg)
  ) %>%
    dplyr::bind_cols(mnlfa_mod_est)
  
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#TEST_DESIGN <- SENSITIVITY_DESIGN %>%
#  filter(
#    popmodel == "1.1",
#    N == 500,
#    moderator == "linear"
#  ) %>%
#  slice_head(n = 3)
#
#test_results <- lapply(
#  seq_len(nrow(TEST_DESIGN)),
#  function(i) {
#    run_one_sens(TEST_DESIGN[i, ])
#  }
#) %>%
#  bind_rows()

#test_results %>%
#  select(
#    popmodel,
#    N,
#    moderator,
#    method,
#    mnlfa_metric_lrt_reject,
#    tree_metric_reject,
#    error_msg,
#    mnlfa_error_msg,
#    semtree_error_msg
#  )
#----------------------------------------
#TEST_DESIGN1 <- SENSITIVITY_DESIGN %>%
#  filter(
#    popmodel == "1.1",
#    N == 500,
#    moderator == "linear",
#    rep_id == 1
#  )

#test_results <- lapply(
#  seq_len(nrow(TEST_DESIGN1)),
#  function(i) {
#    run_one_sens(TEST_DESIGN1[i, ])
#  }
#) %>%
#  bind_rows()

#test_results %>%
#  select(
#    popmodel,
#    N,
#    moderator,
#    method,
#    analysis_form,
#    mnlfa_metric_lrt_reject,
#    mnlfa_scalar_lrt_reject,
#    tree_metric_reject,
#    tree_scalar_reject,
#    error_msg,
#    mnlfa_error_msg,
#    semtree_error_msg
#  )
#----------------------------------------------------
#TEST_DESIGN <- SENSITIVITY_DESIGN %>%
#  filter(
#    popmodel == "1.1",
#    N == 500,
#    rep_id == 1
#  )
#-------------------------------------------------------------------------------

n_workers <- max(1, parallelly::availableCores() - 1)

plan(multisession, workers = n_workers)

# -- Start Sensitivity Simulation --

t1 <- Sys.time()

results_sensitivity <- future.apply::future_sapply(
  seq_len(nrow(SENSITIVITY_DESIGN)),
  function(i) {
    run_one_sens(SENSITIVITY_DESIGN[i, , drop = FALSE])
  },
  simplify = TRUE
)

results_sensitivity <- t(results_sensitivity)

t2 <- Sys.time()

elapsed_total_min <- as.numeric(
  difftime(t2, t1, units = "mins")
)

elapsed_total_min

saveRDS(
  results_sensitivity,
  "results_sensitivity_parallel.rds"
)