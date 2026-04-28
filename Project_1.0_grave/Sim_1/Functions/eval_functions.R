# eval_functions.R
# =============================================================================
# Purpose
#   Unified evaluation helpers for simulation studies with lavaan and mxsem:
#   - Generic kernels: alignment, bias, RMSE, coverage, MCSE.
#   - lavaan extractors: loadings, residual variances (Theta), deltas-by-label.
#   - mxsem/OpenMx extractors: parameter table + delta pull by name.
#   - Backend-agnostic dispatchers for δ (moderation) metrics.
# =============================================================================

# ---- dependencies ------------------------------------------------------------
library(dplyr)
library(purrr)
library(lavaan)   
library(OpenMx)

# Align two named numeric vectors on common names
.align_by_name <- function(est, tru) {
  common <- intersect(names(est), names(tru))
  list(est = est[common], tru = tru[common])
}

# Mean relative bias across parameters (zero-safe: if true==0, denom=1)
mean_relative_bias <- function(est_vec, true_vec) {
  if (length(est_vec) == 0) return(NA_real_)
  al <- .align_by_name(est_vec, true_vec)
  if (length(al$est) == 0) return(NA_real_)
  denom <- ifelse(al$tru == 0, 1, al$tru)
  mean((al$est - al$tru) / denom, na.rm = TRUE)
}

# Relative RMSE (useful for non-zero-anchored targets like loadings/Theta)
relative_rmse <- function(est_vec, true_vec) {
  if (length(est_vec) == 0) return(NA_real_)
  al <- .align_by_name(est_vec, true_vec)
  if (length(al$est) == 0) return(NA_real_)
  sqrt(mean((al$est - al$tru)^2, na.rm = TRUE)) /
    mean(al$tru, na.rm = TRUE)
}

# Absolute RMSE (recommended for deltas, which are often 0 under NULL)
rmse_absolute <- function(est_vec, true_vec) {
  if (length(est_vec) == 0) return(NA_real_)
  al <- .align_by_name(est_vec, true_vec)
  if (length(al$est) == 0) return(NA_real_)
  sqrt(mean((al$est - al$tru)^2, na.rm = TRUE))
}

# Generic coverage calculator for a parameter table that already contains
#   columns: ci.lower, ci.upper, and a key "namekey" mapping to true_vec names.
.coverage_from_tab <- function(tab, name_key, true_vec) {
  if (nrow(tab) == 0) return(NA_real_)
  tab$namekey <- name_key
  ind <- map_lgl(names(true_vec), function(nm) {
    row <- tab %>% filter(namekey == nm)
    if (nrow(row) == 0) return(NA)
    any((row$ci.lower <= true_vec[[nm]]) & (true_vec[[nm]] <= row$ci.upper))
  })
  mean(ind, na.rm = TRUE)
}

# Monte Carlo SEs for a vector of per-rep metrics
calculate_mcse_bias <- function(bias_list) {
  x <- na.omit(bias_list); if (!length(x)) return(NA_real_)
  sqrt(var(x) / length(x))
}
calculate_mcse_rmse <- function(rmse_list) {
  x <- na.omit(rmse_list); if (!length(x)) return(NA_real_)
  K <- length(x); m <- mean(x)
  sqrt(sum((x - m)^2) / (K * (K - 1)))
}

# =============================================================================
# 1) lavaan extractors and wrappers
# =============================================================================

# Extract loadings as a named vector "f1=~y1", ...
get_estimated_loadings <- function(fit) {
  if (is.null(fit) || !lavInspect(fit, "converged")) return(setNames(numeric(0), character(0)))
  pe <- parameterEstimates(fit, ci = TRUE) %>% filter(op == "=~")
  setNames(pe$est, paste0(pe$lhs, "=~", pe$rhs))
}

# Extract residual variances (Theta diagonals) as "y1~~y1", ...
get_estimated_theta_diag <- function(fit, ov_prefix = "y") {
  if (is.null(fit) || !lavInspect(fit, "converged")) return(setNames(numeric(0), character(0)))
  pe <- parameterEstimates(fit, ci = TRUE) %>%
    filter(op == "~~", lhs == rhs, grepl(paste0("^", ov_prefix), lhs))
  setNames(pe$est, paste0(pe$lhs, "~~", pe$rhs))
}

# Extract moderation deltas by label prefix (e.g., "delta_")
# Requires LABELing the relevant parameters in lavaan syntax.
get_estimated_deltas_by_label <- function(fit, label_prefix = "delta_") {
  if (is.null(fit) || !lavInspect(fit, "converged")) return(setNames(numeric(0), character(0)))
  pe <- parameterEstimates(fit, ci = TRUE) %>%
    filter(!is.na(label), startsWith(label, label_prefix))
  setNames(pe$est, pe$label)
}

# Coverage (lavaan): loadings
calculate_coverage_loadings <- function(fit, true_loadings) {
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NA_real_)
  pe <- parameterEstimates(fit, ci = TRUE) %>% filter(op == "=~")
  key <- paste0(pe$lhs, "=~", pe$rhs)
  .coverage_from_tab(pe, key, true_loadings)
}

# Coverage (lavaan): Theta diagonals
calculate_coverage_theta <- function(fit, true_theta, ov_prefix = "y") {
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NA_real_)
  pe <- parameterEstimates(fit, ci = TRUE) %>%
    filter(op == "~~", lhs == rhs, grepl(paste0("^", ov_prefix), lhs))
  key <- paste0(pe$lhs, "~~", pe$rhs)
  .coverage_from_tab(pe, key, true_theta)
}

# Coverage (lavaan): deltas by label
calculate_coverage_deltas <- function(fit, true_deltas, label_prefix = "delta_") {
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NA_real_)
  pe <- parameterEstimates(fit, ci = TRUE) %>%
    filter(!is.na(label), startsWith(label, label_prefix))
  key <- pe$label
  .coverage_from_tab(pe, key, true_deltas)
}

# Bias/RMSE wrappers (lavaan): loadings, Theta, deltas
calculate_relative_bias_loadings <- function(fit, true_loadings) {
  mean_relative_bias(get_estimated_loadings(fit), true_loadings)
}
calculate_relative_rmse_loadings <- function(fit, true_loadings) {
  relative_rmse(get_estimated_loadings(fit), true_loadings)
}
calculate_relative_bias_theta <- function(fit, true_theta, ov_prefix = "y") {
  mean_relative_bias(get_estimated_theta_diag(fit, ov_prefix), true_theta)
}
calculate_relative_rmse_theta <- function(fit, true_theta, ov_prefix = "y") {
  relative_rmse(get_estimated_theta_diag(fit, ov_prefix), true_theta)
}
calculate_relative_bias_deltas <- function(fit, true_deltas, label_prefix = "delta_") {
  mean_relative_bias(get_estimated_deltas_by_label(fit, label_prefix), true_deltas)
}
calculate_relative_rmse_deltas <- function(fit, true_deltas, label_prefix = "delta_") {
  relative_rmse(get_estimated_deltas_by_label(fit, label_prefix), true_deltas)
}

# =============================================================================

# -------------------------------------------------------------------
# Robust extractor for mxsem/OpenMx: estimates, SEs, and CI bounds
# Returns tibble with columns: name, estimate, se, ci.lower, ci.upper
# -------------------------------------------------------------------
mxsem_param_table <- function(fit) {
  sumx <- summary(fit)
  
  # --- 1) parameter estimates (+ SE), standardize names -----------------------
  par <- as.data.frame(sumx$parameters)
  if (!"Estimate"  %in% names(par)) stop("mxsem_param_table: Estimate not found in OpenMx summary.")
  if (!"Std.Error" %in% names(par)) par$Std.Error <- NA_real_
  
  # stable 'name' column
  if (!"name" %in% names(par)) {
    if      ("label"     %in% names(par)) par$name <- par$label
    else if ("Parameter" %in% names(par)) par$name <- par$Parameter
    else                                  par$name <- rownames(par)
  }
  
  par <- dplyr::rename(par, estimate = Estimate, se = Std.Error)
  par <- dplyr::select(par, name, estimate, se)
  
  # --- 2) CI table (may be missing or differently named) ----------------------
  ci_raw <- try(as.data.frame(sumx$CI), silent = TRUE)
  have_ci <- !(inherits(ci_raw, "try-error") || is.null(ci_raw) || !nrow(ci_raw))
  
  if (!have_ci) {
    ci <- data.frame(name = character(0), ci.lower = numeric(0), ci.upper = numeric(0))
  } else {
    ci <- ci_raw
    
    # ensure 'name'
    if (!"name" %in% names(ci)) {
      if      ("param"     %in% names(ci)) ci$name <- ci$param
      else if ("Parameter" %in% names(ci)) ci$name <- ci$Parameter
      else if ("label"     %in% names(ci)) ci$name <- ci$label
      else if (!is.null(rownames(ci)))     ci$name <- rownames(ci)
      else                                  ci$name <- NA_character_
    }
    
    # normalize lower/upper bound names into ci.lower / ci.upper
    has <- function(x) x %in% names(ci)
    
    if (!has("ci.lower")) {
      if      (has("lbound"))   ci$ci.lower <- ci$lbound
      else if (has("lower"))    ci$ci.lower <- ci$lower
      else if (has("Lower"))    ci$ci.lower <- ci$Lower
      else if (has("`2.5%`"))   ci$ci.lower <- ci[["`2.5%`"]]
      else                      ci$ci.lower <- NA_real_
    }
    if (!has("ci.upper")) {
      if      (has("ubound"))   ci$ci.upper <- ci$ubound
      else if (has("upper"))    ci$ci.upper <- ci$upper
      else if (has("Upper"))    ci$ci.upper <- ci$Upper
      else if (has("`97.5%`"))  ci$ci.upper <- ci[["`97.5%`"]]
      else                      ci$ci.upper <- NA_real_
    }
    
    ci <- dplyr::select(ci, name, ci.lower, ci.upper)
    ci$ci.lower <- suppressWarnings(as.numeric(ci$ci.lower))
    ci$ci.upper <- suppressWarnings(as.numeric(ci$ci.upper))
  }
  
  # --- 3) join estimates with CIs ---------------------------------------------
  out <- dplyr::left_join(par, ci, by = "name")
  if (!"ci.lower" %in% names(out)) out$ci.lower <- NA_real_
  if (!"ci.upper" %in% names(out)) out$ci.upper <- NA_real_
  out
}


# Extract only delta parameters by their exact names.
# delta_names: character vector, e.g., names(sim$true_deltas)
extract_deltas_mxsem <- function(fit, delta_names) {
  tab <- mxsem_param_table(fit)
  wanted <- dplyr::filter(tab, .data$name %in% delta_names)
  list(
    est   = setNames(wanted$estimate,  wanted$name),
    se    = setNames(wanted$se,        wanted$name),
    lo    = setNames(wanted$ci.lower,  wanted$name),
    hi    = setNames(wanted$ci.upper,  wanted$name),
    table = wanted
  )
}

# =============================================================================
# 3) Backend-agnostic dispatchers for δ (use in the driver)
# =============================================================================

# Backend tests
.is_lavaan <- function(fit) inherits(fit, "lavaan")
.is_mxsem  <- function(fit) inherits(fit, "MxModel") || inherits(fit, "MxRAMModel")

# Named vector of δ estimates from either backend.
# For lavaan, it pulls by label prefix; for mxsem, it matches exact names.
get_delta_estimates_any <- function(fit, delta_names, label_prefix = "delta_") {
  if (.is_mxsem(fit)) {
    ex <- extract_deltas_mxsem(fit, delta_names)
    return(ex$est)
  } else if (.is_lavaan(fit)) {
    pe  <- lavaan::parameterEstimates(fit, ci = TRUE)
    sel <- pe[!is.na(pe$label) & startsWith(pe$label, label_prefix), ]
    vec <- setNames(sel$est, sel$label)
    return(vec[names(vec) %in% delta_names])
  } else {
    setNames(numeric(0), character(0))
  }
}

# δ coverage from either backend (returns mean coverage; NA if CIs absent).
# For OpenMx, must run with intervals=TRUE to populate CI bounds.
coverage_deltas_any <- function(fit, true_delta, label_prefix = "delta_") {
  if (.is_mxsem(fit)) {
    tab <- mxsem_param_table(fit)
    tab <- tab[tab$name %in% names(true_delta), ]
    if (!nrow(tab)) return(NA_real_)
    mean((tab$ci.lower <= true_delta[tab$name]) &
           (true_delta[tab$name] <= tab$ci.upper), na.rm = TRUE)
  } else if (.is_lavaan(fit)) {
    pe <- lavaan::parameterEstimates(fit, ci = TRUE)
    pe <- pe[!is.na(pe$label) & pe$label %in% names(true_delta), ]
    if (!nrow(pe)) return(NA_real_)
    mean((pe$ci.lower <= true_delta[pe$label]) &
           (true_delta[pe$label] <= pe$ci.upper), na.rm = TRUE)
  } else {
    NA_real_
  }
}

# Convenience: compute δ metrics (bias, abs RMSE, coverage) in one call.
# Returns a named list; coverage is NA if CIs are not available.
delta_metrics_any <- function(fit, true_delta, label_prefix = "delta_") {
  
  dn <- names(true_delta)
  est <- get_delta_estimates_any(fit, dn, label_prefix)
  list(
    bias     = mean_relative_bias(est, true_delta),
    rmse_abs = rmse_absolute(est, true_delta),
    coverage = coverage_deltas_any(fit, true_delta, label_prefix)
  )
}

# ------- MNLFA moderation detector (mxsem/OpenMx) ----------------------------
# Returns:
#   list(any = TRUE/FALSE/NA, n_sig = integer, per_param = named logical vector)
# Uses Wald z-tests when SEs are available; otherwise returns NA.
detect_moderation_mxsem <- function(fit, label_prefix = "delta_", alpha = 0.05) {
  if (!(inherits(fit, "MxModel") || inherits(fit, "MxRAMModel"))) {
    return(list(any = NA, n_sig = NA_integer_, per_param = NULL))
  }
  tab <- tryCatch(mxsem_param_table(fit), error = function(e) NULL)
  if (is.null(tab) || !all(c("name","estimate","se") %in% names(tab))) {
    return(list(any = NA, n_sig = NA_integer_, per_param = NULL))
  }
  del <- dplyr::filter(tab, !is.na(.data$name), startsWith(.data$name, label_prefix))
  if (nrow(del) == 0) {
    return(list(any = NA, n_sig = NA_integer_, per_param = NULL))
  }
  # If SEs exist, do a Wald test; otherwise no detection decision
  if (all(is.finite(del$se))) {
    z   <- abs(del$estimate / del$se)
    cut <- stats::qnorm(1 - alpha/2)
    sig <- z >= cut
    return(list(any = any(sig), n_sig = sum(sig), per_param = setNames(sig, del$name)))
  } else {
    return(list(any = NA, n_sig = NA_integer_, per_param = NULL))
  }
}

