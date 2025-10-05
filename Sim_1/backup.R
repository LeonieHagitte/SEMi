run_one <- function(row) {
  # 1) Case-wise data (contains y1..y4, M, M2)
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
  dat <- sim$data
  truths <- sim$true_deltas  # named δ’s (zeros under NULL/partial as appropriate)
  
  # 2) Build analysis syntaxes
  mnlfa_model <- build_mnlfa_linear_full()  # all δ’s free; aligns with truths’ names
  cfa_base    <- build_cfa_baseline()       # unmoderated CFA for SEM-Tree
  
  # 3) Fit methods
  fit_mnlfa <- run_analysis(dat, mnlfa_model, method = "MNLFA")  # mxsem backend
  fit_tree  <- run_analysis(dat, cfa_base,    method = "SEMTREE")       # semtree backend
  
  # 4) Metrics
  #    δ-metrics for MNLFA (labels must match: delta_lambda_y1..y4, delta_theta_y1..y4)
  out <- list(
    model_type  = row$model_type, N = row$N, rep = row$rep
  )
  if (inherits(fit_mnlfa, "MxModel") || inherits(fit_mnlfa, "MxRAMModel")) {
    out$delta_cov   <- calculate_coverage_deltas(fit_mnlfa, truths)
    out$delta_bias  <- calculate_relative_bias_deltas(fit_mnlfa, truths)
    out$delta_rmse  <- calculate_relative_rmse_deltas(fit_mnlfa, truths)
  } else {
    out$delta_cov <- out$delta_bias <- out$delta_rmse <- NA_real_
  }
  
  #    SEM-Tree moderation detection proxy (did it split on M? how often?)
  det <- semtree_detects_moderation(fit_tree, moderator_name = "M")
  out$tree_split_on_M   <- det$split_on_moderator
  out$tree_n_splits_M   <- det$n_splits_on_moderator
  
  out
}



mxsem_param_table <- function(fit) {
  s <- OpenMx::summary(fit)  # does not require library(OpenMx); namespaced call
  par <- tibble::as_tibble(s$parameters)   # columns include 'name', 'Estimate', 'Std.Error', ...
  ci  <- tryCatch(tibble::as_tibble(s$CI), error = function(e) NULL)
  
  if (!is.null(ci)) {
    ci <- dplyr::rename(ci,
                        ci.lower = dplyr::coalesce(.data$lbound, .data$lower, .data$Lower),
                        ci.upper = dplyr::coalesce(.data$ubound, .data$upper, .data$Upper)
    )
    ci <- dplyr::select(ci, name, ci.lower, ci.upper)
    par <- dplyr::left_join(par, ci, by = "name")
  } else {
    par$ci.lower <- NA_real_; par$ci.upper <- NA_real_
  }
  if (!"Estimate" %in% names(par)) stop("mxsem_param_table: 'Estimate' not found.")
  if (!"Std.Error" %in% names(par)) par$Std.Error <- NA_real_
  
  dplyr::rename(par, estimate = .data$Estimate, se = .data$Std.Error)
}


if (method == "MNLFA") {
  mxmodel <- mxsem::mxsem(model_syntax, data = dat)
  fit <- tryCatch(OpenMx::mxRun(mxmodel, intervals = mx_intervals), error = identity)
  return(fit)
}


run_analysis <- function(data, model_syntax,
                         method = c("MNLFA", "SEMTREE"),
                         tree_predictors = c("M","M2"),
                         tree_control = semtree::semtree.control(method = "score"),
                         mx_intervals = TRUE) {
  method <- match.arg(method)
  dat <- as.data.frame(data)
  
  if (method == "MNLFA") {
    mxmodel <- mxsem::mxsem(model_syntax, data = dat)
    
    delta_names <- c("delta_lambda_y1","delta_lambda_y2","delta_lambda_y3","delta_lambda_y4")
    mxmodel <- OpenMx::mxModel(mxmodel, OpenMx::mxCI(delta_names))
    
    fit <- OpenMx::mxRun(mxmodel, intervals = TRUE)
    return(fit)
  }
  
}

if (method == "SEMTREE") {
  lav_fit <- lavaan::cfa(model_syntax, data = dat)
  tree <- tryCatch(
    semtree::semtree(model = lav_fit, data = dat,
                     control = tree_control,
                     predictors = tree_predictors),
    error = identity
  )
  return(tree)
}

stop("Unknown method in run_analysis().")
}



run_analysis <- function(data, model_syntax,
                         method = c("MNLFA", "SEMTREE"),
                         # SEM Tree knobs
                         tree_predictors = c("M", "M2"),
                         tree_control    = semtree::semtree.control(method = "score"),
                         # MNLFA knobs
                         mx_intervals    = TRUE,
                         delta_names     = c("delta_lambda_y1",
                                             "delta_lambda_y2",
                                             "delta_lambda_y3",
                                             "delta_lambda_y4")) {
  
  method <- match.arg(method)
  dat <- as.data.frame(data)   # defensive: some APIs prefer plain data.frame
  
  if (method == "MNLFA") {
    # --- MNLFA via mxsem/OpenMx ----------------------------------------------
    # Build an OpenMx model from lavaan-like syntax
    mxmodel <- mxsem::mxsem(model_syntax, data = dat)
    
    # Ask OpenMx to compute likelihood-based CIs for the δ parameters.
    # (These labels must exist in your syntax; if not, mxRun() may warn/error.)
    if (!is.null(delta_names) && length(delta_names) > 0) {
      mxmodel <- OpenMx::mxModel(mxmodel, OpenMx::mxCI(delta_names))
    }
    
    # Fit the model; intervals = TRUE triggers CI computation
    fit <- tryCatch(
      OpenMx::mxRun(mxmodel, intervals = isTRUE(mx_intervals)),
      error = identity
    )
    return(fit)
  } else if (method == "SEMTREE") {
    lav_fit <- lavaan::cfa(model_syntax, data = dat)
    tree <- tryCatch(
      semtree::semtree(model = lav_fit, data = dat,
                       control = tree_control, predictors = tree_predictors),
      error = identity
    )
    
    # Try to unwrap to a 'party' object for downstream use
    if (!inherits(tree, "error")) {
      if (inherits(tree, "party")) {
        return(tree)
      } else if (is.list(tree) && !is.null(tree$tree) && inherits(tree$tree, "party")) {
        return(tree$tree)
      } else if (methods::is(tree, "semtree")) {
        p_try <- tryCatch(tree@tree, error = function(e) NULL)
        if (inherits(p_try, "party")) return(p_try)
      }
    }
    return(tree)  # fall back; helper will handle NA gracefully
  }
  
  
  # If method was neither of the above (shouldn’t happen due to match.arg)
  stop("Unknown method in run_analysis().")
}

# Detect splits on a set of moderators from a semtree result

semtree_detects_moderation <- function(tree, moderator = c("M","M2")) {
  # Return NA columns if tree is an error/NULL
  if (inherits(tree, "error") || is.null(tree)) {
    out <- list()
    for (m in moderator) {
      out[[paste0("tree_split_on_", m)]] <- NA
      out[[paste0("tree_n_splits_", m)]] <- NA_integer_
    }
    return(out)
  }
  
  p <- if (inherits(tree, "party")) tree else partykit::as.party(tree)
  ids <- partykit::nodeids(p, terminal = FALSE)
  
  split_vars <- character(0)
  if (length(ids) > 0) {
    split_vars <- unlist(
      partykit::nodeapply(p, ids, FUN = function(nd) {
        sp <- partykit::split_node(nd)
        if (is.null(sp)) return(NULL)
        names(partykit::data_party(p))[partykit::varid_split(sp)]
      })
    )
  }
  
  out <- list()
  for (m in moderator) {
    k <- sum(split_vars == m, na.rm = TRUE)
    out[[paste0("tree_split_on_", m)]] <- (k > 0)
    out[[paste0("tree_n_splits_", m)]] <- as.integer(k)
  }
  out
}





mxsem_param_table <- function(fit) {
  s <- summary(fit)
  # Parameters table
  par <- tibble::as_tibble(s$parameters)
  ci <- tryCatch(tibble::as_tibble(s$CI), error = function(e) NULL)
  
  if (!is.null(ci)) {
    # Standardize names to lower-case to handle lbound/LBound/lower variants
    names(ci) <- tolower(names(ci))
    
    # Some OpenMx versions use 'parameter' instead of 'name'
    if (!"name" %in% names(ci) && "parameter" %in% names(ci)) {
      ci <- dplyr::rename(ci, name = parameter)
    }
    
    # Create unified ci.lower / ci.upper columns with mutate + coalesce (not rename!)
    ci <- ci %>%
      dplyr::mutate(
        ci.lower = dplyr::coalesce(.data$lbound, .data$lower),
        ci.upper = dplyr::coalesce(.data$ubound, .data$upper)
      ) %>%
      dplyr::select(dplyr::any_of(c("name", "ci.lower", "ci.upper")))
    
    # Join CI bounds back to parameter table
    par <- dplyr::left_join(par, ci, by = "name")
  } else {
    par$ci.lower <- NA_real_
    par$ci.upper <- NA_real_
  }
  
  # Ensure standard column names for downstream code
  if (!"Estimate"  %in% names(par)) stop("mxsem_param_table: Estimate not found.")
  if (!"Std.Error" %in% names(par)) par$Std.Error <- NA_real_
  
  par <- dplyr::rename(par, estimate = .data$Estimate, se = .data$`Std.Error`)
  par
}

mxsem_param_table <- function(fit) {
  # 1) Summary (use base generic; OpenMx registers the method)
  s <- summary(fit)
  
  # 2) Parameters (estimates + SEs)
  par <- tryCatch(tibble::as_tibble(s$parameters), error = function(e) NULL)
  if (is.null(par) || nrow(par) == 0) {
    stop("mxsem_param_table: No 'parameters' found in summary(fit).")
  }
  
  # Normalize estimate / se column names
  nm <- names(par)
  est_col <- intersect(nm, c("Estimate","estimate","Est","est"))
  se_col  <- intersect(nm, c("Std.Error","Std.Error.","SE","se","StdError"))
  if (length(est_col) == 0) stop("mxsem_param_table: Could not find Estimate column.")
  if (length(se_col)  == 0) se_col <- NA_character_  # OK: SE may be missing
  
  par <- par %>%
    dplyr::mutate(
      name     = as.character(.data$name),
      estimate = .data[[est_col[1]]],
      se       = if (!is.na(se_col[1])) .data[[se_col[1]]] else NA_real_
    ) %>%
    dplyr::select(name, estimate, se)
  
  # 3) Confidence intervals (may be absent even with intervals=TRUE)
  ci_raw <- tryCatch(tibble::as_tibble(s$CI), error = function(e) NULL)
  if (is.null(ci_raw) || nrow(ci_raw) == 0) {
    ci_raw <- tryCatch(tibble::as_tibble(fit$output$confidenceIntervals),
                       error = function(e) NULL)
  }
  
  if (!is.null(ci_raw) && nrow(ci_raw) > 0) {
    cn <- names(ci_raw)
    name_col  <- intersect(cn, c("name","parameter","param","Parameter","labels"))
    lower_col <- intersect(cn, c("ci.lower","lower","lbound","Lower","lBound"))
    upper_col <- intersect(cn, c("ci.upper","upper","ubound","Upper","uBound"))
    
    if (length(name_col) > 0 && length(lower_col) > 0 && length(upper_col) > 0) {
      ci <- ci_raw %>%
        dplyr::transmute(
          name     = as.character(.data[[name_col[1]]]),
          ci.lower = suppressWarnings(as.numeric(.data[[lower_col[1]]])),
          ci.upper = suppressWarnings(as.numeric(.data[[upper_col[1]]]))
        )
      tab <- dplyr::left_join(par, ci, by = "name")
    } else {
      tab <- par %>% dplyr::mutate(ci.lower = NA_real_, ci.upper = NA_real_)
    }
  } else {
    tab <- par %>% dplyr::mutate(ci.lower = NA_real_, ci.upper = NA_real_)
  }
  
  tab
}

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
  truths <- sim$true_deltas                   # named vector of true δ's
  
  # 2) Analysis syntaxes
  mnlfa_model <- build_mnlfa_linear_full()    # make sure δ labels match eval 
  cfa_base    <- build_cfa_baseline()
  
  # 3) Fit methods
  
  # TODO: AB->LH really hardcoded?
  delta_names <- paste0("delta_lambda_y",1:3)
  
  fit_mnlfa <- run_analysis(
    dat,
    mnlfa_model,
    method       = "MNLFA",
    mx_intervals = TRUE,
    delta_names  = delta_names   # <- asks for CIs for every true δ 
  )
  # <- use exact method key
  fit_tree  <- run_analysis(dat, cfa_base,    method = "SEMTREE")
  
  if (inherits(fit_mnlfa, "simpleError")) {
    print(fit_mnlfa)
    #stop("Error in MNLFA")
  }
  
  # 4) δ metrics (mxsem-aware). Requires delta_metrics_any() from eval_functions.R
  dm <- if (inherits(fit_mnlfa, "MxModel") || inherits(fit_mnlfa, "MxRAMModel")) {
    # If δ labels are not prefixed "delta_", set label_prefix accordingly
    delta_metrics_any(fit_mnlfa, truths, label_prefix = "delta_")
  } else {
    list(bias = NA_real_, rmse_abs = NA_real_, coverage = NA_real_)
  }
  
  # MNLFA detection flag (any δ significant?)
  mn_det <- detect_moderation_mxsem(fit_mnlfa, label_prefix = "delta_", alpha = 0.05)
  mnlfa_detect <- mn_det$any
  mnlfa_n_sig  <- mn_det$n_sig
  
  # SEM-Tree split detection
  tree_count <- list()
  predictors <- getPredictorsFromTree(fit_tree)
  moderator = c("M","M2")
  for (m in moderator) {
    if (!is.null(predictors)) {
      k <- sum(predictors==m)
    } else {
      k <- 0
    }
    tree_count[[paste0("tree_n_splits_", m)]] <- k
    tree_count[[paste0("tree_split_on_", m)]] <- (k > 0)
  }
  
  
  
  # 5) One tidy row
  tibble::tibble(
    model_type        = row$model_type,
    N                 = row$N,
    reliability       = row$reliability,
    lambda            = row$lambda,
    moderator_1_type  = row$moderator_1_type,
    
    mnlfa_detect      = mnlfa_detect,
    mnlfa_n_sig       = mnlfa_n_sig,
    
    tree_split_on_M   = tree_count$tree_split_on_M,
    tree_n_splits_M   = tree_count$tree_n_splits_M,
    tree_split_on_M2  = tree_count$tree_split_on_M2,
    tree_n_splits_M2  = tree_count$tree_n_splits_M2,
    
    delta_lambda_full = row$delta_lambda_full,
    delta_theta_full  = row$delta_theta_full,
    delta_lambda_12   = row$delta_lambda_12,
    rep               = row$rep,
    
    delta_bias        = dm$bias,
    delta_rmse        = dm$rmse_abs,      # absolute RMSE for δ's
    delta_coverage    = dm$coverage      # NA if OpenMx CIs weren’t computed
    
  )
}


mxsem_param_table <- function(fit_mnlfa) {
  
  s <- summary(fit_mnlfa)
  par <- s$parameters
  
  ci  <- tryCatch(tibble::as_tibble(s$CI), error = function(e) NULL)
  ci$name <- rownames(s$CI)
  
  if (!is.null(ci)) {
    ci <- ci %>%
      rename(ci.lower = any_of(c("lbound", "lower", "Lower"))) %>%
      rename(ci.upper = any_of(c("ubound", "upper", "Upper")))
    
    ci <- dplyr::select(ci, name, ci.lower, ci.upper)
    par <- dplyr::left_join(par, ci, by = "name")
  } else {
    par$ci.lower <- NA_real_
    par$ci.upper <- NA_real_
  }
  
  if (!"Estimate"  %in% names(par)) stop("mxsem_param_table: Estimate not found.")
  if (!"Std.Error" %in% names(par)) par$Std.Error <- NA_real_
  
  dplyr::rename(par, estimate = Estimate, se = Std.Error)
}