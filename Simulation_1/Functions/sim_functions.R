# ---- sim_functions.R ----
# Purpose: data generator, model builders, fitters, and per-condition simulation

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(mxsem)
  library(lavaan)   
  library(semtree)  
})


# ---------- Data generation (measurement noninvariance DGP) ----------

# Log-link for residual variance to ensure positivity
inv_link_var <- function(eta) exp(eta)

# Given loading λ, reliability ρ, and Var(η)=v_eta, derive residual variance θ
# ρ = (λ^2 v_eta)/(λ^2 v_eta + θ)  =>  θ = λ^2 v_eta (1-ρ)/ρ
theta_from_reliability <- function(lambda, rho, v_eta = 1) {
  stopifnot(all(rho > 0 & rho < 1))
  lambda^2 * v_eta * (1 - rho) / rho
}

# Optional: sample per-condition DGP parameters (use in sim_script, not required)
sample_dgp_parameters <- function(
    n_items = 4,
    lambda_range = c(0.4, 0.9),
    reliability_range = c(0.6, 0.95),
    intercept_mean_sd = c(0, 1),      # N(mean, sd^2) for ν0
    slope_range = c(-0.2, 0.2),       # moderation slope range
    sampled_items = c(1, 2),
    what = c("lambda0","reliability","nu0","slopes")
) {
  what <- match.arg(what, several.ok = TRUE)
  out <- list()
  
  if ("lambda0" %in% what) {
    lam <- runif(n_items - 1, lambda_range[1], lambda_range[2])
    out$lambda0 <- c(lam[1:3], 1)   # keep item4 ≈1 as anchor
  }
  if ("reliability" %in% what) {
    out$reliability <- runif(n_items, reliability_range[1], reliability_range[2])
  }
  if ("nu0" %in% what) {
    out$nu0 <- rnorm(n_items, mean = intercept_mean_sd[1], sd = intercept_mean_sd[2])
  }
  if ("slopes" %in% what) {
    out$slope_load_items <- setNames(rep(0, n_items), paste0("i", 1:n_items))
    out$slope_int_items  <- setNames(rep(0, n_items), paste0("i", 1:n_items))
    out$slope_res_items  <- setNames(rep(0, n_items), paste0("i", 1:n_items))
    for (j in sampled_items) {
      out$slope_load_items[[paste0("i", j)]] <- runif(1, slope_range[1], slope_range[2])
      out$slope_int_items [[paste0("i", j)]] <- runif(1, slope_range[1], slope_range[2])
      out$slope_res_items [[paste0("i", j)]] <- runif(1, slope_range[1], slope_range[2])
    }
  }
  out
}

# Generate data with selective moderation (λ, ν, θ), optional latent mean/variance moderation,
# and optional reliability-based θ baselines.
generate_data_mni <- function(
    n,
    moderator,                                   # numeric length n
    mod_items_load = integer(0),
    mod_items_int  = integer(0),
    mod_items_res  = integer(0),
    # linear slopes on the moderator (for quadratic DGP, pass moderator^2 if desired)
    slope_load = 0,
    slope_int  = 0,
    slope_logres = 0,                             # on log-variance scale if used with log link
    # BASE parameters (defaults match old function)
    lambda0 = c(0.7, 0.7, 0.7, 1.0),              # item4 ~ 1 for ID
    nu0     = c(0, 0, 0, 0),
    eta_mean0 = 0,
    eta_var0  = 1,
    # latent moment moderation:
    #   η_mean = η_mean0 + α_mean * m
    #   η_var  = exp( log(η_var0) + α_logvar * m )  (keeps variance > 0)
    alpha_mean   = 0,
    alpha_logvar = 0,
    # residual baselines either via logres0 OR via reliability (if provided)
    logres0    = log(c(0.3, 0.3, 0.3, 0.3)),
    reliability = NULL,                            # length-4; if given, overrides logres0
    # optional per-condition draws produced by sample_dgp_parameters()
    sampled_pars = NULL,
    use_log_link_resid = TRUE
) {
  stopifnot(length(moderator) == n, length(lambda0) == 4)
  
  # Optional overrides from sampling
  if (!is.null(sampled_pars$lambda0)) lambda0 <- sampled_pars$lambda0
  if (!is.null(sampled_pars$nu0))     nu0     <- sampled_pars$nu0
  
  # Reliability-based θ baseline (at moderator=0), if provided
  if (!is.null(reliability)) {
    theta0  <- theta_from_reliability(lambda = lambda0, rho = reliability, v_eta = eta_var0)
    logres0 <- log(theta0)
  }
  
  # Person-specific latent moments
  eta_mean <- eta_mean0 + alpha_mean * moderator
  eta_var  <- exp(log(eta_var0) + alpha_logvar * moderator)
  eta      <- rnorm(n, mean = eta_mean, sd = sqrt(eta_var))
  
  # Person-specific measurement parameters
  nu     <- matrix(rep(nu0,     each = n), nrow = n, ncol = 4)
  lambda <- matrix(rep(lambda0, each = n), nrow = n, ncol = 4)
  logres <- matrix(rep(logres0, each = n), nrow = n, ncol = 4)
  
  if (length(mod_items_load)) {
    lambda[, mod_items_load] <- sweep(lambda[, mod_items_load, drop = FALSE], 1, slope_load * moderator, `+`)
  }
  if (length(mod_items_int)) {
    nu[, mod_items_int] <- sweep(nu[, mod_items_int, drop = FALSE], 1, slope_int * moderator, `+`)
  }
  if (length(mod_items_res)) {
    logres[, mod_items_res] <- sweep(logres[, mod_items_res, drop = FALSE], 1, slope_logres * moderator, `+`)
  }
  
  # Realize items
  Y <- matrix(NA_real_, n, 4)
  for (j in 1:4) {
    sd_j <- sqrt(if (use_log_link_resid) inv_link_var(logres[, j]) else pmax(1e-8, logres[, j]))
    e <- rnorm(n, mean = 0, sd = sd_j)
    Y[, j] <- nu[, j] + lambda[, j] * eta + e
  }
  
  df <- as.data.frame(Y)
  names(df) <- paste0("item", 1:4)
  df$moderator <- moderator
  df
}


# ---------- Model builders (mxsem strings; definition-variable moderation) ----------

# Baseline model without moderation (strict invariance across moderator)
build_null_model <- function() {
  '
    latent1 =~ 1*item4 + lambda_item1*item1 + lambda_item2*item2 + lambda_item3*item3

    item1 ~~ item1
    item2 ~~ item2
    item3 ~~ item3
    item4 ~~ item4

    latent1 ~~ latent1

    item1 ~ nu_item1*1
    item2 ~ nu_item2*1
    item3 ~ nu_item3*1
    item4 ~ nu_item4*1

    latent1 ~ 0*1
  '
}

# Moderated loadings for a subset of items (linear)
build_model_moderated_loadings <- function(mod_items_load) {
  stopifnot(all(mod_items_load %in% 1:4))
  load_terms <- purrr::map_chr(1:4, function(i) {
    if (i %in% mod_items_load) {
      glue("{{lambda_item{i} := lambda_item{i}_0 + lambda_item{i}_1*data.moderator}}*item{i}")
    } else if (i == 4) {
      "1*item4"
    } else {
      glue("lambda_item{i}*item{i}")
    }
  }) |> paste(collapse = " + ")
  
  glue('
    latent1 =~ {load_terms}

    item1 ~~ item1
    item2 ~~ item2
    item3 ~~ item3
    item4 ~~ item4

    latent1 ~~ latent1

    item1 ~ nu_item1*1
    item2 ~ nu_item2*1
    item3 ~ nu_item3*1
    item4 ~ nu_item4*1

    latent1 ~ 0*1
  ')
}

# Moderated intercepts for a subset of items (linear)
build_model_moderated_intercepts <- function(mod_items_int) {
  stopifnot(all(mod_items_int %in% 1:4))
  int_terms <- purrr::map_chr(1:4, function(i) {
    if (i %in% mod_items_int) {
      glue("item{i} ~ {{nu_item{i} := nu_item{i}_0 + nu_item{i}_1*data.moderator}}*1")
    } else {
      glue("item{i} ~ nu_item{i}*1")
    }
  }) |> paste(collapse = "\n    ")
  
  '
    latent1 =~ 1*item4 + lambda_item1*item1 + lambda_item2*item2 + lambda_item3*item3

    item1 ~~ item1
    item2 ~~ item2
    item3 ~~ item3
    item4 ~~ item4

    latent1 ~~ latent1

  ' |>
    paste0("    ", int_terms, "\n\n    latent1 ~ 0*1\n")
}

# Moderated residual variances for a subset of items
# Set log_link=TRUE if you want variance := exp(a + b*m) to match the DGP.
build_model_moderated_residuals <- function(mod_items_res, log_link = TRUE) {
  stopifnot(all(mod_items_res %in% 1:4))
  
  res_terms <- purrr::map_chr(1:4, function(i) {
    if (i %in% mod_items_res) {
      if (log_link) {
        glue("item{i} ~~ {{theta_item{i} := exp(theta_item{i}_0 + theta_item{i}_1*data.moderator)}}*item{i}")
      } else {
        glue("item{i} ~~ {{theta_item{i} := theta_item{i}_0 + theta_item{i}_1*data.moderator}}*item{i}")
      }
    } else {
      glue("item{i} ~~ item{i}")
    }
  }) |> paste(collapse = "\n    ")
  
  glue('
    latent1 =~ 1*item4 + lambda_item1*item1 + lambda_item2*item2 + lambda_item3*item3

    {residual_block}

    latent1 ~~ latent1

    item1 ~ nu_item1*1
    item2 ~ nu_item2*1
    item3 ~ nu_item3*1
    item4 ~ nu_item4*1

    latent1 ~ 0*1
  ', residual_block = res_terms)
}

# ---------- lavaan template for SEM Trees (moderation-free; labels match mxsem) ----------
build_null_model_lavaan <- function() {
  '
  latent1 =~ 1*item4 + lambda_item1*item1 + lambda_item2*item2 + lambda_item3*item3

  # residual variances labeled so we can "focus" on them if desired
  item1 ~~ th1*item1
  item2 ~~ th2*item2
  item3 ~~ th3*item3
  item4 ~~ th4*item4

  # identification and means
  latent1 ~~ 1*latent1
  latent1 ~ 0*1

  # intercepts (lavaan expects 1*label)
  item1 ~ 1*nu_item1
  item2 ~ 1*nu_item2
  item3 ~ 1*nu_item3
  item4 ~ 1*nu_item4
  '
}

# ---------- internal: weighted slope of leaf estimates on leaf-mean moderator ----------
.wls_slope <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  if (length(x) <= 1) return(list(beta = NA_real_, se = NA_real_, p = NA_real_))
  xbar <- sum(w * x) / sum(w)
  Sxx  <- sum(w * (x - xbar)^2)
  if (Sxx <= 0) return(list(beta = NA_real_, se = NA_real_, p = NA_real_))
  beta <- sum(w * (x - xbar) * y) / Sxx
  resid  <- y - beta * x
  sigma2 <- sum(w * resid^2) / (sum(w) - 2)
  se <- sqrt(sigma2 / Sxx)
  t  <- beta / se
  p  <- 2 * stats::pt(abs(t), df = max(1, sum(w) - 2), lower.tail = FALSE)
  list(beta = beta, se = se, p = p)
}

# ---------- SEM-tree estimator with MNLFA-like outputs ----------
# target_param_label: e.g., "lambda_item1", "nu_item1", or "th1"
estimate_semtree <- function(
    data,
    target_param_label,
    predictors = "moderator",
    focus_parameters = NULL,                        # optional: restrict score tests to these labels
    control = NULL,
    lavaan_syntax = build_null_model_lavaan()
) {
  if (is.null(control)) {
    control <- semtree::semtree.control(
      method    = "score",                        # score-based SEM trees (Arnold et al.)
      min.N     = 100,
      max.depth = 3,
      alpha     = 0.05,
      linear    = TRUE
    )
  }
  
  # 1) fit moderation-free template
  t0 <- Sys.time()
  fit0 <- lavaan::cfa(lavaan_syntax, data = data, std.lv = TRUE)
  
  # 2) optional focus on specific parameter labels
  constr <- if (length(focus_parameters)) semtree::semtree.constraints(focus.parameters = focus_parameters) else NULL
  
  # 3) grow tree
  tree <- semtree::semtree(model = fit0, data = data, predictors = predictors,
                           control = control, constraints = constr)
  rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # 4) leaves, used predictors
  leaves   <- semtree::getLeafs(tree)
  nL       <- length(leaves)
  used     <- unique(stats::na.omit(tree$tree$splitvar)); used <- used[used != ""]
  any_split <- nL > 1L
  
  # 5) per-leaf estimates of the target parameter label
  par_mat <- try(semtree::parameters(tree), silent = TRUE)
  rix <- if (!inherits(par_mat, "try-error")) match(target_param_label, rownames(par_mat)) else NA_integer_
  leaf_est <- if (!is.na(rix)) as.numeric(par_mat[rix, , drop = TRUE]) else rep(NA_real_, nL)
  
  # 6) leaf mean moderator and sizes
  leaf_sizes <- numeric(nL); leaf_mbar <- numeric(nL)
  if (nL > 0) {
    for (i in seq_len(nL)) {
      idx <- leaves[[i]]$case
      leaf_sizes[i] <- length(idx)
      leaf_mbar[i]  <- mean(data$moderator[idx], na.rm = TRUE)
    }
  }
  
  # 7) single slope analog (so your evaluation metrics apply directly)
  slope <- if (nL >= 2 && any(is.finite(leaf_est))) .wls_slope(leaf_mbar, leaf_est, leaf_sizes)
  else list(beta = NA_real_, se = NA_real_, p = NA_real_)
  
  list(
    estimate        = slope$beta,
    se              = slope$se,
    p               = slope$p,
    any_split       = any_split,
    split_on_pred   = any(predictors %in% used),
    n_splits        = max(0L, nL - 1L),
    runtime_seconds = rt
  )
}



# ---------- Fitting and extraction ----------
estimate_mxsem <- function(data, model_string) {
  fit <- mxsem(
    model = model_string,
    data  = data,
    scale_loadings = FALSE,
    scale_latent_variances = FALSE
  )
  res <- mxRun(fit)
  s <- summary(res)
  
  pars <- s$parameters
  z2   <- (pars$Estimate / pars$Std.Error)^2
  pvals <- pchisq(z2, df = 1, lower.tail = FALSE)
  
  est_tab <- tibble::tibble(
    name     = pars$name,
    estimate = pars$Estimate,
    se       = pars$Std.Error,
    p        = pvals
  )
  
  # --- Fit indices via OpenMx reference models (CFI/RMSEA) ---
  ref <- try(OpenMx::mxRefModels(res, run = TRUE), silent = TRUE)
  cfi <- rmsea <- NA_real_
  if (!inherits(ref, "try-error")) {
    chi_model <- max(0, s$Minus2LogLikelihood - ref$Saturated$fitfunction$result)
    df_model  <- max(0, s$degreesOfFreedom      - ref$Saturated$output$degreesOfFreedom)
    chi_indep <- max(0, ref$Independence$fitfunction$result - ref$Saturated$fitfunction$result)
    df_indep  <- max(0, ref$Independence$output$degreesOfFreedom - ref$Saturated$output$degreesOfFreedom)
    
    num <- max(chi_indep - df_indep, 0) - max(chi_model - df_model, 0)
    den <- max(chi_indep - df_indep, 0)
    cfi <- if (den > 0) num / den else NA_real_
    
    N <- nrow(data)
    rmsea <- if (df_model > 0) sqrt(max((chi_model - df_model), 0) / (df_model * N)) else NA_real_
  }
  
  # --- Robust SRMR from standardized covariance residuals ---
  # Uses observed variances to standardize, excludes diagonals, handles NA, and guards zeros.
  srmr <- NA_real_
  srmr_try <- try({
    # model-implied covariance for the manifest block
    Sigma <- OpenMx::mxGetExpected(res, "covariance")
    vars  <- colnames(Sigma)
    # restrict to variables present in the data (exclude moderators/aux)
    vars  <- vars[vars %in% colnames(data)]
    if (length(vars) >= 2) {
      X <- data[, vars, drop = FALSE]
      # observed covariance (pairwise for robustness to missing)
      S <- stats::cov(X, use = "pairwise.complete.obs")
      # reorder S to match Sigma (defensive)
      S <- S[vars, vars, drop = FALSE]
      Sigma <- Sigma[vars, vars, drop = FALSE]
      # standardize by observed SDs
      sd_obs <- sqrt(pmax(diag(S), 1e-12))
      denom  <- outer(sd_obs, sd_obs)
      Rres   <- (S - Sigma) / denom
      # off-diagonals only
      r <- Rres[upper.tri(Rres, diag = FALSE)]
      if (length(r)) srmr <- sqrt(mean(r^2, na.rm = TRUE))
    }
  }, silent = TRUE)
  if (inherits(srmr_try, "try-error")) srmr <- NA_real_
  
  list(
    estimates   = est_tab,
    fit_indices = c(CFI = cfi, RMSEA = rmsea, SRMR = srmr, AIC = s$AIC, BIC = s$BIC),
    mx          = res
  )
}


#----------- pre-helper ----------------------------------
# Choose a single target parameter name for this condition
# Priority: first moderated loading -> intercept -> residual.
choose_target_param_name <- function(mod_items_load, mod_items_int, mod_items_res) {
  if (length(mod_items_load)) return(glue::glue("lambda_item{min(mod_items_load)}_1"))
  if (length(mod_items_int))  return(glue::glue("nu_item{min(mod_items_int)}_1"))
  if (length(mod_items_res))  return(glue::glue("theta_item{min(mod_items_res)}_1"))
  # Fallback: a loading slope that should be 0 if no moderation is present
  "lambda_item1_1"
}

# Map target name to its true value given the condition
true_value_for_target <- function(target_param_name, slope_load, slope_int, slope_logres) {
  if (startsWith(target_param_name, "lambda_item")) return(slope_load)
  if (startsWith(target_param_name, "nu_item"))     return(slope_int)
  if (startsWith(target_param_name, "theta_item"))  return(slope_logres) # on log scale if log-link is used
  0
}
# ---------- Single-condition simulation driver ----------
run_condition <- function(
    n, moderator_vec,
    mod_items_load = integer(0),
    mod_items_int  = integer(0),
    mod_items_res  = integer(0),
    slope_load = 0, slope_int = 0, slope_logres = 0,
    n_reps = 100,
    null_model_string = build_null_model(),
    alt_model_string,
    log_link_resid_in_analysis = TRUE,
    moderator_type = NA_character_   # <-- NEW ARG
) {
  
  # Select one target parameter for full metrics
  target_param_name <- choose_target_param_name(mod_items_load, mod_items_int, mod_items_res)
  target_true_value <- true_value_for_target(target_param_name, slope_load, slope_int, slope_logres)
  
  # true nonzero slope name sets for family-level detection (unchanged)
  true_loading_slopes   <- paste0("lambda_item", mod_items_load, "_1")
  true_intercept_slopes <- paste0("nu_item",     mod_items_int,  "_1")
  true_residual_slopes  <- paste0("theta_item",  mod_items_res,  "_1")
  
  detect_family <- function(est_tab, targets) {
    if (!length(targets) || is.null(est_tab) || !nrow(est_tab)) return(FALSE)
    any(est_tab$name %in% targets & est_tab$p < 0.05, na.rm = TRUE)
  }
  
  # Storage
  det_load <- logical(n_reps)
  det_int  <- logical(n_reps)
  det_res  <- logical(n_reps)
  det_any  <- logical(n_reps)
  fit_null_idx <- matrix(NA_real_, n_reps, 5, dimnames = list(NULL, c("CFI","RMSEA","SRMR","AIC","BIC")))
  fit_alt_idx  <- fit_null_idx
  
  rep_records <- vector("list", n_reps)
  
  for (r in seq_len(n_reps)) {
    # ----- generate data for this replication
    dat <- generate_data_mni(
      n = n,
      moderator = moderator_vec,
      mod_items_load = mod_items_load,
      mod_items_int  = mod_items_int,
      mod_items_res  = mod_items_res,
      slope_load = slope_load,
      slope_int  = slope_int,
      slope_logres = slope_logres
    )
    
    # ----- fit null
    t0 <- Sys.time()
    fit_null <- try(estimate_mxsem(dat, null_model_string), silent = TRUE)
    rt_null <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    conv_null <- !inherits(fit_null, "try-error")
    
    # ----- fit alternative
    t1 <- Sys.time()
    fit_alt <- try(estimate_mxsem(dat, alt_model_string), silent = TRUE)
    rt_alt <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    conv_alt <- !inherits(fit_alt, "try-error")
    
    # extract indices
    if (conv_null) fit_null_idx[r, ] <- fit_null$fit_indices
    if (conv_alt)  fit_alt_idx[r, ]  <- fit_alt$fit_indices
    
    # detection flags (family-level)
    est_alt <- if (conv_alt) fit_alt$estimates else NULL
    det_load[r] <- detect_family(est_alt, true_loading_slopes)
    det_int[r]  <- detect_family(est_alt, true_intercept_slopes)
    det_res[r]  <- detect_family(est_alt, true_residual_slopes)
    det_any[r]  <- any(det_load[r], det_int[r], det_res[r])
    
    # ----- per-replication record for the single target parameter (from alternative fit)
    row <- if (conv_alt && !is.null(est_alt) && nrow(est_alt)) {
      est_alt[match(target_param_name, est_alt$name), , drop = FALSE]
    } else NULL
    
    estimate <- if (!is.null(row) && nrow(row) == 1) row$estimate else NA_real_
    se       <- if (!is.null(row) && nrow(row) == 1) row$se       else NA_real_
    pval     <- if (!is.null(row) && nrow(row) == 1) row$p        else NA_real_
    
    # alt-fit row
    rep_records[[r]] <- tibble::tibble(
      rep = r,
      method = "mxsem",
      target_param_name = target_param_name,
      estimate = estimate,
      se = se,
      p = pval,
      converged = conv_alt,
      runtime_seconds = rt_alt,
      true_value = target_true_value,
      n = n,
      moderator_type = moderator_type,   
      mod_items_load = list(mod_items_load),
      mod_items_int  = list(mod_items_int),
      mod_items_res  = list(mod_items_res),
      slope_load = slope_load,
      slope_int  = slope_int,
      slope_logres = slope_logres
    )
    
    # ---- APPEND null-fit row *inside* the loop ----
    rep_records[[r]] <- dplyr::bind_rows(
      rep_records[[r]],
      tibble::tibble(
        rep = r,
        method = "mxsem_null",
        target_param_name = NA_character_,
        estimate = NA_real_,
        se = NA_real_,
        p = NA_real_,
        converged = conv_null,
        runtime_seconds = rt_null,
        true_value = NA_real_,
        n = n,
        moderator_type = moderator_type,   
        mod_items_load = list(mod_items_load),
        mod_items_int  = list(mod_items_int),
        mod_items_res  = list(mod_items_res),
        slope_load = slope_load,
        slope_int  = slope_int,
        slope_logres = slope_logres
      )
    )
    # ---- SEM-tree (score-based)----
    # Map mxsem target name (e.g., "lambda_item1_1") to lavaan label ("lambda_item1")
    target_label <- sub("_1$", "", target_param_name)
    
    # Focus set = the exact families/items moderated in the DGP
    focus_params <- unique(c(
      paste0("lambda_item", mod_items_load),
      paste0("nu_item",     mod_items_int),
      paste0("th",          mod_items_res)   # residuals are th1..th4 in the lavaan template
    ))
    
    st <- try(estimate_semtree(
      data = dat,
      target_param_label = target_label,
      predictors = "moderator",
      focus_parameters = focus_params
    ), silent = TRUE)
    
    if (inherits(st, "try-error")) {
      rep_records[[r]] <- dplyr::bind_rows(
        rep_records[[r]],
        tibble::tibble(
          rep = r,
          method = "semtree",
          target_param_name = target_label,
          estimate = NA_real_,
          se = NA_real_,
          p = NA_real_,
          converged = NA,                         # not applicable to trees
          runtime_seconds = NA_real_,
          true_value = target_true_value,
          n = n,
          moderator_type = moderator_type,
          mod_items_load = list(mod_items_load),
          mod_items_int  = list(mod_items_int),
          mod_items_res  = list(mod_items_res),
          slope_load = slope_load,
          slope_int  = slope_int,
          slope_logres = slope_logres,
          any_split = NA,
          split_on_moderator = NA,
          n_splits = NA_integer_
        )
      )
    } else {
      rep_records[[r]] <- dplyr::bind_rows(
        rep_records[[r]],
        tibble::tibble(
          rep = r,
          method = "semtree",
          target_param_name = target_label,
          estimate = st$estimate,                 # tree-based slope analog
          se = st$se,
          p = st$p,
          converged = NA,                         # not applicable
          runtime_seconds = st$runtime_seconds,
          true_value = target_true_value,
          n = n,
          moderator_type = moderator_type,
          mod_items_load = list(mod_items_load),
          mod_items_int  = list(mod_items_int),
          mod_items_res  = list(mod_items_res),
          slope_load = slope_load,
          slope_int  = slope_int,
          slope_logres = slope_logres,
          any_split = st$any_split,
          split_on_moderator = st$split_on_pred,
          n_splits = st$n_splits
        )
      )
    }
    
  }  # end for-loop
  
  # ----- condition-level summary 
  summary_tbl <- tibble::tibble(
    n = n,
    mod_items_load = paste(mod_items_load, collapse = ","),
    mod_items_int  = paste(mod_items_int,  collapse = ","),
    mod_items_res  = paste(mod_items_res,  collapse = ","),
    slope_load = slope_load,
    slope_int = slope_int,
    slope_logres = slope_logres,
    reps = n_reps,
    det_rate_load = mean(det_load),
    det_rate_int  = mean(det_int),
    det_rate_res  = mean(det_res),
    det_rate_any  = mean(det_any),
    AIC_null = mean(fit_null_idx[, "AIC"], na.rm = TRUE),
    AIC_alt  = mean(fit_alt_idx[,  "AIC"], na.rm = TRUE),
    BIC_null = mean(fit_null_idx[, "BIC"], na.rm = TRUE),
    BIC_alt  = mean(fit_alt_idx[,  "BIC"], na.rm = TRUE),
    moderator_type = moderator_type
  )
  
  # ----- add SEM-tree summaries for this condition -----
  # bind per-rep rows once
  reps_tbl <- dplyr::bind_rows(rep_records)
  
  # subset to SEM-tree rows (may be zero if tree failed every time)
  semt_rows <- dplyr::filter(reps_tbl, .data$method == "semtree")
  
  # helper to compute a mean with NA handling that returns NA if no rows
  .mean_or_na <- function(x) if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE)
  
  summary_tbl <- summary_tbl %>%
    dplyr::mutate(
      det_rate_tree_any     = .mean_or_na(semt_rows$any_split),          # proportion of reps with ≥1 split
      det_rate_tree_on_mod  = .mean_or_na(semt_rows$split_on_moderator), # proportion of reps where moderator used
      mean_tree_splits      = .mean_or_na(semt_rows$n_splits)            # avg number of splits
    )
  
  # ----- return -----
  list(
    summary = summary_tbl,
    rep_records = reps_tbl
  )
}
