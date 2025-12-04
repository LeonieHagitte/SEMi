# analysis.R
# -----------------------------------------------------------------------------
# Purpose
#   Provide analysis-model builders used by your driver:
#     - build_cfa_baseline(): unmoderated CFA for SEM Tree partitioning
#     - build_mnlfa_syntax(): MNLFA-style model with loading moderation via
#       a latent-by-observed interaction (labels "delta_lambda_y?").
#
# Notes (important):
#   1) In lavaan, *loading moderation by an observed M* can be implemented by
#      creating a latent version of M (no measurement error), then defining a
#      latent interaction:  int_f1M =~ f1 XWITH m
#      and regressing indicators on f1 and the interaction term.
#   2) *Residual variance moderation* (Θ moderated by M) is not natively
#      supported in lavaan. Options: (i) detect via SEM Trees (they will split),
#      (ii) approximate with grouping, or (iii) use OpenMx/Mplus for true
#      variance regression. Below we implement *loading moderation* directly.
# -----------------------------------------------------------------------------
# Load necessary libraries
library(lavaan)
library(dplyr)
library(purrr)
library(mxsem)
library(OpenMx)
library(semtree)
library(lavaan)
library(partykit)

# -----------------------------------------------------------------------------
# Baseline (for SEM Tree) and MNLFA (mxsem) builders for a 1-factor, 4-item model
# -----------------------------------------------------------------------------

# ---------------------- unmoderated CFA (SEM Tree base) ----------------------
build_cfa_baseline <- function() {
  "
  # Measurement: marker loading on y4
  f1 =~ y1 + y2 + y3 + 1*y4

  # Residual variances (Θ): free
  y1 ~~ y1
  y2 ~~ y2
  y3 ~~ y3
  y4 ~~ y4

  # Latent variance (Ψ): free; latent mean fixed for ID
  f1 ~~ f1
  f1 ~ 0*1

  # Intercepts (means) of observed variables: free
  y1 ~ 1
  y2 ~ 1
  y3 ~ 1
  y4 ~ 1
  "
}

# ---------------------- MNLFA via mxsem (linear link) ------------------------
# Full moderation (loadings + residual variances), linear in M.
# Labels 'delta_lambda_y*' and 'delta_theta_y*' need to match δ-metrics.
build_mnlfa_linear_full <- function() {
  "
  # Loadings: λ_j(M) = λ0_j + δ_j * M
  f1 =~ {l1 := l10 + delta_lambda_y1*data.M}*y1 +
        {l2 := l20 + delta_lambda_y2*data.M}*y2 +
        {l3 := l30 + delta_lambda_y3*data.M}*y3 +
        1*y4

  # Residual variances (log link): θ_j(M) = exp( θ0_j + δθ_j * M )
  y1 ~~ {v1 := exp(v10 + delta_theta_y1*data.M)}*y1
  y2 ~~ {v2 := exp(v20 + delta_theta_y2*data.M)}*y2
  y3 ~~ {v3 := exp(v30 + delta_theta_y3*data.M)}*y3
  y4 ~~ {v4 := exp(v40 + delta_theta_y4*data.M)}*y4

  # Latent variance free; latent mean fixed for ID
  f1 ~~ f1
  f1 ~ 0*1

  # Intercepts free
  y1 ~ 1
  y2 ~ 1
  y3 ~ 1
  y4 ~ 1
  "
}

# Partial moderation of loadings only for y1,y2 (no residual moderation).
build_mnlfa_linear_partial <- function() {
  "
  f1 =~ {l1 := l10 + delta_lambda_y1*data.M}*y1 +
        {l2 := l20 + delta_lambda_y2*data.M}*y2 +
        l30*y3 + 1*y4

  # Residual variances unmoderated (constant, positive)
  y1 ~~ {v1 := exp(v10)}*y1
  y2 ~~ {v2 := exp(v20)}*y2
  y3 ~~ {v3 := exp(v30)}*y3
  y4 ~~ {v4 := exp(v40)}*y4

  f1 ~~ f1
  f1 ~ 0*1

  y1 ~ 1
  y2 ~ 1
  y3 ~ 1
  y4 ~ 1
  "
}
# ---------------------------------------------------
# ----------------------------------------------------
run_analysis <- function(data, model_syntax,
                         method = c("MNLFA", "SEMTREE"),
                         # SEM Tree knobs
                         tree_predictors = c("M", "M2"),
                         tree_control    = semtree::semtree.control(method = "score"),
                         # MNLFA knobs
                         mx_intervals    = TRUE,
                         delta_names     = c("delta_lambda_y1",
                                             "delta_lambda_y2",
                                             "delta_lambda_y3"
                         )) {
  
  method <- match.arg(method)
  dat <- as.data.frame(data)
  
  if (method == "MNLFA") {
    mxmodel <- mxsem::mxsem(model_syntax, data = dat)
    if (!is.null(delta_names) && length(delta_names) > 0) {
      mxmodel <- OpenMx::mxModel(mxmodel, OpenMx::mxCI(delta_names))
    }
    fit <- tryCatch(OpenMx::mxRun(mxmodel, intervals = isTRUE(mx_intervals)),
                    error = identity)
    return(fit)
  }
  
  if (method == "SEMTREE") {
    # Fit unmoderated CFA in lavaan (baseline)
    lav_fit <- lavaan::cfa(model_syntax, data = dat)
    # Build constraints: only evaluate splits w.r.t. these loadings
    cnst <- semtree::semtree.constraints(
      focus.parameters = c("f1__y1","f1__y2","f1__y3","f1__y4")
    )
    # Grow SEM tree (score-based splitting)
    st <- tryCatch(
      semtree::semtree(
        model      = lav_fit,
        data       = dat,
        control    = tree_control,
        constraints = cnst,
        predictors = tree_predictors
      ),
      error = identity
    )
    return(st)  # either an S4 'semtree' or an 'error'
  }
  
  stop("Unknown method in run_analysis().")
}

# -----------------------
# helpers-semtree.R
# Expectation: 'st' is S4 object of class "semtree" with slot @tree 
semtree_detects_moderation <- function(st, moderator = c("M","M2")) {
  # Prebuild NA outputs (one boolean + one count per moderator)
  out <- list()
  for (m in moderator) {
    out[[paste0("tree_split_on_", m)]] <- NA
    out[[paste0("tree_n_splits_", m)]] <- NA_integer_
  }
  
  # If the build failed or is NULL, keep NA
  if (inherits(st, "error") || is.null(st)) return(out)
  
  # Enforce the semtree 
  if (!methods::is(st, "semtree")) return(out)
  
  # Extract the internal 'party' tree 
  pt <- tryCatch(methods::slot(st, "tree"), error = function(e) NULL)
  if (!inherits(pt, "party")) return(out)
  
  # Traverse internal (non-terminal) nodes, collect split variable names
  ids <- partykit::nodeids(pt, terminal = FALSE)
  split_vars <- character(0)
  if (length(ids) > 0) {
    split_vars <- unlist(partykit::nodeapply(pt, ids, FUN = function(nd) {
      sp <- partykit::split_node(nd)
      if (is.null(sp)) return(NULL)
      names(partykit::data_party(pt))[partykit::varid_split(sp)]
    }))
  }
  
  # Tally per moderator
  for (m in moderator) {
    k <- sum(split_vars == m, na.rm = TRUE)
    out[[paste0("tree_split_on_", m)]] <- (k > 0)
    out[[paste0("tree_n_splits_", m)]] <- as.integer(k)
  }
  out
}


#getPredictorsFromTree <- function(tree) {
#  if (tree$caption == "TERMINAL") {
#    return()
#  } else {
#    l <- getPredictorsFromTree(tree$left_child)
#    r <- getPredictorsFromTree(tree$right_child)
#    return(c(l, r, tree$rule$name))
#  }
#}

getPredictorsFromTree <- function(st) {
  if (inherits(st, "error") || is.null(st)) return(NULL)
  if (!methods::is(st, "semtree")) return(NULL)
  
  pt <- tryCatch(methods::slot(st, "tree"), error = function(e) NULL)
  if (!inherits(pt, "party")) return(NULL)
  
  ids <- partykit::nodeids(pt, terminal = FALSE)
  if (length(ids) == 0L) return(NULL)
  
  split_vars <- unlist(partykit::nodeapply(pt, ids, FUN = function(nd) {
    sp <- partykit::split_node(nd)
    if (is.null(sp)) return(NULL)
    names(partykit::data_party(pt))[partykit::varid_split(sp)]
  }))
  
  split_vars
}
