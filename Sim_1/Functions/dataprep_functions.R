# dataprep_functions.R
# -----------------------------------------------------------------------------
# Purpose
#   Generate SEM population matrices for a single latent factor with four
#   indicators and explicit moderator handling, in the coding style of gen_mat.R.
#   Supported model_type:
#     - "NULL" : null model (no moderation), both moderators present and non-informative
#     - "1.1"  : full moderation of ALL loadings and ALL residual variances (THETA)
#     - "1.2"  : partial moderation of loadings for items 1 and 2 only (THETA fixed)
#
#   Always returns two moderators in the result:
#     - moderator_1: can be linear / sigmoid / quadratic / noise
#     - moderator_2: always non-informative (noise; has no effect on parameters)
#
# Matrices returned:
#   - lambda (LAMBDA) : 4 x 1 loadings; reference loading is item4 = 1
#   - theta  (THETA)  : 4 x 4 diagonal of residual variances
#   - psi    (PSI)    : 1 x 1 latent variance (Var(latent1))
#   - beta   (BETA)   : 1 x 1 structural matrix (always 0 for single factor)
# -----------------------------------------------------------------------------

# Load necessary libraries
library(lavaan)
library(Matrix)

# --- internal: transform h(M) for moderator_1
.mod_h <- function(M, type = c("linear","sigmoid","quadratic","noise"), slope = 2.0) {
  type <- match.arg(type)
  switch(type,
         linear    = M,
         sigmoid   = 2 * plogis(slope * M) - 1,   # R → (-1,1)
         quadratic = M^2,
         noise     = 0)
}

# Generate matrices
gen_mat <- function(model_type,
                    # fixed single-factor, four indicators:
                    nfactors = 1, nvar.factor = 4,
                    # base parameters
                    lambda = 0.70,            # base loading for items 1–3; item4 is reference = 1
                    reliability = 0.80,       # per-indicator reliability target at base
                    psi_var = 1.0,            # Var(latent1)
                    # moderation controls (only moderator_1 can affect parameters)
                    moderator_1_value = 0,
                    moderator_1_type  = c("linear","sigmoid","quadratic","noise"),
                    sigmoid_slope_1   = 2.0,
                    # effect magnitudes
                    delta_lambda_full = 0.20, # for model 1.1, adds to ALL loadings: λ_i(M)=λ_i0 + delta* h(M)
                    delta_theta_full  = 0.20, # for model 1.1, scales THETA diag multiplicatively: θ_i(M)=θ_i0*(1+delta*h)
                    delta_lambda_12   = 0.20, # for model 1.2, adds to items 1 & 2 only
                    # always-present non-informative moderator_2
                    moderator_2_value = 0) {
  
  stopifnot(nfactors == 1, nvar.factor == 4)
  
  moderator_1_type <- match.arg(moderator_1_type)
  
  # ---------------- Base (unmoderated) model ----------------
  # Order matches your lavaan syntax: latent1 =~ item4 + item1 + item2 + item3
  # We keep matrix rows in item order (item1..item4), but set the reference on item4.
  # LAMBDA_base: [item1,item2,item3,item4]^T
  LAMBDA_base <- matrix(c(lambda, lambda, lambda, 1.0), nrow = 4, ncol = 1)
  
  # Latent variance
  PSI <- matrix(psi_var, nrow = 1, ncol = 1)
  
  # Base THETA from reliability at the base (unmoderated) model:
  # Var_true_i = λ_i^2 * Var(η); θ_ii = Var_true_i / rel - Var_true_i
  true_var_base <- as.numeric(LAMBDA_base^2) * as.numeric(PSI)
  theta_diag_base <- true_var_base / reliability - true_var_base
  stopifnot(all(theta_diag_base > 0))
  THETA_base <- diag(theta_diag_base, nrow = 4, ncol = 4)
  
  # Structural (single factor): no within-latent regressions
  BETA <- matrix(0, nrow = 1, ncol = 1)
  
  # ---------------- Moderators (both are always returned) ----------------
  # moderator_1 may be informative; moderator_2 is fixed noise (non-informative).
  hM1 <- .mod_h(moderator_1_value, moderator_1_type, slope = sigmoid_slope_1)
  hM2 <- 0  # by design, moderator_2 is always "noise"
  
  # ---------------- Apply per model_type ----------------
  if (model_type == "NULL") {
    # Null model: no moderation; both moderators present but non-informative.
    LAMBDA <- LAMBDA_base
    THETA  <- THETA_base
    
  } else if (model_type == "1.1") {
    # Full moderation: ALL loadings and ALL residual variances depend on moderator_1.
    # Loadings: additive shift by delta_lambda_full * hM1 (keeps reference free to vary).
    LAMBDA <- LAMBDA_base + delta_lambda_full * hM1
    
    # Residual variances: multiplicative scaling of base θ by (1 + delta_theta_full * hM1)
    # This keeps θ > 0 if |delta*hM1| < 1.
    scale_theta <- (1 + delta_theta_full * hM1)
    if (scale_theta <= 0) stop("Residual variance scaling became non-positive; reduce delta_theta_full or hM1.")
    THETA <- THETA_base * as.numeric(scale_theta)
    
  } else if (model_type == "1.2") {
    # Partial moderation: ONLY loadings for items 1 and 2 are moderated; THETA fixed at base.
    LAMBDA <- LAMBDA_base
    LAMBDA[1, 1] <- LAMBDA[1, 1] + delta_lambda_12 * hM1  # item1
    LAMBDA[2, 1] <- LAMBDA[2, 1] + delta_lambda_12 * hM1  # item2
    THETA <- THETA_base
    
  } else {
    stop("Unsupported model_type. Use 'NULL', '1.1', or '1.2'.")
  }
  
  # ---------------- Assemble return ----------------
  MLIST <- list(
    lambda = LAMBDA,
    theta  = THETA,
    psi    = PSI,
    beta   = BETA,
    moderators = list(
      moderator_1 = list(value = moderator_1_value,
                         type  = moderator_1_type,
                         hM    = hM1),
      moderator_2 = list(value = moderator_2_value,
                         type  = "noise",
                         hM    = hM2)
    )
  )
  
  return(MLIST)
}

# Utility (kept for style consistency; not used here)
lav_matrix_diag_idx <- function(n) {
  seq(1, n^2, by = n + 1)
}

# =============================================================================
# Case-wise simulator 
# Purpose:
#   Generate N observations for the single-factor model with per-case moderation.
#   Uses the same population equations as gen_mat(), but varies M per subject,
#   so methods like MNLFA and SEM Tree have a real moderation signal to detect.
# =============================================================================

# Build a named vector of true δ's for reporting/evaluation
.build_true_deltas <- function(model_type, delta_lambda_full, delta_theta_full, delta_lambda_12) {
  if (model_type == "1.1") {
    c(delta_lambda_y1 = delta_lambda_full,
      delta_lambda_y2 = delta_lambda_full,
      delta_lambda_y3 = delta_lambda_full,
      delta_lambda_y4 = delta_lambda_full,
      delta_theta_y1  = delta_theta_full,
      delta_theta_y2  = delta_theta_full,
      delta_theta_y3  = delta_theta_full,
      delta_theta_y4  = delta_theta_full)
  } else if (model_type == "1.2") {
    c(delta_lambda_y1 = delta_lambda_12,
      delta_lambda_y2 = delta_lambda_12,
      delta_lambda_y3 = 0,
      delta_lambda_y4 = 0)
  } else { # "NULL"
    c(delta_lambda_y1 = 0, delta_lambda_y2 = 0, delta_lambda_y3 = 0, delta_lambda_y4 = 0)
  }
}

# Case-wise generator
simulate_moderated_onefactor <- function(
    # --- design knobs (same naming as in gen_mat for consistency) --------------
    model_type = c("NULL", "1.1", "1.2"),
    N = 500,
    lambda = 0.70,                 # base loading for y1..y3; y4 is reference 1.0
    reliability = 0.80,            # target reliability at base (unmoderated)
    psi_var = 1.0,                 # Var(f1)
    moderator_1_type = c("linear", "sigmoid", "quadratic", "noise"),
    sigmoid_slope_1 = 2.0,
    delta_lambda_full = 0.20,      # used in 1.1 for ALL loadings
    delta_theta_full  = 0.20,      # used in 1.1 for ALL residual variances
    delta_lambda_12   = 0.20,      # used in 1.2 for y1 & y2 only
    moderate_reference_loading = TRUE,  # if FALSE: do not moderate y4’s loading
    # --- optional: provide your own moderators (else drawn as N(0,1)) ----------
    M  = NULL,   # informative moderator; if NULL, rnorm(N)
    M2 = NULL,   # noise moderator; always non-informative; if NULL, rnorm(N)
    # --- numerical guard --------------------------------------------------------
    min_theta = 1e-8
) {
  model_type       <- match.arg(model_type)
  moderator_1_type <- match.arg(moderator_1_type)
  
  # 0) Get the base (unmoderated) matrices using your gen_mat at h(M)=0
  base_pop <- gen_mat(
    model_type        = model_type,
    nfactors          = 1, nvar.factor = 4,
    lambda            = lambda,
    reliability       = reliability,
    psi_var           = psi_var,
    moderator_1_value = 0,                 # <-- key: base (h=0)
    moderator_1_type  = moderator_1_type,
    sigmoid_slope_1   = sigmoid_slope_1,
    delta_lambda_full = delta_lambda_full,
    delta_theta_full  = delta_theta_full,
    delta_lambda_12   = delta_lambda_12
  )
  lambda0 <- as.numeric(base_pop$lambda[, 1])  # y1..y4
  theta0  <- diag(base_pop$theta)              # y1..y4
  psi0    <- as.numeric(base_pop$psi)
  
  # 1) Draw moderators per case if not supplied
  if (is.null(M))  M  <- rnorm(N, 0, 1)
  if (is.null(M2)) M2 <- rnorm(N, 0, 1)       # explicit “noise” moderator
  hM <- .mod_h(M, type = moderator_1_type, slope = sigmoid_slope_1)  # same link as gen_mat()
  
  # 2) Build case-wise Λ(M) and Θ(M)
  lambda_mat <- matrix(rep(lambda0, each = N), nrow = N, ncol = 4)
  theta_mat  <- matrix(rep(theta0,  each = N), nrow = N, ncol = 4)
  
  if (model_type == "1.1") {
    # Loadings: additive shift on all indicators (optionally skip reference y4)
    add <- delta_lambda_full * hM
    if (moderate_reference_loading) {
      lambda_mat <- sweep(lambda_mat, 1, add, "+")
    } else {
      lambda_mat[, 1:3] <- sweep(lambda_mat[, 1:3, drop = FALSE], 1, add, "+")
    }
    # Residuals: multiplicative scaling θ(M) = θ0 * (1 + δθ * hM)
    scale <- pmax(1 + delta_theta_full * hM, .Machine$double.eps)
    theta_mat <- sweep(theta_mat, 1, scale, "*")
  } else if (model_type == "1.2") {
    add <- delta_lambda_12 * hM
    lambda_mat[, 1] <- lambda_mat[, 1] + add  # y1
    lambda_mat[, 2] <- lambda_mat[, 2] + add  # y2
    # Θ unchanged
  } # NULL: unchanged
  
  # Enforce positivity of residual variances
  theta_mat <- pmax(theta_mat, min_theta)
  
  # 3) Simulate latent η and residuals e, then y = Λ(M)*η + e
  eta <- rnorm(N, 0, sqrt(psi0))
  E   <- matrix(rnorm(N * 4), nrow = N, ncol = 4) * sqrt(theta_mat)
  Y   <- lambda_mat * eta + E
  colnames(Y) <- paste0("y", 1:4)
  
  # 4) Bundle outputs
  list(
    data = data.frame(Y, M = M, M2 = M2),
    # base values and truths (useful for metrics/sanity)
    base = list(lambda0 = lambda0, theta0 = theta0, psi = psi0),
    true_deltas = .build_true_deltas(model_type, delta_lambda_full, delta_theta_full, delta_lambda_12),
    knobs = list(
      model_type = model_type, N = N, lambda = lambda, reliability = reliability, psi_var = psi_var,
      moderator_1_type = moderator_1_type, sigmoid_slope_1 = sigmoid_slope_1,
      delta_lambda_full = delta_lambda_full, delta_theta_full = delta_theta_full, delta_lambda_12 = delta_lambda_12,
      moderate_reference_loading = moderate_reference_loading
    )
  )
}


# ------------- gen_pop_lavsyntax --------------------

# Function to round numbers to 3 digits
round_three <- function(x) {
  return(format(round(x, 3), nsmall = 3))
}

# Function to generate lavaan model syntax from model matrices
gen_pop_model_syntax <- function(MLIST, ov.prefix = "y", lv.prefix = "f", include.values = TRUE) {
  
  LAMBDA <- MLIST$lambda
  THETA  <- MLIST$theta
  PSI    <- MLIST$psi
  BETA   <- MLIST$beta
  
  if (ov.prefix == lv.prefix) {
    stop("lavaan ERROR: ov.prefix can not be the same as lv.prefix")
  }
  
  header <- "# syntax generated by gen_pop_model_syntax()"
  
  # LAMBDA
  if (!is.null(LAMBDA)) {
    IDXV <- row(LAMBDA)[(LAMBDA != 0)]
    IDXF <- col(LAMBDA)[(LAMBDA != 0)]
    
    unique_factors <- unique(IDXF)
    IDXV <- as.integer(unlist(sapply(unique_factors, function(j) {
      ji <- IDXV[which(IDXF == j)]
      j1 <- which(abs(LAMBDA[ji, j] - 1) < .Machine$double.eps)
      if (length(j1) > 0) {
        ji[c(1, j1)] <- ji[c(j1, 1)]
      }
      return(ji)
    })))
    
    IDXF <- rep(unique_factors, times = sapply(unique_factors, function(j) length(which(IDXF == j))))
    
    nel <- length(IDXF)
    lambda.txt <- character(nel)
    for (i in seq_len(nel)) {
      value <- LAMBDA[IDXV[i], IDXF[i]]
      if (include.values) {
        lambda.txt[i] <- sprintf("%s%d =~ %s*%s%d", lv.prefix, IDXF[i], round_three(value), ov.prefix, IDXV[i])
      } else {
        lambda.txt[i] <- sprintf("%s%d =~ %s%d", lv.prefix, IDXF[i], ov.prefix, IDXV[i])
      }
    }
  } else {
    lambda.txt <- character(0L)
  }
  
  # THETA
  if (!is.null(THETA)) {
    IDX1 <- row(THETA)[(THETA != 0) & upper.tri(THETA, diag = TRUE)]
    IDX2 <- col(THETA)[(THETA != 0) & upper.tri(THETA, diag = TRUE)]
    nel <- length(IDX1)
    theta.txt <- character(nel)
    for (i in seq_len(nel)) {
      value <- THETA[IDX1[i], IDX2[i]]
      if (include.values) {
        theta.txt[i] <- sprintf("%s%d ~~ %s*%s%d", ov.prefix, IDX1[i], round_three(value), ov.prefix, IDX2[i])
      } else {
        theta.txt[i] <- sprintf("%s%d ~~ %s%d", ov.prefix, IDX1[i], ov.prefix, IDX2[i])
      }
    }
  } else {
    theta.txt <- character(0L)
  }
  
  # PSI
  if (!is.null(PSI)) {
    IDX1 <- row(PSI)[(PSI != 0) & upper.tri(PSI, diag = TRUE)]
    IDX2 <- col(PSI)[(PSI != 0) & upper.tri(PSI, diag = TRUE)]
    nel <- length(IDX1)
    psi.txt <- character(nel)
    for (i in seq_len(nel)) {
      value <- PSI[IDX1[i], IDX2[i]]
      if (include.values) {
        psi.txt[i] <- sprintf("%s%d ~~ %s*%s%d", lv.prefix, IDX1[i], round_three(value), lv.prefix, IDX2[i])
      } else {
        psi.txt[i] <- sprintf("%s%d ~~ %s%d", lv.prefix, IDX1[i], lv.prefix, IDX2[i])
      }
    }
  } else {
    psi.txt <- character(0L)
  }
  
  # BETA
  if (!is.null(BETA)) {
    IDX1 <- row(BETA)[(BETA != 0)]
    IDX2 <- col(BETA)[(BETA != 0)]
    nel <- length(IDX1)
    beta.txt <- character(nel)
    for (i in seq_len(nel)) {
      value <- BETA[IDX1[i], IDX2[i]]
      if (include.values) {
        beta.txt[i] <- sprintf("%s%d ~ %s*%s%d", lv.prefix, IDX1[i], round_three(value), lv.prefix, IDX2[i])
      } else {
        beta.txt[i] <- sprintf("%s%d ~ %s%d", lv.prefix, IDX1[i], lv.prefix, IDX2[i])
      }
    }
  } else {
    beta.txt <- character(0L)
  }
  
  # Assemble
  syntax <- paste(c(header, lambda.txt, theta.txt, psi.txt, beta.txt, ""),
                  collapse = "\n")
  
  return(syntax)
}
