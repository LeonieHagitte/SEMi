# dataprep_functions.R
# -----------------------------------------------------------------------------

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
                    moderator_1_value = 1,
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

# =============================================================================
# Case-wise simulator 
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
    # --- optional: provide own moderators  ----------
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
    moderator_1_value = 1,                 # <-- key: base (h=0)
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
