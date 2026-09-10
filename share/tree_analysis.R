# =========================================================================
# CASP-12 SEM-tree analysis
# =========================================================================

casp_tree_analysis <- function(
    data,
    predictors = c("age_z", "female", "cultural_cluster"),
    alpha = 0.05,
    control = NULL,
    verbose = FALSE
) {
  
  # -----------------------------------------------------------------------
  # Data and variables
  # -----------------------------------------------------------------------
  
  dat <- as.data.frame(data)
  
  manVars <- c(
    "C1", "C2", "C3",
    "A1", "A2", "A3",
    "P1", "P2", "P3",
    "S1", "S2", "S3"
  )
  
  latVars <- c(
    "Control",
    "Autonomy",
    "Pleasure",
    "SelfRealization"
  )
  
  # Check required variables
  miss_man <- setdiff(manVars, names(dat))
  miss_pred <- setdiff(predictors, names(dat))
  
  if (length(miss_man) > 0) {
    stop(
      "Missing manifest variables: ",
      paste(miss_man, collapse = ", ")
    )
  }
  
  if (length(miss_pred) > 0) {
    stop(
      "Missing SEM-tree predictors: ",
      paste(miss_pred, collapse = ", ")
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Tree control
  # -----------------------------------------------------------------------
  
  if (is.null(control)) {
    
    control <- semtree::semtree_control(
      method = "score",
      alpha = alpha,
      max.depth = 3,
      bonferroni = TRUE,
      min.N = 100
    )
  }
  
  
  # -----------------------------------------------------------------------
  # OpenMx data
  #
  # Baseline CFA only needs manifest variables.
  # Predictors remain in `dat` for semtree().
  # -----------------------------------------------------------------------
  
  mx_dat <- dat[, manVars, drop = FALSE]
  
  mx_dat[manVars] <- lapply(
    mx_dat[manVars],
    as.numeric
  )
  
  mxdata <- mxData(
    observed = mx_dat,
    type = "raw"
  )
  
  
  # -----------------------------------------------------------------------
  # Manifest intercepts
  # -----------------------------------------------------------------------
  
  #path_nu <- mxPath(
   # from = "one",
   #to = manVars,
   # arrows = 1,
  #  free = TRUE,
  #  values = 3,
  #  labels = paste0("nu_", manVars)
  #)
  
  path_nu_free <- mxPath(
    from = "one",
    to = c(
      "C2", "C3",
      "A2", "A3",
      "P2", "P3",
      "S2", "S3"
    ),
    arrows = 1,
    free = TRUE,
    values = 3,
    labels = paste0(
      "nu_",
      c("C2", "C3", "A2", "A3", "P2", "P3", "S2", "S3")
    )
  )
  
  path_nu_anchor <- mxPath(
    from = "one",
    to = c("C1", "A1", "P1", "S1"),
    arrows = 1,
    free = FALSE,
    values = 0
  )
  
  # -----------------------------------------------------------------------
  # Factor loadings
  # -----------------------------------------------------------------------
  
  anchor_items <- c("C1", "A1", "P1", "S1")
  
  # Anchor loadings fixed to 1
  path_lambda_anchor <- mxPath(
    from = c(
      "Control",
      "Autonomy",
      "Pleasure",
      "SelfRealization"
    ),
    to = anchor_items,
    arrows = 1,
    free = FALSE,
    values = 1
  )
  
  # Non-anchor loadings freely estimated
  path_lambda_control <- mxPath(
    from = "Control",
    to = c("C2", "C3"),
    arrows = 1,
    free = TRUE,
    values = c(1.08, .99),
    labels = paste0("lambda_", c("C2", "C3"))
  )
  
  path_lambda_autonomy <- mxPath(
    from = "Autonomy",
    to = c("A2", "A3"),
    arrows = 1,
    free = TRUE,
    values = c(.51, 1.08),
    labels = paste0("lambda_", c("A2", "A3"))
  )
  
  path_lambda_pleasure <- mxPath(
    from = "Pleasure",
    to = c("P2", "P3"),
    arrows = 1,
    free = TRUE,
    values = c(1.08, .86),
    labels = paste0("lambda_", c("P2", "P3"))
  )
  
  path_lambda_selfrealization <- mxPath(
    from = "SelfRealization",
    to = c("S2", "S3"),
    arrows = 1,
    free = TRUE,
    values = c(1.14, 1.17),
    labels = paste0("lambda_", c("S2", "S3"))
  )
  
  
  # -----------------------------------------------------------------------
  # Residual variances
  # -----------------------------------------------------------------------
  
  path_resid <- mxPath(
    from = manVars,
    arrows = 2,
    free = TRUE,
    values = .5,
    labels = paste0("resid_", manVars)
  )
  
  
  # -----------------------------------------------------------------------
  # Latent variances free
  # -----------------------------------------------------------------------
  
  path_latvar <- mxPath(
    from = latVars,
    arrows = 2,
    free = TRUE,
    values = 1,
    labels = paste0("var_", latVars)
  )
  
  
  # -----------------------------------------------------------------------
  # Latent covariances freely estimated
  # -----------------------------------------------------------------------
  
  latent_pairs <- combn(latVars, 2)
  
  path_latcov <- mxPath(
    from = latent_pairs[1, ],
    to = latent_pairs[2, ],
    arrows = 2,
    free = TRUE,
    values = .3,
    labels = c(
      "cov_Control_Autonomy",
      "cov_Control_Pleasure",
      "cov_Control_SelfRealization",
      "cov_Autonomy_Pleasure",
      "cov_Autonomy_SelfRealization",
      "cov_Pleasure_SelfRealization"
    )
  )
  
  
  # -----------------------------------------------------------------------
  # Latent means freely estimated
  # -----------------------------------------------------------------------
  
  path_latmean <- mxPath(
    from = "one",
    to = latVars,
    arrows = 1,
    free = TRUE,
    values = 0,
    labels = paste0("mean_", latVars)
  )
  
  
  # -----------------------------------------------------------------------
  # Baseline four-factor CASP model
  # -----------------------------------------------------------------------
  
  fitF <- mxFitFunctionML()
  
  modbase <- mxModel(
    model = "CASP12_baseline",
    type = "RAM",
    
    manifestVars = manVars,
    latentVars = latVars,
    
    path_nu_free,
    path_nu_anchor,
    
    path_lambda_anchor,
    path_lambda_control,
    path_lambda_autonomy,
    path_lambda_pleasure,
    path_lambda_selfrealization,
    
    path_resid,
    
    path_latvar,
    path_latcov,
    path_latmean,
    
    fitF,
    mxdata
  )
  
  
  # -----------------------------------------------------------------------
  # Fit baseline model
  # -----------------------------------------------------------------------
  
  fitbase <- mxRun(modbase)
  
  
  # -----------------------------------------------------------------------
  # Check convergence before attempting trees
  # -----------------------------------------------------------------------
  
  if (
    is.null(fitbase$output$status$code) ||
    fitbase$output$status$code != 0
  ) {
    stop(
      "Baseline OpenMx CASP model did not converge. Status code: ",
      fitbase$output$status$code
    )
  }
  #return(fitbase)
  
  # -----------------------------------------------------------------------
  # Metric MNI: factor loadings
  # -----------------------------------------------------------------------
  anchor_items <- c("C1", "A1", "P1", "S1")
  
  metric_items <- c(
    "C2", "C3",
    "A2", "A3",
    "P2", "P3",
    "S2", "S3"
  )
  
  metric_constraints <- semtree::semtree.constraints(
    focus.parameters = paste0("lambda_", metric_items)
  )
  
  metric_tree <- tryCatch(
    semtree::semtree(
      model = fitbase,
      data = dat,
      predictors = predictors,
      control = control,
      constraints = metric_constraints,
      verbose = verbose
    ),
    error = identity
  )
  
  metric_test <- extract_tree_test(
    metric_tree,
    alpha = alpha
  )
  
  
  # -----------------------------------------------------------------------
  # Scalar MNI: item intercepts
  # -----------------------------------------------------------------------
  
  scalar_items <- c(
    "C2", "C3",
    "A2", "A3",
    "P2", "P3",
    "S2", "S3"
  )
  
  scalar_constraints <- semtree::semtree.constraints(
    focus.parameters = paste0("nu_", scalar_items)
  )
  
  scalar_tree <- tryCatch(
    semtree::semtree(
      model = fitbase,
      data = dat,
      predictors = predictors,
      control = control,
      constraints = scalar_constraints,
      verbose = verbose
    ),
    error = identity
  )
  
  scalar_test <- extract_tree_test(
    scalar_tree,
    alpha = alpha
  )
  
  
  # -----------------------------------------------------------------------
  # Output
  # -----------------------------------------------------------------------
  
  return(
    list(
      baseline_model = modbase,
      baseline_fit = fitbase,
      
      metric_tree = metric_tree,
      metric_test = metric_test,
      
      scalar_tree = scalar_tree,
      scalar_test = scalar_test
    )
  )
}




fit_casp_openmx2 <- casp_tree_analysis(
  data = tree_data
)
#--------------------------------------------------
fit_casp_openmx2$output$status$code
summary(fit_casp_openmx2)
omxGetParameters(fit_casp_openmx2, free = TRUE)
fit_casp_openmx2$output$Minus2LogLikelihood


fit_casp_openmx2$metric_tree
fit_casp_openmx2$metric_test

fit_casp_openmx2$scalar_tree
fit_casp_openmx2$scalar_test

# Metric tree
plot(fit_casp_openmx2$metric_tree)

# Scalar tree
plot(fit_casp_openmx2$scalar_tree)

# OpenMx convergence status
fit_casp_openmx2$baseline_fit$output$status$code

# -2 Log Likelihood
fit_casp_openmx2$baseline_fit$output$Minus2LogLikelihood

# Parameter estimates
omxGetParameters(
  fit_casp_openmx2$baseline_fit,
  free = TRUE
)

# Full OpenMx summary
summary(fit_casp_openmx2$baseline_fit)

#--------------
tree_baseline_summary <- tibble::tibble(
  N = fit_casp_openmx2$baseline_fit$data$numObs,
  Parameters = summary(
    fit_casp_openmx2$baseline_fit
  )$estimatedParameters,
  minus2LL =
    fit_casp_openmx2$baseline_fit$output$Minus2LogLikelihood,
  AIC = AIC(fit_casp_openmx2$baseline_fit),
  BIC = BIC(fit_casp_openmx2$baseline_fit),
  Status =
    fit_casp_openmx2$baseline_fit$output$status$code
)

tree_baseline_summary
#-----------------------------------------------------------------------------

