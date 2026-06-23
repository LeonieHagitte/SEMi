library(lavaan)
library(dplyr)
library(purrr)
library(mxsem)
library(OpenMx)
library(semtree)
library(partykit)

# -----------------------------------------------------------------------------


compare_models_lrt <- function(fit_less, fit_more, alpha = 0.05) {
  cmp <- tryCatch(OpenMx::mxCompare(fit_less, fit_more), error = identity)
  
  if (inherits(cmp, "error")) {
    return(list(
      chisq_diff = NA_real_,
      df_diff = NA_real_,
      p_value = NA_real_,
      reject_h0 = NA,
      lrt_status = "error",
      lrt_note = conditionMessage(cmp)
    ))
  }
  
  # Coerce safely to data.frame if possible
  cmp_df <- tryCatch(as.data.frame(cmp), error = function(e) NULL)
  
  if (is.null(cmp_df) || nrow(cmp_df) < 2) {
    return(list(
      chisq_diff = NA_real_,
      df_diff = NA_real_,
      p_value = NA_real_,
      reject_h0 = NA,
      lrt_status = "unavailable",
      lrt_note = "mxCompare did not return a usable comparison table."
    ))
  }
  
  # second row is the more constrained model in a 2-model comparison
  chisq_diff <- suppressWarnings(as.numeric(cmp_df$diffLL[2]))
  df_diff    <- suppressWarnings(as.numeric(cmp_df$diffdf[2]))
  p_value    <- suppressWarnings(as.numeric(cmp_df$p[2]))
  
  list(
    chisq_diff = chisq_diff,
    df_diff = df_diff,
    p_value = p_value,
    reject_h0 = if (!is.na(p_value)) p_value <= alpha else NA,
    lrt_status = "ok",
    lrt_note = NA_character_
  )
}

extract_tree_test <- function(st, alpha = 0.05) {
  if (inherits(st, "error") || is.null(st) || !methods::is(st, "semtree")) {
    return(list(
      p_value = NA_real_,
      p_uncorrected = NA_real_,
      reject_h0 = NA,
      tree_test_status = "error"
    ))
  }
  
  p_corr <- if (!is.null(st$p)) as.numeric(st$p) else NA_real_
  p_unc  <- if (!is.null(st$p.uncorrected)) as.numeric(st$p.uncorrected) else NA_real_
  
  list(
    p_value = p_corr,
    p_uncorrected = p_unc,
    reject_h0 = if (!is.na(p_corr)) p_corr <= alpha else NA,
    tree_test_status = "ok"
  )
}

# -------------extract moderation parameters MNLFA ----------------------------


extract_mnlfa_moderation <- function(fit, p) {
  if (is.null(fit) || inherits(fit, "error")) return(NULL)
  
  est_free <- OpenMx::omxGetParameters(fit, free = TRUE)
  est_all  <- OpenMx::omxGetParameters(fit, free = FALSE)
  
  get_param <- function(prefix, row, col) {
    nm <- paste0(prefix, "_", row, "_", col)
    
    if (nm %in% names(est_free)) return(unname(est_free[nm]))
    if (nm %in% names(est_all))  return(unname(est_all[nm]))
    
    NA_real_
  }
  
  tibble::tibble(
    item = paste0("x", seq_len(p)),
    
    mnlfa_est_dnu_am1      = purrr::map_dbl(seq_len(p), ~ get_param("matB1",  1, .x)),
    mnlfa_est_dnu_am2      = purrr::map_dbl(seq_len(p), ~ get_param("matB2",  1, .x)),
    mnlfa_est_dnu_am12     = purrr::map_dbl(seq_len(p), ~ get_param("matB12", 1, .x)),
    
    mnlfa_est_dlambda_am1  = purrr::map_dbl(seq_len(p), ~ get_param("matC1",  .x, 1)),
    mnlfa_est_dlambda_am2  = purrr::map_dbl(seq_len(p), ~ get_param("matC2",  .x, 1)),
    mnlfa_est_dlambda_am12 = purrr::map_dbl(seq_len(p), ~ get_param("matC12", .x, 1))
  )
}


# ---------------------- unrestricted CFA (SEM Tree base) NULL -----------------
mnlfa_analysis <- function(data, p = 4, nfactors = 1, alpha = 0.05) {
  
  manVars <- grep("^x\\d+$", names(data), value = TRUE)
  p <- length(manVars)
  
  mxdata <- mxData(observed = data, type = "raw")
  
  # Intercept Matrices
  matT0 <- mxMatrix(type = "Full", 
                    nrow = 1, 
                    ncol = p, 
                    free = TRUE, 
                    values = 1,
                    name = "matT0") #The matrix matT0 is a full matrix containing the baseline intercepts.
  #All baseline intercepts are freely estimated with starting values of one by setting free=TRUE and values=1. 
 
  matB1 <- mxMatrix(type="Full", 
                    nrow = 1, 
                    ncol=p,#Matrix matB1 and matB2 are full matrices containing the direct effects of the background variables M1 and M2, respectively, on the intercepts.
                    free=TRUE, 
                    values = 0, 
                    labels = paste0("matB1_1_", seq_len(p)),
                    name="matB1")
  
  matB2 <- mxMatrix(type="Full", 
                    nrow = 1, 
                    ncol=p,
                    free=TRUE, # shows that effects of m2 in this case are freely estimated
                    values = 0,
                    labels = paste0("matB2_1_", seq_len(p)),
                    name="matB2")
  
  matB12 <- mxMatrix(type = "Full", #Matrix for interaction effects on intercepts
                     nrow = 1,
                     ncol = p,
                     free = TRUE,
                     values = 0,
                     labels = paste0("matB12_1_", seq_len(p)),
                     name = "matB12")
  
  # Loading Matrices
  matL0 <- mxMatrix(type="Full", #matL0 is a full matrix containing the baseline factor loadings
                    nrow=p, 
                    ncol = nfactors, 
                    free= TRUE, 
                    values= rep(1, p),
                    name="matL0")
  
  matC1 <- mxMatrix(type="Full", # Matrices matC1 and matC2 are full matrices containing the direct effects of M1 and M2, respectively, on the factor loadings
                    nrow=p,
                    ncol = nfactors,
                    free= TRUE, 
                    values = 0, # value sets the moderation effect
                    labels = paste0("matC1_", seq_len(p), "_1"),
                    name="matC1")
  
  matC2 <- mxMatrix(type="Full", 
                    nrow=p,
                    ncol = nfactors,
                    free= TRUE, 
                    values = 0,
                    labels = paste0("matC2_", seq_len(p), "_1"),
                    name="matC2")
  
  matC12 <- mxMatrix(type = "Full",
                     nrow = p,
                     ncol = nfactors,
                     free = TRUE,
                     values = 0,
                     labels = paste0("matC12_", seq_len(p), "_1"),
                     name = "matC12")
  
  # Residual Variances 
  matE0 <- mxMatrix(type="Diag", # matE0 is a diagonal matrix containing the baseline residual variances 
                    nrow=p, 
                    ncol=p,
                    free=TRUE,
                    values = 1,
                    name="matE0")
  
  matD1 <- mxMatrix(type="Diag", #matD1 is a diagonal matrix containing the effects of M1 on the residual variances
                    nrow=p, 
                    ncol=p,
                    free=TRUE,
                    values = 0,
                    name="matD1")
  
  matD2 <- mxMatrix(type="Diag",  #and matD2 is a diagonal matrix containing the effects of M2 on the residual variances
                    nrow=p, 
                    ncol=p,
                    free=TRUE,
                    values = 0,
                    name="matD2")
  
  # Latent factor variance
  matP0 <- mxMatrix(type="Symm", 
                    nrow = 1, #because only one latent factor 
                    ncol = 1,
                    free= FALSE,
                    values= 1,
                    name="matP0")
  
  # latent factor mean
  matA0 <- mxMatrix(type="Full", #matA0 is a matrix containing the baseline common-factor means 
                    nrow = 1, 
                    ncol = 1,
                    free=FALSE, #fixed
                    values = 0, # to 0
                    name="matA0")
  
  matG1 <- mxMatrix(type="Full", #The matG1 and matG2 matrices contain the direct effects of M1 and M2, respectively, on the common-factor means 
                    nrow = 1, 
                    ncol = 1,
                    free=FALSE, #fixed
                    values = 0, #to 0
                    name="matG1")
  
  matG2 <- mxMatrix(type="Full", 
                    nrow = 1, 
                    ncol = 1,
                    free=FALSE, #fixed
                    values = 0, #to 0
                    name="matG2")
  
  # ---------- Matrices for Matrix algebra
  # Moderators as Definition Variables
  matV1 <- mxMatrix(type="Full", 
                    nrow = 1, 
                    ncol = 1,
                    free = FALSE,
                    labels ="data.am1",
                    name ="am1")
  
  matV2 <- mxMatrix(type="Full", 
                    nrow = 1, 
                    ncol = 1,
                    free = FALSE,
                    labels ="data.am2",
                    name ="am2")
  
  matV12 <- mxMatrix(type = "Full",
                     nrow = 1,
                     ncol = 1,
                     free = FALSE,
                     labels = "data.am12",
                     name = "am12")
  
  # -----------------------------------------------------------------------------
  #The mxAlgebra() function can be used to define a matrix of model parameters as
  #a function of background variables. The first argument of this function, 
  #expression, should be used for specifying an R expression of one or more 
  #MxMatrix objects. Most R operators like +, ∗, and %∗%, an general R functions 
  #like mean(), log(), and exp() are supported in this argument of the mxAlgebra() 
  #function. A name for the defined matrix can be assigned with the name argument
  
  # -----------------------------------------------------------------------------
  matT <- mxAlgebra(expression = matT0 +  matB1*am1 + matB2*am2 + matB12*am12, # linear moderation function for the indicator intercepts
                    name="matT") #intercepts
  
  matL <- mxAlgebra(expression = matL0 + matC1*am1 + matC2*am2 + matC12*am12, # linear moderation function for the factor loadings
                    name="matL") #loadings
  
  matE <- mxAlgebra(expression = matE0*exp (matD1*am1 + matD2*am2), # log-linear function for the residual variances
                    name="matE") #residual variances
  
  matA <- mxAlgebra(expression = matA0 + matG1*am1 + matG2*am2, #common-factor means modeled as a linear function of the background variables
                    name="matA") #matrix of common-factor means
  
  matP <- mxAlgebra(matP0, name="matP") #latent covariance matrix in the covariance algebra - for a single-factor CFA, that matrix is simply 1×1 and is the factor variance
  
  # ------------ Matrices for model implied moments -------
  matM <- mxAlgebra(matT + matA %*% t(matL), name = "matM")
  
  matC <- mxAlgebra(expression = matL %*% matP%*%t(matL) + matE ,
                    name="matC")
  
  # ------------ Model expectations and fit-function -------
  expF <- mxExpectationNormal(covariance="matC", means = "matM", dimnames=manVars)
  
  fitF <- mxFitFunctionML() #mxFitFunctionML() function stored in fitF is used to 
  #indicate that the free parameters of the configural 
  #model should be estimated using full-information maximum likelihood.
  
  modConfig <- mxModel(model="Configural",
                        matT, matT0, matB1, 
                        matB2, matB12, matL, matL0,
                        matC1, matC2, matC12, matE, 
                        matE0, matD1, matD2,
                        matP, matP0, 
                        matA, matA0, 
                        matG1, matG2,
                        matV1, matV2, matV12,
                        matM, matC, expF, 
                        fitF, mxdata)
  
  #The model can be fitted to the data using the mxRun() function. 
  fitConfig <- mxRun(modConfig)
  
  configural_moderation_estimates <- NULL
  if (!is.null(fitConfig$output$status$code) && fitConfig$output$status$code == 0) {
    configural_moderation_estimates <- extract_mnlfa_moderation(
      fit = fitConfig,
      p = p
    )
  }
  
  if (!is.null(fitConfig$output$status$code) && fitConfig$output$status$code == 0) {
    # ---------------------- metric moderated model --------------------------
    
    # Intercept Matrices
    
    matB1 <- mxMatrix(type="Full", 
                      nrow = 1, 
                      ncol=p,#Matrix matB1 and matB2 are full matrices containing the direct effects of the background variables M1 and M2, respectively, on the intercepts.
                      free=TRUE, 
                      values = 0,
                      labels = paste0("matB1_1_", seq_len(p)),
                      name="matB1")
    
    matB2 <- mxMatrix(type="Full", 
                      nrow = 1, 
                      ncol=p,
                      free=TRUE, 
                      values = 0,
                      labels = paste0("matB2_1_", seq_len(p)),
                      name="matB2")
    
    matB12 <- mxMatrix(type = "Full", #Matrix for interaction effects on intercepts
                       nrow = 1,
                       ncol = p,
                       free = TRUE,
                       values = 0,
                       labels = paste0("matB12_1_", seq_len(p)),
                       name = "matB12")
    
    # Loading Matrices
    
    matC1 <- mxMatrix(type="Full", # Matrices matC1 and matC2 are full matrices containing the direct effects of M1 and M2, respectively, on the factor loadings
                      nrow=p,
                      ncol = nfactors,
                      free= FALSE, 
                      values = 0, # value sets the moderation effect
                      labels = paste0("matC1_", seq_len(p), "_1"),
                      name="matC1")
    
    matC2 <- mxMatrix(type="Full", 
                      nrow=p,
                      ncol = nfactors,
                      free= FALSE, 
                      values = 0,
                      labels = paste0("matC2_", seq_len(p), "_1"),
                      name="matC2")
    
    matC12 <- mxMatrix(type = "Full",
                       nrow = p,
                       ncol = nfactors,
                       free = FALSE,
                       values = 0,
                       labels = paste0("matC12_", seq_len(p), "_1"),
                       name = "matC12")
    
    # -----------------------------------------------------------------------------
    matT <- mxAlgebra(expression = matT0 +  matB1*am1 + matB2*am2 + matB12*am12, # linear moderation function for the indicator intercepts
                      name="matT") #intercepts
    
    matL <- mxAlgebra(expression = matL0 + matC1*am1 + matC2*am2 + matC12*am12, # linear moderation function for the factor loadings
                      name="matL") #loadings
    
    matE <- mxAlgebra(expression = matE0*exp (matD1*am1 + matD2*am2), # log-linear function for the residual variances
                      name="matE") #residual variances
    
    matA <- mxAlgebra(expression = matA0 + matG1*am1 + matG2*am2, #common-factor means modeled as a linear function of the background variables
                      name="matA") #matrix of common-factor means
    
    matP <- mxAlgebra(matP0, name="matP") #latent covariance matrix in the covariance algebra - for a single-factor CFA, that matrix is simply 1×1 and is the factor variance
    
    # ------------ Matrices for model implied moments -------
    matM <- mxAlgebra(matT + matA %*% t(matL), name = "matM")
    
    matC <- mxAlgebra(expression = matL %*% matP%*%t(matL) + matE ,
                      name="matC")
    
    # ------------ Model expectations and fit-function -------
    expF <- mxExpectationNormal(covariance="matC", means = "matM", dimnames=manVars)
    
    fitF <- mxFitFunctionML() #mxFitFunctionML() function stored in fitF is used to 
    #indicate that the free parameters of the configural 
    #model should be estimated using full-information maximum likelihood.
    
    modmetric <- mxModel(model="Metric",
                          matT, matT0, matB1, 
                          matB2, matB12, matL, matL0,
                          matC1, matC2, matC12, matE, 
                          matE0, matD1, matD2,
                          matP, matP0, 
                          matA, matA0, 
                          matG1, matG2,
                          matV1, matV2, matV12,
                          matM, matC, expF, 
                          fitF, mxdata)
    
    #The model can be fitted to the data using the mxRun() function. 
    fitmetric <- mxRun(modmetric)
    if (is.null(fitmetric$output$status$code) || fitmetric$output$status$code != 0) {
      stop("scalar not run because metric model failed to converge.")
    }
    
    metric_lrt <- compare_models_lrt(fitConfig, fitmetric, alpha = alpha)
    
    if (!is.null(fitmetric) &&
        !is.null(fitmetric$output$status$code) &&
        fitmetric$output$status$code == 0) {
      # ---------------------- Scalar moderated model --------------------------
      
      # Intercept Matrices
      
      matB1 <- mxMatrix(type="Full", 
                        nrow = 1, 
                        ncol=p,#Matrix matB1 and matB2 are full matrices containing the direct effects of the background variables M1 and M2, respectively, on the intercepts.
                        free=FALSE, 
                        values = 0, 
                        labels = paste0("matB1_1_", seq_len(p)),
                        name="matB1")
      
      matB2 <- mxMatrix(type="Full", 
                        nrow = 1, 
                        ncol=p,
                        free=FALSE, 
                        values = 0, 
                        labels = paste0("matB2_1_", seq_len(p)),
                        name="matB2")
      
      matB12 <- mxMatrix(type = "Full", #Matrix for interaction effects on intercepts
                         nrow = 1,
                         ncol = p,
                         free = FALSE,
                         values = 0,
                         labels = paste0("matB12_1_", seq_len(p)),
                         name = "matB12")
      
      # Loading Matrices
      
      matC1 <- mxMatrix(type="Full", # Matrices matC1 and matC2 are full matrices containing the direct effects of M1 and M2, respectively, on the factor loadings
                        nrow=p,
                        ncol = nfactors,
                        free= FALSE, 
                        values = 0, # value sets the moderation effect
                        labels = paste0("matC1_", seq_len(p), "_1"),
                        name="matC1")
      
      matC2 <- mxMatrix(type="Full", 
                        nrow=p,
                        ncol = nfactors,
                        free= FALSE, 
                        values = 0,
                        labels = paste0("matC2_", seq_len(p), "_1"),
                        name="matC2")
      
      matC12 <- mxMatrix(type = "Full",
                         nrow = p,
                         ncol = nfactors,
                         free = FALSE,
                         values = 0,
                         labels = paste0("matC12_", seq_len(p), "_1"),
                         name = "matC12")
      
      matG1 <- mxMatrix(type="Full", #The matG1 and matG2 matrices contain the direct effects of M1 and M2, respectively, on the common-factor means 
                        nrow = 1, 
                        ncol = 1,
                        free=TRUE, 
                        values = 0, 
                        name="matG1")
      
      matG2 <- mxMatrix(type="Full", 
                        nrow = 1, 
                        ncol = 1,
                        free=TRUE, 
                        values = 0, 
                        name="matG2")
      
      # -----------------------------------------------------------------------------
      matT <- mxAlgebra(expression = matT0 +  matB1*am1 + matB2*am2 + matB12*am12, # linear moderation function for the indicator intercepts
                        name="matT") #intercepts
      
      matL <- mxAlgebra(expression = matL0 + matC1*am1 + matC2*am2 + matC12*am12, # linear moderation function for the factor loadings
                        name="matL") #loadings
      
      matE <- mxAlgebra(expression = matE0*exp (matD1*am1 + matD2*am2), # log-linear function for the residual variances
                        name="matE") #residual variances
      
      matA <- mxAlgebra(expression = matA0 + matG1*am1 + matG2*am2, #common-factor means modeled as a linear function of the background variables
                        name="matA") #matrix of common-factor means
      
      matP <- mxAlgebra(matP0, name="matP") #latent covariance matrix in the covariance algebra - for a single-factor CFA, that matrix is simply 1×1 and is the factor variance
      
      # ------------ Matrices for model implied moments -------
      matM <- mxAlgebra(matT + matA %*% t(matL), name = "matM")
      
      matC <- mxAlgebra(expression = matL %*% matP%*%t(matL) + matE ,
                        name="matC")
      
      # ------------ Model expectations and fit-function -------
      expF <- mxExpectationNormal(covariance="matC", means = "matM", dimnames=manVars)
      
      fitF <- mxFitFunctionML() #mxFitFunctionML() function stored in fitF is used to 
      #indicate that the free parameters of the configural 
      #model should be estimated using full-information maximum likelihood.
      
      modscalar <- mxModel(model="Scalar",
                            matT, matT0, matB1, 
                            matB2, matB12, matL, matL0,
                            matC1, matC2, matC12, matE, 
                            matE0, matD1, matD2,
                            matP, matP0, 
                            matA, matA0, 
                            matG1, matG2,
                            matV1, matV2, matV12,
                            matM, matC, expF, 
                            fitF, mxdata)
      
      #The model can be fitted to the data using the mxRun() function. 
      fitscalar <- mxRun(modscalar)
      if (is.null(fitscalar$output$status$code) || fitscalar$output$status$code != 0) {
        stop("Scalar model did not converge.")
      }
      
      scalar_lrt  <- compare_models_lrt(fitmetric, fitscalar, alpha = alpha)
      omnibus_lrt <- compare_models_lrt(fitConfig, fitscalar, alpha = alpha)
      
    }
  }
  return(list(
    fitConfig = fitConfig,
    fitMetric = if (exists("fitmetric")) fitmetric else NULL,
    fitScalar = if (exists("fitscalar")) fitscalar else NULL,
    metric_lrt = if (exists("metric_lrt")) metric_lrt else NULL,
    scalar_lrt = if (exists("scalar_lrt")) scalar_lrt else NULL,
    omnibus_lrt = if (exists("omnibus_lrt")) omnibus_lrt else NULL,
    configural_moderation_estimates = configural_moderation_estimates
    
  ))
  }

# -----------------------------------------------------#####################################
tree_analysis_ram <- function(data, p = 4, alpha = 0.05, nfactors = 1,
                              predictors = c("am1", "am2","m0"),
                              control = NULL, verbose = FALSE){
  
  if (is.null(control)) {
    control <- semtree::semtree.control(
      method = "score",
      alpha = alpha,
      max.depth = 3,
      bonferroni = TRUE,
      min.N = 50
    )}
  
    stopifnot(nfactors == 1)  # TODO: extend latent covariance structure for nfactors > 1
  
  dat <- as.data.frame(data)
  
  manVars <- grep("^x\\d+$", names(dat), value = TRUE)
  p <- length(manVars)
  
  if (p == 0) {
    stop("No manifest variables found. Expected columns like x1, x2, ...")
  }
  
  miss_p <- setdiff(predictors, names(dat))
  if (length(miss_p) > 0) {
    stop("Missing SEMTREE predictor columns: ", paste(miss_p, collapse = ", "))
  }
  
  latVars <- "F1"
  
  mxdata <- mxData(observed = dat, type = "raw")
  
  # ---------------- RAM paths ----------------
  # Manifest intercepts (nu_i)
  path_nu <- mxPath(
    from = "one",
    to = manVars,
    arrows = 1,
    free = TRUE,
    values = rep(0.6, p),
    labels = paste0("nu_", 1:p)
  )
  
  # Factor loadings (lambda_i)
  path_lambda <- mxPath(
    from = latVars,
    to = manVars,
    arrows = 1,
    free = TRUE,
    values = rep(1, p),
    labels = paste0("lambda_", 1:p)
  )
  
  # Residual variances
  path_resid <- mxPath(
    from = manVars,
    arrows = 2,
    free = TRUE,
    values = rep(1, p)
  )
  
  # Latent variance fixed to 1
  path_latvar <- mxPath(
    from = latVars,
    arrows = 2,
    free = FALSE,
    values = 1
  )
  
  # Latent mean fixed to 0
  path_latmean <- mxPath(
    from = "one",
    to = latVars,
    arrows = 1,
    free = FALSE,
    values = 0
  )
  
  fitTr <- mxFitFunctionML()
  
  modbase <- mxModel(
    model = "baseline",
    type = "RAM",
    manifestVars = manVars,
    latentVars = latVars,
    path_nu,
    path_lambda,
    path_resid,
    path_latvar,
    path_latmean,
    fitTr,
    mxdata
  )
  
  # The model can be fitted to the data using mxRun()
  fitbase <- mxRun(modbase)
  
  # ---------------- Metric stage: loadings ----------------
  metric_constraints <- semtree::semtree.constraints(
    focus.parameters = paste0("lambda_", 1:p) #c("lambda_1", "lambda_2", "lambda_3", "lambda_4")#########################flag
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
  
  #####
 # if (!inherits(metric_tree, "error")) {
  #  cat("\n--- metric_tree diagnostics ---\n")
  #  print(class(metric_tree))
    
  # if (isS4(metric_tree)) {
  #    cat("slotNames(metric_tree):\n")
  #    print(methods::slotNames(metric_tree))
  #  }
    
  #  cat("names(metric_tree):\n")
  #  print(names(metric_tree))
    
  #  str(metric_tree, max.level = 2)
  #}
  #####
  metric_test <- extract_tree_test(metric_tree, alpha = alpha)
  metric_split <- NA
  if (!inherits(metric_tree, "error") && methods::is(metric_tree, "semtree")) {
    metric_split <- !is.null(metric_tree$caption) &&
      !identical(metric_tree$caption, "TERMINAL")
  }
  
  scalar_tree <- NULL
  scalar_split <- NA
  scalar_test <- NULL
  
  # ---------------- Scalar stage: intercepts ----------------
 if (!isTRUE(metric_split)) {
    scalar_constraints <- semtree::semtree.constraints(
      focus.parameters = paste0("nu_", 1:p)
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
    
    #####
   # if (!inherits(scalar_tree, "error")) {
  #    cat("\n--- scalar_tree diagnostics ---\n")
  #    print(class(scalar_tree))
  #    
  #    if (isS4(scalar_tree)) {
  #      cat("slotNames(scalar_tree):\n")
  #      print(methods::slotNames(scalar_tree))
  #    }
      
  #   cat("names(scalar_tree):\n")
  #    print(names(scalar_tree))
  #    
  #    str(scalar_tree, max.level = 2)
  #  }
    #####
    
    if (!inherits(scalar_tree, "error") && methods::is(scalar_tree, "semtree")) {
      scalar_split <- !is.null(scalar_tree$caption) &&
        !identical(scalar_tree$caption, "TERMINAL")
    }
    scalar_test <- extract_tree_test(scalar_tree, alpha = alpha)
  }
  
  return(list(
    baseline_model = modbase,
    baseline_fit = fitbase,
    metric_tree = metric_tree,
    metric_split = metric_split,
    metric_test = metric_test,
    scalar_tree = scalar_tree,
    scalar_split = scalar_split,
    scalar_test = if (exists("scalar_test")) scalar_test else NULL
  ))
}
# ------------------------------------------------------------------------

run_analysis <- function(data,
                         methods = c("MNLFA", "SEMTREE"),
                         nfactors = 1,
                         alpha = 0.05,
                         predictors = c("am1", "am2", "m0")) {
  
  methods <- match.arg(methods, choices = c("MNLFA", "SEMTREE"), several.ok = TRUE)
  dat <- as.data.frame(data)
  
  out <- list(methods = methods)
  
  if ("MNLFA" %in% methods) {
    out$mnlfa <- tryCatch(
      mnlfa_analysis(data = dat, nfactors = nfactors, alpha = alpha),
      error = identity
    )
  }
  
  if ("SEMTREE" %in% methods) {
    out$semtree <- tryCatch(
      tree_analysis_ram(data = dat, nfactors = nfactors, alpha = alpha, predictors = predictors),
      error = identity
    )
  }
  
  return(out)
}
# ----------------------------------------------------------------------

#semtree_detects_moderation <- function(st, moderators = c("am1", "am2", "m0")) {
#  
#  out <- list()
#  for (m in moderators) {
#    out[[paste0("tree_split_on_", m)]] <- NA
#    out[[paste0("tree_n_splits_", m)]] <- NA_integer_
#  }
#  
#  if (inherits(st, "error") || is.null(st)) return(out)
#  if (!methods::is(st, "semtree")) return(out)
#  
#  # No split occurred
#  if (!is.null(st$caption) && identical(st$caption, "TERMINAL")) {
#    for (m in moderators) {
#      out[[paste0("tree_split_on_", m)]] <- FALSE
#      out[[paste0("tree_n_splits_", m)]] <- 0L
#    }
#    return(out)
#  }
#  
#  # For the current semtree object structure, a non-terminal root indicates a split.
#  split_var <- NULL
#  if (!is.null(st$result) && !is.null(st$result$name.max)) {
#    split_var <- st$result$name.max
#  }
#  
#  for (m in moderators) {
#    out[[paste0("tree_split_on_", m)]] <- identical(split_var, m)
#    out[[paste0("tree_n_splits_", m)]] <- as.integer(identical(split_var, m))
#  }
#  
#  out
#}
semtree_detects_moderation <- function(st, moderators = c("am1", "am2", "m0")) {
  
  out <- list()
  for (m in moderators) {
    out[[paste0("tree_split_on_", m)]] <- NA
    out[[paste0("tree_n_splits_", m)]] <- NA_integer_
  }
  
  if (inherits(st, "error") || is.null(st)) return(out)
  if (!methods::is(st, "semtree")) return(out)
  
  collect_splits <- function(node) {
    if (is.null(node)) return(character(0))
    
    splits <- character(0)
    
    if (methods::is(node, "semtree")) {
      if (!is.null(node$caption) && identical(node$caption, "TERMINAL")) {
        return(character(0))
      }
      
      if (!is.null(node$result) && !is.null(node$result$name.max)) {
        splits <- c(splits, as.character(node$result$name.max))
      }
    }
    
    if (is.list(node)) {
      child_splits <- unlist(
        lapply(node, collect_splits),
        use.names = FALSE
      )
      splits <- c(splits, child_splits)
    }
    
    splits
  }
  
  split_vars <- collect_splits(st)
  
  for (m in moderators) {
    n_m <- sum(split_vars == m, na.rm = TRUE)
    out[[paste0("tree_split_on_", m)]] <- n_m > 0
    out[[paste0("tree_n_splits_", m)]] <- as.integer(n_m)
  }
  
  out
}
# -----------------------------------------------------------------------
#getPredictorsFromTree <- function(st) {
#  if (inherits(st, "error") || is.null(st)) return(NULL)
#  if (!methods::is(st, "semtree")) return(NULL)
#  
#  if (!is.null(st$caption) && identical(st$caption, "TERMINAL")) {
#    return(character(0))
#  }
#  
#  if (!is.null(st$result) && !is.null(st$result$name.max)) {
#    return(st$result$name.max)
#  }
#  
#  NULL
#}
getPredictorsFromTree <- function(st) {
  if (inherits(st, "error") || is.null(st)) return(NULL)
  if (!methods::is(st, "semtree")) return(NULL)
  
  collect_splits <- function(node) {
    if (is.null(node)) return(character(0))
    
    splits <- character(0)
    
    if (methods::is(node, "semtree")) {
      if (!is.null(node$caption) && identical(node$caption, "TERMINAL")) {
        return(character(0))
      }
      
      if (!is.null(node$result) && !is.null(node$result$name.max)) {
        splits <- c(splits, as.character(node$result$name.max))
      }
    }
    
    if (is.list(node)) {
      child_splits <- unlist(
        lapply(node, collect_splits),
        use.names = FALSE
      )
      splits <- c(splits, child_splits)
    }
    
    splits
  }
  
  unique(collect_splits(st))
}
# -----------------------------------------------------------------------

append_results <- function(out, results_path) {
  existing_header <- names(read.csv(results_path, nrows = 0, check.names = FALSE))
  
  # add missing columns as NA
  missing_cols <- setdiff(existing_header, names(out))
  for (nm in missing_cols) {
    out[[nm]] <- NA
  }
  
  # check for unexpected extra columns
  extra_cols <- setdiff(names(out), existing_header)
  if (length(extra_cols) > 0) {
    stop("Output has extra columns not present in results file: ",
         paste(extra_cols, collapse = ", "))
  }
  
  # reorder to match the file
  out <- out[, existing_header, drop = FALSE]
  
  write.table(
    out,
    file = results_path,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
}


#0 or 1 at the root level, not guaranteed total counts across a full tree
