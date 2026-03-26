# analysis_rt.R
# -----------------------------------------------------------------------------
# Purpose
#   Provide analysis-model builders used by your driver:
#     - build_cfa_baseline(): unmoderated CFA for SEM Tree partitioning
#     - build_mnlfa_syntax(): MNLFA-style model with loading moderation via
#       a latent-by-observed interaction .
#

#mxMatrix(type = , #type can be "Diag", "Full" or "Symm"
#         nrow = , #nrow: refers to the number of rows of the matrix.
#         ncol = , #ncol: refers to the number of columns of the matrix.
#        free = , #free: indicates which elements of the matrix can be freely estimated (TRUE or T) or are fixed parameters (FALSE or F).
#         values = , #reflects the values of the elements in the matrix. If an element is freely estimated, it reflects the starting value. If an element is not freely estimated, it reflects the fixed value.
#         name = , #refers to the user-specified name of the matrix which is used within OpenMx when performing an operation on this matrix.
#         ) 

#p = number of manifest variables
#nfactors = number of latent factors

library(lavaan)
library(dplyr)
library(purrr)
library(mxsem)
library(OpenMx)
library(semtree)
library(partykit)

# -----------------------------------------------------------------------------
get_fit_indices <- function(fit) {
  refs <- OpenMx::mxRefModels(fit, run = TRUE)
  summ <- summary(fit, refModels = refs)
  
  list(
    cfi = unname(summ$CFI),
    rmsea = unname(summ$RMSEA)
  )
}

compare_fit_deterioration <- function(fit_less, fit_more,
                                      cfi_cut = 0.01,
                                      rmsea_cut = 0.015) {
  idx_less <- get_fit_indices(fit_less)
  idx_more <- get_fit_indices(fit_more)
  
  delta_cfi <- idx_more$cfi - idx_less$cfi
  delta_rmsea <- idx_more$rmsea - idx_less$rmsea
  
  retain <- !is.na(delta_cfi) &&
    !is.na(delta_rmsea) &&
    delta_cfi >= -cfi_cut &&
    delta_rmsea <= rmsea_cut
  
  list(
    less_cfi = idx_less$cfi,
    less_rmsea = idx_less$rmsea,
    more_cfi = idx_more$cfi,
    more_rmsea = idx_more$rmsea,
    delta_cfi = delta_cfi,
    delta_rmsea = delta_rmsea,
    retain = retain
  )
}
# ---------------------- unrestricted CFA (SEM Tree base) NULL -----------------
mnlfa_analysis <- function(data, p = 4, nfactors = 1) {
  
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
                    name="matB1")
  
  matB2 <- mxMatrix(type="Full", 
                    nrow = 1, 
                    ncol=p,
                    free=TRUE, # shows that effects of m2 in this case are freely estimated
                    values = 0, 
                    name="matB2")
  
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
                    name="matC1")
  
  matC2 <- mxMatrix(type="Full", 
                    nrow=p,
                    ncol = nfactors,
                    free= TRUE, 
                    values = 0,
                    name="matC2")
  
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
  matV1 <- mxMatrix(type="Full", nrow = 1, 
                    ncol = 1,
                    free=FALSE,
                    labels="data.m1",
                    name="m1")
  
  matV2 <- mxMatrix(type="Full", nrow = 1, 
                    ncol = 1,
                    free=FALSE,
                    labels="data.m2",
                    name="m2")
  # -----------------------------------------------------------------------------
  #The mxAlgebra() function can be used to define a matrix of model parameters as
  #a function of background variables. The first argument of this function, 
  #expression, should be used for specifying an R expression of one or more 
  #MxMatrix objects. Most R operators like +, ∗, and %∗%, an general R functions 
  #like mean(), log(), and exp() are supported in this argument of the mxAlgebra() 
  #function. A name for the defined matrix can be assigned with the name argument
  
  # -----------------------------------------------------------------------------
  matT <- mxAlgebra(expression = matT0 +  matB1*m1 + matB2*m2, # linear moderation function for the indicator intercepts
                    name="matT") #intercepts
  
  matL <- mxAlgebra(expression = matL0 + matC1*m1 + matC2*m2, # linear moderation function for the factor loadings
                    name="matL") #loadings
  
  matE <- mxAlgebra(expression = matE0*exp (matD1*m1 + matD2*m2), # log-linear function for the residual variances
                    name="matE") #residual variances
  
  matA <- mxAlgebra(expression = matA0 + matG1*m1 + matG2*m2, #common-factor means modeled as a linear function of the background variables
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
                        matB2,matL, matL0,
                        matC1, matC2, matE, 
                        matE0, matD1, matD2,
                        matP, matP0, 
                        matA, matA0, 
                        matG1, matG2,
                        matV1, matV2,
                        matM, matC, expF, 
                        fitF, mxdata)
  
  #The model can be fitted to the data using the mxRun() function. 
  fitConfig <- mxRun(modConfig)
  
  if (!is.null(fitConfig$output$status$code) && fitConfig$output$status$code == 0) {
    # ---------------------- metric moderated model --------------------------
    
    # Intercept Matrices
    
    matB1 <- mxMatrix(type="Full", 
                      nrow = 1, 
                      ncol=p,#Matrix matB1 and matB2 are full matrices containing the direct effects of the background variables M1 and M2, respectively, on the intercepts.
                      free=TRUE, 
                      values = 0, 
                      name="matB1")
    
    matB2 <- mxMatrix(type="Full", 
                      nrow = 1, 
                      ncol=p,
                      free=TRUE, 
                      values = 0, 
                      name="matB2")
    
    # Loading Matrices
    
    matC1 <- mxMatrix(type="Full", # Matrices matC1 and matC2 are full matrices containing the direct effects of M1 and M2, respectively, on the factor loadings
                      nrow=p,
                      ncol = nfactors,
                      free= FALSE, 
                      values = 0, # value sets the moderation effect
                      name="matC1")
    
    matC2 <- mxMatrix(type="Full", 
                      nrow=p,
                      ncol = nfactors,
                      free= FALSE, 
                      values = 0,
                      name="matC2")
    
    # -----------------------------------------------------------------------------
    matT <- mxAlgebra(expression = matT0 +  matB1*m1 + matB2*m2, # linear moderation function for the indicator intercepts
                      name="matT") #intercepts
    
    matL <- mxAlgebra(expression = matL0 + matC1*m1 + matC2*m2, # linear moderation function for the factor loadings
                      name="matL") #loadings
    
    matE <- mxAlgebra(expression = matE0*exp (matD1*m1 + matD2*m2), # log-linear function for the residual variances
                      name="matE") #residual variances
    
    matA <- mxAlgebra(expression = matA0 + matG1*m1 + matG2*m2, #common-factor means modeled as a linear function of the background variables
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
                          matB2,matL, matL0,
                          matC1, matC2, matE, 
                          matE0, matD1, matD2,
                          matP, matP0, 
                          matA, matA0, 
                          matG1, matG2,
                          matV1, matV2,
                          matM, matC, expF, 
                          fitF, mxdata)
    
    #The model can be fitted to the data using the mxRun() function. 
    fitmetric <- mxRun(modmetric)
    if (is.null(fitmetric$output$status$code) || fitmetric$output$status$code != 0) {
      stop("Metric model did not converge; scalar not run.")
    }
    metric_fit <- compare_fit_deterioration(fitConfig, fitmetric)
    
    metric_retained <- isTRUE(metric_fit$retain)
    
    if (metric_retained) {
      # ---------------------- Scalar moderated model --------------------------
      
      # Intercept Matrices
      
      matB1 <- mxMatrix(type="Full", 
                        nrow = 1, 
                        ncol=p,#Matrix matB1 and matB2 are full matrices containing the direct effects of the background variables M1 and M2, respectively, on the intercepts.
                        free=FALSE, 
                        values = 0, 
                        name="matB1")
      
      matB2 <- mxMatrix(type="Full", 
                        nrow = 1, 
                        ncol=p,
                        free=FALSE, 
                        values = 0, 
                        name="matB2")
      
      # Loading Matrices
      
      matC1 <- mxMatrix(type="Full", # Matrices matC1 and matC2 are full matrices containing the direct effects of M1 and M2, respectively, on the factor loadings
                        nrow=p,
                        ncol = nfactors,
                        free= FALSE, 
                        values = 0, # value sets the moderation effect
                        name="matC1")
      
      matC2 <- mxMatrix(type="Full", 
                        nrow=p,
                        ncol = nfactors,
                        free= FALSE, 
                        values = 0,
                        name="matC2")
      
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
      matT <- mxAlgebra(expression = matT0 +  matB1*m1 + matB2*m2, # linear moderation function for the indicator intercepts
                        name="matT") #intercepts
      
      matL <- mxAlgebra(expression = matL0 + matC1*m1 + matC2*m2, # linear moderation function for the factor loadings
                        name="matL") #loadings
      
      matE <- mxAlgebra(expression = matE0*exp (matD1*m1 + matD2*m2), # log-linear function for the residual variances
                        name="matE") #residual variances
      
      matA <- mxAlgebra(expression = matA0 + matG1*m1 + matG2*m2, #common-factor means modeled as a linear function of the background variables
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
                            matB2,matL, matL0,
                            matC1, matC2, matE, 
                            matE0, matD1, matD2,
                            matP, matP0, 
                            matA, matA0, 
                            matG1, matG2,
                            matV1, matV2,
                            matM, matC, expF, 
                            fitF, mxdata)
      
      #The model can be fitted to the data using the mxRun() function. 
      fitscalar <- mxRun(modscalar)
      if (is.null(fitscalar$output$status$code) || fitscalar$output$status$code != 0) {
        stop("Scalar model did not converge.")
      }
      
      scalar_fit <- compare_fit_deterioration(fitmetric, fitscalar)
      scalar_retained <- isTRUE(scalar_fit$retain)
    }
  }
  return(list(
    fitConfig = fitConfig,
    fitMetric = if (exists("fitmetric")) fitmetric else NULL,,
    fitScalar = if (exists("fitscalar")) fitscalar else NULL,
    metric_fit = if (exists("metric_fit")) metric_fit else NULL,
    scalar_fit = if (exists("scalar_fit")) scalar_fit else NULL
  ))
  }

# -----------------------------------------------------#####################################
tree_analysis_ram <- function(data, p = 4, alpha = 0.05, nfactors = 1,
                              predictors = c("m1", "m2","m0"),
                              control = NULL, verbose = FALSE){
  
  if (is.null(control)) {
    control <- semtree::semtree_control(
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
    focus.parameters = c("lambda_1", "lambda_2", "lambda_3", "lambda_4") #paste0("lambda_", 1:p) #########################flag
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
  
  metric_split <- !(inherits(metric_tree, "error")) &&
    methods::is(metric_tree, "semtree") &&
    length(partykit::nodeids(methods::slot(metric_tree, "tree"), terminal = FALSE)) > 0
  
  scalar_tree <- NULL
  scalar_split <- NA
  
  # ---------------- Scalar stage: intercepts ----------------
  if (!metric_split) {
    scalar_constraints <- semtree::semtree.constraints(
      focus.parameters = c("nu_1", "nu_2", "nu_3", "nu_4") #paste0("nu_", 1:p)
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
    
    scalar_split <- !(inherits(scalar_tree, "error")) &&
      methods::is(scalar_tree, "semtree") &&
      length(partykit::nodeids(methods::slot(scalar_tree, "tree"), terminal = FALSE)) > 0
  }
  
  return(list(
    baseline_model = modbase,
    baseline_fit = fitbase,
    
    metric_tree = metric_tree,
    metric_split = metric_split,
    
    scalar_tree = scalar_tree,
    scalar_split = scalar_split
  ))
}
# ------------------------------------------------------------------------

run_analysis <- function(data,
                         methods = c("MNLFA", "SEMTREE"),
                         nfactors = 1,
                         alpha = 0.05,
                         predictors = c("m1", "m2", "m0")) {
  
  methods <- match.arg(methods, choices = c("MNLFA", "SEMTREE"), several.ok = TRUE)
  dat <- as.data.frame(data)
  
  out <- list(methods = methods)
  
  if ("MNLFA" %in% methods) {
    out$mnlfa <- tryCatch(
      mnlfa_analysis(data = dat, nfactors = nfactors),
      error = identity
    )
  }
  
  if ("SEMTREE" %in% methods) {
    out$semtree <- tryCatch(
      tree_analysis_ram(data = dat,
                    nfactors = nfactors,
                    predictors = predictors),
      error = identity
    )
  }
  
  return(out)
}
# ----------------------------------------------------------------------

semtree_detects_moderation <- function(st, moderators = c("m1", "m2", "m0")) {
  
  out <- list()
  for (m in moderators) {
    out[[paste0("tree_split_on_", m)]] <- NA
    out[[paste0("tree_n_splits_", m)]] <- NA_integer_
  }
  
  if (inherits(st, "error") || is.null(st)) return(out)
  if (!methods::is(st, "semtree")) return(out)
  
  pt <- tryCatch(methods::slot(st, "tree"), error = function(e) NULL)
  if (!inherits(pt, "party")) return(out)
  
  ids <- partykit::nodeids(pt, terminal = FALSE)
  split_vars <- character(0)
  
  if (length(ids) > 0) {
    split_vars <- unlist(partykit::nodeapply(pt, ids, FUN = function(nd) {
      sp <- partykit::split_node(nd)
      if (is.null(sp)) return(NULL)
      names(partykit::data_party(pt))[partykit::varid_split(sp)]
    }))
  }
  
  for (m in moderators) {
    k <- sum(split_vars == m, na.rm = TRUE)
    out[[paste0("tree_split_on_", m)]] <- (k > 0)
    out[[paste0("tree_n_splits_", m)]] <- as.integer(k)
  }
  
  out
}

# -----------------------------------------------------------------------
getPredictorsFromTree <- function(st) {
  if (inherits(st, "error") || is.null(st)) return(NULL)
  if (!methods::is(st, "semtree")) return(NULL)
  
  pt <- tryCatch(methods::slot(st, "tree"), error = function(e) NULL)
  if (!inherits(pt, "party")) return(NULL)
  
  ids <- partykit::nodeids(pt, terminal = FALSE)
  if (length(ids) == 0L) return(character(0))
  
  split_vars <- unlist(partykit::nodeapply(pt, ids, FUN = function(nd) {
    sp <- partykit::split_node(nd)
    if (is.null(sp)) return(NULL)
    names(partykit::data_party(pt))[partykit::varid_split(sp)]
  }))
  
  split_vars
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
