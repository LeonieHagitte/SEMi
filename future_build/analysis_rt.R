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
# ---------------------- unrestricted CFA (SEM Tree base) NULL -----------------
mnlfa_analysis <- function(data, p = 4, alpha=0.05, nfactors = 1) {
  
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
  matM <- mxAlgebra(matT + t(matL %*% matA), name="matM")
  
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
  
  ##################################################################################
  print(fitConfig$output$status)
  str(fitConfig$output$status)
  #####################################################################################
  
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
    matM <- mxAlgebra(matT + t(matL %*% matA), name="matM")
    
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
    
    miTest1 <- mxCompare(fitConfig, fitmetric)
    p_metric <- miTest1$p[2]  # safest for two-model compare
    
    if (!is.na(p_metric) && p_metric >= alpha) {
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
      matM <- mxAlgebra(matT + t(matL %*% matA), name="matM")
      
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
      miTest2 <- mxCompare(fitmetric, fitscalar)
    }
  }
  return(list(
    fitConfig = fitConfig,
    fitMetric = fitmetric,
    fitScalar = if (exists("fitscalar")) fitscalar else NULL,
    miTest_metric = if (exists("miTest1")) miTest1 else NULL,
    miTest_scalar = if (exists("miTest2")) miTest2 else NULL
  ))
  }
# ------------------------------------------------------------------------

tree_analysis <- function(data, p = 4, alpha=0.05, nfactors = 1, predictors = c("m1","m2","m0"),
                          control = semtree::semtree_control(method = "score"), verbose = FALSE){
  
  stopifnot(nfactors==1) #ToDO: latent covariance matrix is not populated correctly
  
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
  
  mxdata <- mxData(observed = dat, type = "raw")
  
  # Intercept Matrices
  matT0 <- mxMatrix(type = "Full", 
                    nrow = 1, 
                    ncol = p, 
                    free = TRUE, 
                    values = 0.6,
                    labels = paste0("nu_", 1:p),
                    name = "matT0") #The matrix matT0 is a full matrix containing the baseline intercepts.
  
  # Loading Matrices
  matL0 <- mxMatrix(type="Full", #matL0 is a full matrix containing the baseline factor loadings
                    nrow=p, 
                    ncol = nfactors, 
                    free= TRUE, 
                    values= rep(1, p),
                    labels = paste0("lambda_", 1:p),
                    name="matL0")
  
  # Residual Variances 
  matE0 <- mxMatrix(type="Diag", # matE0 is a diagonal matrix containing the baseline residual variances 
                    nrow=p, 
                    ncol=p,
                    free=TRUE,
                    values = 1,
                    name="matE0")
  
  # Latent factor variance
  matP0 <- mxMatrix(type="Symm", 
                    nrow = nfactors, #because only one latent factor 
                    ncol = nfactors,
                    free= FALSE,
                    values= 1, #todo - assumes nfactors is 1
                    name="matP0")
  
  # latent factor mean
  matA0 <- mxMatrix(type="Full", #matA0 is a matrix containing the baseline common-factor means 
                    nrow = 1, 
                    ncol = nfactors,
                    free=FALSE, #fixed
                    values = 0, # to 0
                    name="matA0")
  
  # -----------------------------------------------------------------------------
  matT <- mxAlgebra(expression = matT0, name="matT") #intercepts baseline
  
  matL <- mxAlgebra(expression = matL0, name="matL") #loadings baseline
  
  matE <- mxAlgebra(expression = matE0, name="matE") #residual variances
  
  matA <- mxAlgebra(expression = matA0, name="matA") #matrix of common-factor means
  
  matP <- mxAlgebra(matP0, name="matP") #latent covariance matrix in the covariance algebra - for a single-factor CFA, that matrix is simply 1×1 and is the factor variance
  
  # ------------ Matrices for model implied moments -------
  matM <- mxAlgebra(matT + t(matL %*% matA), name="matM")
  
  matC <- mxAlgebra(expression = matL %*% matP%*%t(matL) + matE ,
                    name="matC")
  
  # ------------ Model expectations and fit-function -------
  expTr <- mxExpectationNormal(covariance="matC", means = "matM", dimnames=manVars)
  
  fitTr <- mxFitFunctionML() #mxFitFunctionML() function stored in fitF is used to 
  #indicate that the free parameters of the configural 
  #model should be estimated using full-information maximum likelihood.
  
  modbase <- mxModel(model="baseline",
                       matT, matT0, 
                       matL, matL0,
                       matE, matE0, 
                       matP, matP0, 
                       matA, matA0, 
                       matM, matC, expTr, 
                       fitTr, mxdata)
  
  #The model can be fitted to the data using the mxRun() function. 
  fitbase <- mxRun(modbase)
  
  # ---------------- Metric stage: loadings ----------------
  metric_constraints <- semtree::semtree.constraints(
    focus.parameters = paste0("lambda_", 1:p)
  )
  
  metric_tree <- tryCatch(
    semtree::semtree(
      model = modbase,
      data = dat,
      predictors = predictors,
      control = control,
      constraints = metric_constraints,
      verbose = verbose
    ),
    error = identity
  )
  
  metric_split <- !(inherits(metric_tree, "error")) &&
    length(partykit::nodeids(metric_tree, terminal = FALSE)) > 0
  
  scalar_tree <- NULL
  scalar_split <- NA
  
  # ---------------- Scalar stage: intercepts ----------------
  if (!metric_split) {
    scalar_constraints <- semtree::semtree.constraints(
      focus.parameters = paste0("nu_", 1:p)
    )
    
    scalar_tree <- tryCatch(
      semtree::semtree(
        model = modbase,
        data = dat,
        predictors = predictors,
        control = control,
        constraints = scalar_constraints,
        verbose = verbose
      ),
      error = identity
    )
    
    scalar_split <- !(inherits(scalar_tree, "error")) &&
      length(partykit::nodeids(scalar_tree, terminal = FALSE)) > 0
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
      mnlfa_analysis(data = dat, nfactors = nfactors, alpha = alpha),
      error = identity
    )
  }
  
  if ("SEMTREE" %in% methods) {
    out$semtree <- tryCatch(
      tree_analysis(data = dat,
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

