library(dplyr)
library(tidyr)
library(stringr)
library(tibble)


mnlfa_mods <- c(
  "age_z",
  "female",
  "culture_nordic",
  "culture_southern",
  "culture_eastern"
)

mnlfa_data <- df_complete |>
  dplyr::select(
    dplyr::all_of(casp_items),
    age_z,
    female,
    culture_nordic,
    culture_southern,
    culture_eastern
  ) |>
  as.data.frame()
str(mnlfa_data)
#--------------------------------------------
casp_mnlfa_analysis <- function(data, alpha = .05) {
  
  dat <- as.data.frame(data)
  
  manVars <- c(
    "C1","C2","C3",
    "A1","A2","A3",
    "P1","P2","P3",
    "S1","S2","S3"
  )
  
  latVars <- c(
    "Control",
    "Autonomy",
    "Pleasure",
    "SelfRealization"
  )
  
  moderators <- c(
    "age_z",
    "female",
    "culture_nordic",
    "culture_southern",
    "culture_eastern"
  )
  
  anchor_items <- c("C1","A1","P1","S1")
  
  nonanchor_items <- c(
    "C2","C3",
    "A2","A3",
    "P2","P3",
    "S2","S3"
  )
  
  
  # -----------------------------------------------------------------------
  # Loading structure
  # -----------------------------------------------------------------------
  
  loading_free <- matrix(
    FALSE,
    nrow = 12,
    ncol = 4,
    dimnames = list(manVars, latVars)
  )
  
  loading_free["C2","Control"] <- TRUE
  loading_free["C3","Control"] <- TRUE
  
  loading_free["A2","Autonomy"] <- TRUE
  loading_free["A3","Autonomy"] <- TRUE
  
  loading_free["P2","Pleasure"] <- TRUE
  loading_free["P3","Pleasure"] <- TRUE
  
  loading_free["S2","SelfRealization"] <- TRUE
  loading_free["S3","SelfRealization"] <- TRUE
  
  
  loading_values <- matrix(
    0,
    nrow = 12,
    ncol = 4,
    dimnames = list(manVars, latVars)
  )
  
  # marker loadings
  loading_values["C1","Control"] <- 1
  loading_values["A1","Autonomy"] <- 1
  loading_values["P1","Pleasure"] <- 1
  loading_values["S1","SelfRealization"] <- 1
  
  # starting values from fitted marker CFA
  loading_values["C2","Control"] <- 1.084
  loading_values["C3","Control"] <- .987
  
  loading_values["A2","Autonomy"] <- .507
  loading_values["A3","Autonomy"] <- 1.078
  
  loading_values["P2","Pleasure"] <- 1.077
  loading_values["P3","Pleasure"] <- .862
  
  loading_values["S2","SelfRealization"] <- 1.140
  loading_values["S3","SelfRealization"] <- 1.168
  
  
  # -----------------------------------------------------------------------
  # Function constructing one MNLFA model
  # -----------------------------------------------------------------------
  
  build_mnlfa <- function(
    model_name,
    loading_mod = TRUE,
    intercept_mod = TRUE
  ) {
    
    mxdata <- mxData(
      observed = dat,
      type = "raw"
    )
    
    
    # =====================================================================
    # Baseline item intercepts
    # =====================================================================
    
    intercept_free <- !manVars %in% anchor_items
    
    matT0 <- mxMatrix(
      type = "Full",
      nrow = 1,
      ncol = 12,
      free = intercept_free,
      values = c(
        0, .098, .763,
        0, 1.555, -.706,
        0, -.236, .384,
        0, -.432, -.582
      ),
      labels = ifelse(
        intercept_free,
        paste0("nu_", manVars),
        NA
      ),
      name = "matT0"
    )
    
    
    # =====================================================================
    # Intercept moderation
    # =====================================================================
    
    B_free <- rep(intercept_mod, 12)
    B_free[manVars %in% anchor_items] <- FALSE
    
    matB_age <- mxMatrix(
      "Full", 1, 12,
      free = B_free,
      values = 0,
      labels = ifelse(B_free, paste0("B_age_", manVars), NA),
      name = "matB_age"
    )
    
    matB_female <- mxMatrix(
      "Full", 1, 12,
      free = B_free,
      values = 0,
      labels = ifelse(B_free, paste0("B_female_", manVars), NA),
      name = "matB_female"
    )
    
    matB_nordic <- mxMatrix(
      "Full", 1, 12,
      free = B_free,
      values = 0,
      labels = ifelse(B_free, paste0("B_nordic_", manVars), NA),
      name = "matB_nordic"
    )
    
    matB_southern <- mxMatrix(
      "Full", 1, 12,
      free = B_free,
      values = 0,
      labels = ifelse(B_free, paste0("B_southern_", manVars), NA),
      name = "matB_southern"
    )
    
    matB_eastern <- mxMatrix(
      "Full", 1, 12,
      free = B_free,
      values = 0,
      labels = ifelse(B_free, paste0("B_eastern_", manVars), NA),
      name = "matB_eastern"
    )
    
    
    # =====================================================================
    # Baseline loadings
    # =====================================================================
    
    matL0 <- mxMatrix(
      type = "Full",
      nrow = 12,
      ncol = 4,
      free = loading_free,
      values = loading_values,
      name = "matL0"
    )
    
    
    # =====================================================================
    # Loading moderation
    # =====================================================================
    
    C_free <- loading_free
    
    if (!loading_mod) {
      C_free[,] <- FALSE
    }
    
    make_C_labels <- function(prefix) {
      
      labs <- matrix(
        NA_character_,
        12,
        4,
        dimnames = list(manVars, latVars)
      )
      
      for (i in manVars) {
        for (f in latVars) {
          if (C_free[i, f]) {
            labs[i, f] <- paste0(prefix, "_", i)
          }
        }
      }
      
      labs
    }
    
    matC_age <- mxMatrix(
      "Full", 12, 4,
      free = C_free,
      values = 0,
      labels = make_C_labels("C_age"),
      name = "matC_age"
    )
    
    matC_female <- mxMatrix(
      "Full", 12, 4,
      free = C_free,
      values = 0,
      labels = make_C_labels("C_female"),
      name = "matC_female"
    )
    
    matC_nordic <- mxMatrix(
      "Full", 12, 4,
      free = C_free,
      values = 0,
      labels = make_C_labels("C_nordic"),
      name = "matC_nordic"
    )
    
    matC_southern <- mxMatrix(
      "Full", 12, 4,
      free = C_free,
      values = 0,
      labels = make_C_labels("C_southern"),
      name = "matC_southern"
    )
    
    matC_eastern <- mxMatrix(
      "Full", 12, 4,
      free = C_free,
      values = 0,
      labels = make_C_labels("C_eastern"),
      name = "matC_eastern"
    )
    
    
    # =====================================================================
    # Residual variances
    # =====================================================================
    
    matE0 <- mxMatrix(
      type = "Diag",
      nrow = 12,
      ncol = 12,
      free = TRUE,
      values = c(
        .640, .483, .445,
        .642, .812, .995,
        .303, .209, .341,
        .345, .269, .280
      ),
      name = "matE0"
    )
    
    
    # Residual variance moderation, as in analysis_rt
    
    matD_age <- mxMatrix(
      "Diag", 12, 12,
      free = TRUE,
      values = 0,
      name = "matD_age"
    )
    
    matD_female <- mxMatrix(
      "Diag", 12, 12,
      free = TRUE,
      values = 0,
      name = "matD_female"
    )
    
    matD_nordic <- mxMatrix(
      "Diag", 12, 12,
      free = TRUE,
      values = 0,
      name = "matD_nordic"
    )
    
    matD_southern <- mxMatrix(
      "Diag", 12, 12,
      free = TRUE,
      values = 0,
      name = "matD_southern"
    )
    
    matD_eastern <- mxMatrix(
      "Diag", 12, 12,
      free = TRUE,
      values = 0,
      name = "matD_eastern"
    )
    
    
    # =====================================================================
    # Latent covariance matrix
    # =====================================================================
    
    matP0 <- mxMatrix(
      type = "Symm",
      nrow = 4,
      ncol = 4,
      free = TRUE,
      values = matrix(
        c(
          .394, .216, .160, .235,
          .216, .162, .152, .193,
          .160, .152, .266, .248,
          .235, .193, .248, .366
        ),
        4, 4,
        byrow = TRUE
      ),
      name = "matP0"
    )
    
    
    # =====================================================================
    # Latent means
    # =====================================================================
    
    matA0 <- mxMatrix(
      type = "Full",
      nrow = 1,
      ncol = 4,
      free = TRUE,
      values = c(
        2.526,
        3.225,
        3.505,
        3.088
      ),
      name = "matA0"
    )
    
    
    # IMPORTANT:
    # latent mean moderation remains free in ALL models
    
    matG_age <- mxMatrix(
      "Full", 1, 4,
      free = TRUE,
      values = 0,
      name = "matG_age"
    )
    
    matG_female <- mxMatrix(
      "Full", 1, 4,
      free = TRUE,
      values = 0,
      name = "matG_female"
    )
    
    matG_nordic <- mxMatrix(
      "Full", 1, 4,
      free = TRUE,
      values = 0,
      name = "matG_nordic"
    )
    
    matG_southern <- mxMatrix(
      "Full", 1, 4,
      free = TRUE,
      values = 0,
      name = "matG_southern"
    )
    
    matG_eastern <- mxMatrix(
      "Full", 1, 4,
      free = TRUE,
      values = 0,
      name = "matG_eastern"
    )
    
    
    # =====================================================================
    # Definition variables
    # =====================================================================
    
    age <- mxMatrix(
      "Full", 1, 1,
      free = FALSE,
      labels = "data.age_z",
      name = "age"
    )
    
    female <- mxMatrix(
      "Full", 1, 1,
      free = FALSE,
      labels = "data.female",
      name = "female"
    )
    
    nordic <- mxMatrix(
      "Full", 1, 1,
      free = FALSE,
      labels = "data.culture_nordic",
      name = "nordic"
    )
    
    southern <- mxMatrix(
      "Full", 1, 1,
      free = FALSE,
      labels = "data.culture_southern",
      name = "southern"
    )
    
    eastern <- mxMatrix(
      "Full", 1, 1,
      free = FALSE,
      labels = "data.culture_eastern",
      name = "eastern"
    )
    
    
    # =====================================================================
    # MNLFA algebras
    # =====================================================================
    
    matT <- mxAlgebra(
      matT0 +
        matB_age * age +
        matB_female * female +
        matB_nordic * nordic +
        matB_southern * southern +
        matB_eastern * eastern,
      name = "matT"
    )
    
    
    matL <- mxAlgebra(
      matL0 +
        matC_age * age +
        matC_female * female +
        matC_nordic * nordic +
        matC_southern * southern +
        matC_eastern * eastern,
      name = "matL"
    )
    
    
    matE <- mxAlgebra(
      matE0 * exp(
        matD_age * age +
          matD_female * female +
          matD_nordic * nordic +
          matD_southern * southern +
          matD_eastern * eastern
      ),
      name = "matE"
    )
    
    
    matA <- mxAlgebra(
      matA0 +
        matG_age * age +
        matG_female * female +
        matG_nordic * nordic +
        matG_southern * southern +
        matG_eastern * eastern,
      name = "matA"
    )
    
    
    # latent covariance matrix currently allowed to be freely estimated,
    # but constant across moderators
    matP <- mxAlgebra(
      matP0,
      name = "matP"
    )
    
    
    # =====================================================================
    # Implied moments
    # =====================================================================
    
    matM <- mxAlgebra(
      matT + matA %*% t(matL),
      name = "matM"
    )
    
    matCOV <- mxAlgebra(
      matL %*% matP %*% t(matL) + matE,
      name = "matCOV"
    )
    
    
    expF <- mxExpectationNormal(
      covariance = "matCOV",
      means = "matM",
      dimnames = manVars
    )
    
    fitF <- mxFitFunctionML()
    
    
    mxModel(
      model_name,
      
      matT0,
      matB_age,
      matB_female,
      matB_nordic,
      matB_southern,
      matB_eastern,
      
      matL0,
      matC_age,
      matC_female,
      matC_nordic,
      matC_southern,
      matC_eastern,
      
      matE0,
      matD_age,
      matD_female,
      matD_nordic,
      matD_southern,
      matD_eastern,
      
      matP0,
      
      matA0,
      matG_age,
      matG_female,
      matG_nordic,
      matG_southern,
      matG_eastern,
      
      age,
      female,
      nordic,
      southern,
      eastern,
      
      matT,
      matL,
      matE,
      matA,
      matP,
      matM,
      matCOV,
      
      expF,
      fitF,
      mxdata
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Configural / unrestricted MNLFA
  # -----------------------------------------------------------------------
  
  modConfig <- build_mnlfa(
    "CASP_Configural",
    loading_mod = TRUE,
    intercept_mod = TRUE
  )
  
  fitConfig <- mxRun(modConfig)
  
  
  if (fitConfig$output$status$code != 0) {
    stop(
      "Configural MNLFA failed. Status: ",
      fitConfig$output$status$code
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Metric model
  #
  # loading moderation fixed to zero
  # intercept moderation still free
  # -----------------------------------------------------------------------
  
  modMetric <- build_mnlfa(
    "CASP_Metric",
    loading_mod = FALSE,
    intercept_mod = TRUE
  )
  
  fitMetric <- mxRun(modMetric)
  
  
  if (fitMetric$output$status$code != 0) {
    stop(
      "Metric MNLFA failed. Status: ",
      fitMetric$output$status$code
    )
  }
  
  
  metric_lrt <- compare_models_lrt(
    fitConfig,
    fitMetric,
    alpha = alpha
  )
  
  
  # -----------------------------------------------------------------------
  # Scalar model
  #
  # loading moderation = 0
  # intercept moderation = 0
  #
  # latent-mean moderation remains free
  # -----------------------------------------------------------------------
  
  modScalar <- build_mnlfa(
    "CASP_Scalar",
    loading_mod = FALSE,
    intercept_mod = FALSE
  )
  
  fitScalar <- mxRun(modScalar)
  
  
  if (fitScalar$output$status$code != 0) {
    stop(
      "Scalar MNLFA failed. Status: ",
      fitScalar$output$status$code
    )
  }
  
  
  scalar_lrt <- compare_models_lrt(
    fitMetric,
    fitScalar,
    alpha = alpha
  )
  
  
  return(
    list(
      fitConfig = fitConfig,
      fitMetric = fitMetric,
      fitScalar = fitScalar,
      metric_lrt = metric_lrt,
      scalar_lrt = scalar_lrt
    )
  )
}
#------------------------------------------
casp_mnlfa <- casp_mnlfa_analysis(
  data = mnlfa_data
)

#----------------------------------------
casp_mnlfa$fitConfig$output$status$code
casp_mnlfa$fitMetric$output$status$code
casp_mnlfa$fitScalar$output$status$code

casp_mnlfa$metric_lrt
casp_mnlfa$scalar_lrt
#----------------------------------------

# ------------------------------------------------------------
# MNLFA model-fit information
# ------------------------------------------------------------

get_mnlfa_fit <- function(fit, model_name) {
  
  npar <- length(
    OpenMx::omxGetParameters(
      fit,
      free = TRUE
    )
  )
  
  minus2LL <- fit$output$fit
  
  tibble::tibble(
    Model = model_name,
    npar = npar,
    minus2LL = minus2LL,
    AIC = minus2LL + 2 * npar,
    BIC = minus2LL + log(nrow(mnlfa_data)) * npar
  )
}

mnlfa_fit_summary <- dplyr::bind_rows(
  get_mnlfa_fit(
    casp_mnlfa$fitConfig,
    "Unrestricted"
  ),
  get_mnlfa_fit(
    casp_mnlfa$fitMetric,
    "Metric"
  ),
  get_mnlfa_fit(
    casp_mnlfa$fitScalar,
    "Scalar"
  )
)

mnlfa_fit_summary
#------------------------------------------------------------------------------
mnlfa_fit_table <- mnlfa_fit_summary %>%
  mutate(
    minus2LL = round(minus2LL, 2),
    AIC = round(AIC, 2),
    BIC = round(BIC, 2)
  )

mnlfa_fit_table

knitr::kable(
  mnlfa_fit_table,
  caption = "Fit statistics for the empirical MNLFA models."
)
#------------------------------------------------------------------------------
mnlfa_pars <- OpenMx::omxGetParameters(
  casp_mnlfa$fitConfig,
  free = TRUE
)

# Loading moderation
loading_mod <- mnlfa_pars[
  grepl("^C_(age|female|nordic|southern|eastern)_", names(mnlfa_pars))
]

# Intercept moderation
intercept_mod <- mnlfa_pars[
  grepl("^B_(age|female|nordic|southern|eastern)_", names(mnlfa_pars))
]

loading_mod
intercept_mod

#--------------------------------------------------------------------------
#==========================================================================

mnlfa_lrt_summary <- bind_rows(
  tibble(
    Test = "Metric invariance",
    Comparison = "Unrestricted vs. metric",
    Chisq = casp_mnlfa$metric_lrt$chisq_diff,
    df = casp_mnlfa$metric_lrt$df_diff,
    p = casp_mnlfa$metric_lrt$p_value,
    Reject = casp_mnlfa$metric_lrt$reject_h0
  ),
  tibble(
    Test = "Scalar invariance",
    Comparison = "Metric vs. scalar",
    Chisq = casp_mnlfa$scalar_lrt$chisq_diff,
    df = casp_mnlfa$scalar_lrt$df_diff,
    p = casp_mnlfa$scalar_lrt$p_value,
    Reject = casp_mnlfa$scalar_lrt$reject_h0
  )
) %>%
  mutate(
    p_report = ifelse(
      is.na(p),
      NA_character_,
      ifelse(p < .001, "< .001", sprintf("%.3f", p))
    )
  )

mnlfa_lrt_summary

mnlfa_lrt_table <- mnlfa_lrt_summary %>%
  transmute(
    Test,
    Comparison,
    `Chi-square` = round(Chisq, 2),
    df,
    p = p_report
  )

mnlfa_lrt_table
#----------------------------------------------------------------------------

loading_summary <- enframe(
  loading_mod,
  name = "parameter",
  value = "estimate"
) %>%
  mutate(
    type = "Loading",
    parameter = sub("^C_", "", parameter)
  ) %>%
  separate(
    parameter,
    into = c("moderator", "item"),
    sep = "_",
    extra = "merge"
  )

intercept_summary <- enframe(
  intercept_mod,
  name = "parameter",
  value = "estimate"
) %>%
  mutate(
    type = "Intercept",
    parameter = sub("^B_", "", parameter)
  ) %>%
  separate(
    parameter,
    into = c("moderator", "item"),
    sep = "_",
    extra = "merge"
  )

mnlfa_coef_summary <- bind_rows(
  loading_summary,
  intercept_summary
) %>%
  mutate(
    moderator = recode(
      moderator,
      age = "Age",
      female = "Female",
      nordic = "Nordic",
      southern = "Southern",
      eastern = "Eastern"
    ),
    abs_estimate = abs(estimate)
  ) %>%
  arrange(type, moderator, item)

mnlfa_coef_summary

mnlfa_loading_table <- mnlfa_coef_summary %>%
  filter(type == "Loading") %>%
  select(item, moderator, estimate) %>%
  pivot_wider(
    names_from = moderator,
    values_from = estimate
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

mnlfa_intercept_table <- mnlfa_coef_summary %>%
  filter(type == "Intercept") %>%
  select(item, moderator, estimate) %>%
  pivot_wider(
    names_from = moderator,
    values_from = estimate
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

mnlfa_loading_table
mnlfa_intercept_table
#------------------------------------------------
knitr::kable(
  mnlfa_lrt_table,
  caption = "Likelihood-ratio tests of metric and scalar measurement invariance in the MNLFA."
)

knitr::kable(
  mnlfa_loading_table,
  digits = 3,
  caption = "MNLFA loading-moderation coefficients."
)

knitr::kable(
  mnlfa_intercept_table,
  digits = 3,
  caption = "MNLFA intercept-moderation coefficients."
)

#---------------------
mnlfa_largest_effects <- mnlfa_coef_summary %>%
  group_by(type, moderator) %>%
  slice_max(
    order_by = abs_estimate,
    n = 2,
    with_ties = FALSE
  ) %>%
  arrange(type, moderator, desc(abs_estimate))

mnlfa_largest_effects
