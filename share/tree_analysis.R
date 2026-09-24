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

tree_data <- df_complete %>%
  dplyr::select(
    dplyr::all_of(casp_items),
    age_z,
    female,
    cultural_cluster
  ) %>%
  as.data.frame()

# Checks
dim(tree_data)
anyNA(tree_data)
colSums(is.na(tree_data))


fit_casp_openmx2 <- casp_tree_analysis(
  data = tree_data
)
#------------------------------------------------------------------------------
fit_casp_openmx2$output$status$code
summary(fit_casp_openmx2)
omxGetParameters(fit_casp_openmx2, free = TRUE)
fit_casp_openmx2$output$Minus2LogLikelihood


fit_casp_openmx2$metric_tree
fit_casp_openmx2$metric_test

fit_casp_openmx2$scalar_tree
fit_casp_openmx2$scalar_test

# Metric tree -------------------------------------------------Plot-------------
plot(fit_casp_openmx2$metric_tree)
plotTreeStructure(fit_casp_openmx2$metric_tree)

# labels on the nodes rather than on top 

# Scalar tree
plot(fit_casp_openmx2$scalar_tree)
plotTreeStructure(fit_casp_openmx2$scalar_tree)
#-------------------------------------------------------------------------------
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

metric_items <- c(
  "C2", "C3",
  "A2", "A3",
  "P2", "P3",
  "S2", "S3"
)

metric_constraints <- semtree::semtree.constraints(
  focus.parameters = paste0("lambda_", metric_items)
)

forest_control_test <- semtree::semforest_control(
  num.trees = 5 # TO DO: increase number
)

set.seed(20260915)

metric_forest_test <- semtree::semforest(
  model = fit_casp_openmx2$baseline_fit,
  data = tree_data,
  predictors = c(
    "age_z",
    "female",
    "cultural_cluster"
  ),
  control = forest_control_test,
  constraints = metric_constraints,
  seeds = TRUE
)

#---------------------------------------------
metric_forest_test
summary(metric_forest_test)

metric_forest_test$param.names
#-----------------------------------------------------------------------------

metric_pars_test <- predict(
  metric_forest_test,
  data = tree_data[1:10, , drop = FALSE],
  type = "pars"
)

class(metric_pars_test)
dim(metric_pars_test)
str(metric_pars_test)
metric_pars_test

pd_metric_test <- semtree::partialDependence(
  metric_forest_test,
  data = tree_data,
  reference.var = "age_z",
  support = 20
)

class(pd_metric_test)
dim(pd_metric_test)
str(pd_metric_test)
pd_metric_test
#------------------------------------------------------------------------------
metric_pars <- c(
  "lambda_C2", "lambda_C3",
  "lambda_A2", "lambda_A3",
  "lambda_P2", "lambda_P3",
  "lambda_S2", "lambda_S3"
)

pd <- as.data.frame(pd_metric_test$samples)

par(mfrow = c(2, 4))

for (p in metric_pars) {
  plot(
    pd$age_z,
    pd[[p]],
    type = "b",
    xlab = "Age (z)",
    ylab = "Loading",
    main = p
  )
}

par(mfrow = c(1, 1))
#------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

metric_pars <- c(
  "lambda_C2", "lambda_C3",
  "lambda_A2", "lambda_A3",
  "lambda_P2", "lambda_P3",
  "lambda_S2", "lambda_S3"
)
age_mean <- 69.12062
age_sd   <- 9.737268

# Convert PDP x-axis to age in years
pd_long <- pd_long %>%
  mutate(age = age_mean + age_z * age_sd)

# Get empirical 2.5% and 97.5% limits from age_z,
# then convert those limits to years
age_z_limits <- quantile(
  tree_data$age_z,
  probs = c(.025, .975),
  na.rm = TRUE
)

age_limits <- age_mean + age_z_limits * age_sd

age_limits

pd_long %>%
  filter(
    age >= age_limits[1],
    age <= age_limits[2]
  ) %>%
  count(parameter)

ggplot(
  pd_long %>%
    filter(
      age >= age_limits[1],
      age <= age_limits[2]
    ),
  aes(x = age, y = loading)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 4) +
  labs(
    x = "Age (years)",
    y = "Forest-predicted loading"
  ) +
  theme_minimal()
#-------------------------------------------------------------------------------
# One-tree PDP
#set.seed(20260916)

#t_fast <- system.time({
  
#  metric_forest_fast1 <- semtree::semforest(
#    model = fit_casp_openmx2$baseline_fit,
#    data = tree_data,
#    predictors = c(
#      "age_z",
#      "female",
#      "cultural_cluster"
#    ),
#    control = forest_control_fast_test,
#    constraints = metric_constraints,
#    seeds = TRUE
#  )
#  
#})
#
#t_fast
#
#metric_forest_fast1
#summary(metric_forest_fast1)
#
#plot(metric_forest_fast1$forest[[1]])
#
#
#pd_metric_fast1 <- semtree::partialDependence(
#  metric_forest_fast1,
#  data = tree_data,
#  reference.var = "age_z",
#  support = 20
#)
#
#pd_fast <- as.data.frame(pd_metric_fast1$samples)
#----------------------------------------------------------------------------

pd_long_fast <- pd_fast %>%
  select(age_z, all_of(metric_pars)) %>%
  pivot_longer(
    cols = all_of(metric_pars),
    names_to = "parameter",
    values_to = "loading"
  )

age_z_limits <- quantile(
  tree_data$age_z,
  probs = c(.025, .975),
  na.rm = TRUE
)

age_limits <- age_mean + age_z_limits * age_sd

############################################################# 22/09
# Constrained metric SEM forest:
# depth 3, 10 trees
#########################################

# Keep individual trees shallow
forest_tree_control <- semtree::semtree_control(
  method = "naive",
  alpha = 1,
  bonferroni = FALSE,
  max.depth = 3,
  min.N = 500,
  exclude.heywood = FALSE
)

# Increase number of trees for greater
# ensemble averaging / PDP stability
forest_control_10 <- semtree::semforest_control(
  num.trees = 10,
  sampling = "subsample",
  control = forest_tree_control,
  mtry = 2
)

forest_control_10


#########################################
# Fit metric forest
#########################################

set.seed(20260916)

t_metric_10 <- system.time({
  
  metric_forest_10 <- semtree::semforest(
    model = fit_casp_openmx2$baseline_fit,
    data = tree_data,
    predictors = c(
      "age_z",
      "female",
      "cultural_cluster"
    ),
    control = forest_control_10,
    constraints = metric_constraints,
    seeds = TRUE
  )
  
})

t_metric_10

metric_forest_10
summary(metric_forest_10)

# Save immediately
saveRDS(
  metric_forest_10,
  "metric_forest_depth3_min500_10trees.rds"
)

#########################################
# Partial dependence for age
#########################################

pd_metric_10 <- semtree::partialDependence(
  metric_forest_10,
  data = tree_data,
  reference.var = "age_z",
  support = 20
)

saveRDS(
  pd_metric_10,
  "pd_metric_depth3_min500_10trees.rds"
)


#########################################
# Prepare PDP for plotting
#########################################

metric_pars <- c(
  "lambda_C2", "lambda_C3",
  "lambda_A2", "lambda_A3",
  "lambda_P2", "lambda_P3",
  "lambda_S2", "lambda_S3"
)

pd_10 <- as.data.frame(pd_metric_10$samples)

names(pd_10)

pd_long_10 <- pd_10 %>%
  select(age_z, all_of(metric_pars)) %>%
  pivot_longer(
    cols = all_of(metric_pars),
    names_to = "parameter",
    values_to = "loading"
  ) %>%
  mutate(
    age = age_mean + age_z * age_sd
  )


#########################################
# Restrict displayed age range
#########################################

age_z_limits <- quantile(
  tree_data$age_z,
  probs = c(.025, .975),
  na.rm = TRUE
)

age_limits <- age_mean + age_z_limits * age_sd


#########################################
# Plot
#########################################

ggplot(
  pd_long_10 %>%
    filter(
      age >= age_limits[1],
      age <= age_limits[2]
    ),
  aes(x = age, y = loading)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(
    ~ parameter,
    scales = "free_y",
    ncol = 4
  ) +
  labs(
    x = "Age (years)",
    y = "Forest-predicted loading"
  ) +
  theme_minimal()

#--------------------------------------------
pd_plot <- pd_long_10 %>%
  filter(
    age >= age_limits[1],
    age <= age_limits[2]
  ) %>%
  mutate(
    factor = case_when(
      grepl("^lambda_C", parameter) ~ "C",
      grepl("^lambda_A", parameter) ~ "A",
      grepl("^lambda_P", parameter) ~ "P",
      grepl("^lambda_S", parameter) ~ "S"
    ),
    indicator = case_when(
      grepl("2$", parameter) ~ "Indicator 2",
      grepl("3$", parameter) ~ "Indicator 3"
    )
  )



ggplot(
  pd_plot,
  aes(
    x = age,
    y = loading,
    linetype = indicator,
    group = parameter
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  facet_wrap(
    ~ factor,
    ncol = 2,
    scales = "fixed",
    axes = "all",
    axis.labels = "all"
  ) +
  scale_x_continuous(
    breaks = seq(50, 90, by = 10)
  ) +
  scale_y_continuous(
    breaks = seq(0.4, 1.8, by = 0.2),
    labels = function(x) {
      ifelse(
        seq_along(x) %% 2 == 1,
        sprintf("%.1f", x),
        ""
      )
    }
  ) +
  labs(
    x = "Age (years)",
    y = "Forest-predicted factor loading",
    linetype = NULL
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    strip.background = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.65)
  )
