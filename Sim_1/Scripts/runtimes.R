library(dplyr)

runtimes_summary2 <- evaluation_runtimes |>
  dplyr::group_by(
    N,
    moderator_1_type,
    variant_label,
    method
  ) |>
  dplyr::summarise(
    mean_rt = mean(runtime_seconds, na.rm = TRUE),
    .groups = "drop"
  )

## 2) Recover number of replications per cell from the original table
reps_per_cell <- nrow(evaluation_runtimes) / nrow(runtimes_summary2)
reps_per_cell


## 3) Total runtime of the *current* simulation design (all cells, all reps)
total_runtime_current <- runtimes_summary2 |>
  summarise(
    total_seconds = sum(mean_rt * reps_per_cell),
    total_minutes = total_seconds / 60,
    total_hours   = total_minutes / 60
  )

total_runtime_current

# --------------------------------------------------

runtime_df <- runtimes_summary2 |>
  mutate(
    # map variant -> p
    p = case_when(
      variant_label == "MNLFA_none"           ~ 12,
      variant_label == "MNLFA_linear_partial" ~ 14,
      variant_label == "MNLFA_linear_full"    ~ 19,
      TRUE ~ NA_real_
    ),
    logN = log(N),
    logp = log(p)
  )

fit_runtime <- lm(
  log(mean_rt) ~ logN + logp + method + variant_label,
  data = runtime_df
)

summary(fit_runtime)

# ----------------------------------------------------------
# Approximation of final runtime
# ----------------------------------------------------------

future_cells <- tidyr::expand_grid(
  N            = c(250, 500, 700, 1000, 5000),              
  moderator_1_type = c("linear", "sigmoid", "quadratic"),     
  variant_label    = c("MNLFA_none", "MNLFA_linear_partial", "MNLFA_linear_full"),
  method           = c("mxsem", "semtree"),
  reps             = 1000                          
)



# ----------------------------------------------------------
future_cells <- future_cells |>
  mutate(
    p = case_when(
      variant_label == "MNLFA_none"           ~ 12,
      variant_label == "MNLFA_linear_partial" ~ 14,
      variant_label == "MNLFA_linear_full"    ~ 19,
      TRUE ~ NA_real_
    ),
    logN = log(N),
    logp = log(p),
    pred_log_rt = predict(fit_runtime, newdata = cur_data()),
    pred_rt      = exp(pred_log_rt),        # seconds per fit
    cell_seconds = pred_rt * reps,
    cell_minutes = cell_seconds / 60,
    cell_hours   = cell_seconds / 3600
  )

# total predicted runtime for the design
future_total <- future_cells |>
  summarise(
    total_seconds = sum(cell_seconds),
    total_minutes = total_seconds / 60,
    total_hours   = total_seconds / 3600
  )

