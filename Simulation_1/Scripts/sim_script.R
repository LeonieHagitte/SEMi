# ---- sim_script.R ----
# Purpose: main simulation loop (optionally parallel), saving raw per-condition outputs

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(doParallel)
  library(foreach)
  library(doRNG)
  library(lavaan)
  library(semtree)
})

# source functions
source(here("Simulation_1", "Functions", "dataprep_functions.R"))
source(here("Simulation_1", "Functions", "sim_functions.R"))

results_dir <- here("Simulation_1", "Results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# reproducibility
set.seed(20251)

# load conditions
conditions <- readRDS(file.path(results_dir, "conditions.rds"))

# repetitions per condition
n_reps <- 300

# parallel setup
n_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# run simulations over conditions
start_time <- Sys.time()

results_list <- foreach(
  i = 1:nrow(conditions),
  .packages = c("tidyverse", "mxsem", "glue", "lavaan", "semtree"),
  
  .export   = c("generate_moderator", "generate_data_mni", "build_null_model",
                               "build_model_moderated_loadings", "build_model_moderated_intercepts",
                               "build_model_moderated_residuals", "estimate_mxsem", "run_condition",
                               "choose_target_param_name", "true_value_for_target",
                               "estimate_semtree", "build_null_model_lavaan")
) %dorng% {
  cond <- conditions[i, ]
  mvec <- generate_moderator(cond$n, cond$moderator_type)
  
  out <- run_condition(
    n = cond$n,
    moderator_vec = mvec,
    mod_items_load = cond$mod_items_load[[1]],
    mod_items_int  = cond$mod_items_int[[1]],
    mod_items_res  = cond$mod_items_res[[1]],
    slope_load = cond$slope_load,
    slope_int  = cond$slope_int,
    slope_logres = cond$slope_logres,
    n_reps = n_reps,
    null_model_string = build_null_model(),
    alt_model_string  = cond$alt_model_string[[1]],
    log_link_resid_in_analysis = TRUE,
    moderator_type = cond$moderator_type   # <-- PASS THROUGH
  )
  
  out
}

stopCluster(cl)

simtime_total <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

# Split the results into summaries and rep-level records
summaries_list <- lapply(results_list, `[[`, "summary")
rep_records_list <- lapply(results_list, `[[`, "rep_records")

# Persist: keep previous contract for summaries (list of tibbles) AND write the reps
saveRDS(summaries_list, file = file.path(results_dir, "sim_raw_results.rds"))

rep_records <- dplyr::bind_rows(rep_records_list)
saveRDS(rep_records, file = file.path(results_dir, "sim_rep_records.rds"))

# log
writeLines(c(
  paste("n_conditions:", length(summaries_list)),
  paste("n_reps:", n_reps),
  paste("n_rep_rows:", nrow(rep_records)),
  paste("simtime_total_seconds:", round(simtime_total))
), con = file.path(results_dir, "sim_log.txt"))

message("Saved simulation results to: ",
        file.path(results_dir, "sim_raw_results.rds"), " (summaries, list) and ",
        file.path(results_dir, "sim_rep_records.rds"), " (per-rep records).")

