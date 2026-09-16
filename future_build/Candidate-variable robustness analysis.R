#Candidate-variable robustness analysis
# SEM-tree variable selection with irrelevant candidate moderators
# ============================================================================
library(dplyr)
library(tibble)
library(future)
library(future.apply)

source("dataprep_rt.R")
source("analysis_rt.R")
#------------------------------------------------------------------------------

set.seed(42)

N_REP <- 1000L
N_NOISE <- 10L

ROBUSTNESS_DESIGN <- tibble(
  job_id = seq_len(N_REP),
  popmodel = "1.1",
  N = 500L,
  reliability = 0.75,
  lambda = 0.70,
  intercepts = 1,
  delta_lambda = 0.20,
  delta_nu = 0.50,
  moderator = "linear",
  method = "SEMTREE",
  rep_id = seq_len(N_REP),
  num_noisy_predictors = N_NOISE,
  seed = 500000L + seq_len(N_REP)
)
#-------------------------------------
test_result <- run_one(
  ROBUSTNESS_DESIGN[1, , drop = FALSE]
)

test_result %>%
  select(
    tree_metric_reject,
    tree_metric_split,
    tree_metric_split_on_am1,
    tree_metric_n_splits_am1,
    tree_metric_split_on_noisy,
    semtree_error_msg,
    error_msg
  )
#-------------------------------------------------------------------------------
n_workers <- max(
  1,
  parallelly::availableCores() - 1
)

future::plan(
  future::multisession,
  workers = n_workers
)

t1 <- Sys.time()

robustness_results_1000 <- future.apply::future_lapply(
  seq_len(nrow(ROBUSTNESS_DESIGN)),
  function(i) {
    run_one(
      ROBUSTNESS_DESIGN[i, , drop = FALSE]
    )
  },
  future.seed = TRUE
) %>%
  dplyr::bind_rows()

t2 <- Sys.time()

elapsed_total_min <- as.numeric(
  difftime(t2, t1, units = "mins")
)

elapsed_total_min

future::plan(future::sequential)

saveRDS(
  robustness_results_1000,
  "results_robustness_10_noise_1000.rds"
)
#-----------------------------------------------------------------------------
#robustness_results <- results %>%
#  dplyr::filter(num_noisy_predictors == 10)
#
#nrow(robustness_results)
#-------------------------
robustness_summary_1000 <- robustness_results_1000 %>%
  summarise(
    n = n(),
    
    metric_rejection_rate =
      mean(tree_metric_reject, na.rm = TRUE),
    
    # PRIMARY SUCCESS CRITERION:
    # true moderator selected at least once
    successful_moderator_identification_rate =
      mean(tree_metric_split_on_am1, na.rm = TRUE),
    
    # Secondary diagnostics
    any_noise_selected_rate =
      mean(tree_metric_split_on_noisy, na.rm = TRUE),
    
    true_and_noise_selected_rate =
      mean(
        tree_metric_split_on_am1 &
          tree_metric_split_on_noisy,
        na.rm = TRUE
      ),
    
    true_only_selected_rate =
      mean(
        tree_metric_split_on_am1 &
          !tree_metric_split_on_noisy,
        na.rm = TRUE
      ),
    
    noise_without_true_rate =
      mean(
        !tree_metric_split_on_am1 &
          tree_metric_split_on_noisy,
        na.rm = TRUE
      ),
    
    neither_selected_rate =
      mean(
        !tree_metric_split_on_am1 &
          !tree_metric_split_on_noisy,
        na.rm = TRUE
      ),
    
    n_errors =
      sum(!is.na(semtree_error_msg))
  )

robustness_summary_1000
# A tibble: 1 × 9
#n metric_rejection_rate successful_moderator_identification_rate any_noise_selected_rate true_and_noise_selected_rate true_only_selected_rate noise_without_true_rate neither_selected_rate n_errors
#<int>                 <dbl>                                    <dbl>                   <dbl>                        <dbl>                   <dbl>                   <dbl>                 <dbl>    <int>
#  1  1000                 0.791                                    0.774                   0.109                        0.092                   0.682                   0.017                 0.209        0
> 