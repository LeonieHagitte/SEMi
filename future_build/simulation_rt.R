library(tidyverse)
library(furrr)
library(dplyr)
library(tibble)

source(dataprep_rt.R)


manifest_path <- "progress_manifest.csv"
results_path  <- "results.csv"
lock_dir <- "locks"


MOD_TYPES <- c("linear","sigmoid","quadratic","noise")

DESIGN <- tidyr::expand_grid(
  popmodel     = c("NULL","1.1", "1.11", "1.12","1.2","1.21","1.22","1.3","1.32"),
  N            = c(300, 500, 700, 1000),
  reliability  = c(0.60, 0.70, 0.80, 0.95),
  lambda       = 0.70,
  intercepts   = c(-1, 0, 1),
  # latentmean  = 0,
  # delta_eta   = c(-1, -0.5, 0.5, 1),
  delta_lambda = c(-0.3, -0.2, 0.2, 0.3),
  delta_nu     = c(-1, -0.5, 0.5, 1),
  moderator    = MOD_TYPES,
  rep          = 1:10
)

# Optional: enforces a fixed ordering before row-number ids
DESIGN <- dplyr::arrange(
  DESIGN,
  popmodel, N, reliability, lambda, intercepts, delta_lambda, delta_nu, moderator, rep
)

DESIGN$job_id <- seq_len(nrow(DESIGN))
DESIGN$seed   <- DESIGN$job_id  #seed is set equal to job_id -simple deterministic mapping

manifest <- transform(
  DESIGN,
  status = "PENDING",           # PENDING | RUNNING | DONE | ERROR
  started_at = NA_character_,
  finished_at = NA_character_,
  worker = NA_character_,
  error_msg = NA_character_
)

if (!file.exists(manifest_path)) {
  write.csv(manifest, manifest_path, row.names = FALSE)
}
if (!dir.exists(lock_dir)) dir.create(lock_dir, recursive = TRUE)


if (!file.exists(results_path)) {
  write.csv(
    data.frame(
      job_id = integer(),
      popmodel = character(),
      N = integer(),
      reliability = numeric(),
      lambda = numeric(),
      intercepts = integer(),
      delta_lambda = numeric(),
      delta_nu = numeric(),
      moderator = character(),
      rep = integer(),
      seed = integer(),
      # add outcome columns here:
      result = numeric() #placeholder
    ),
    results_path,
    row.names = FALSE
  )
}

#######
run_one <- function(row) { #run_one <- function(seed, N, popmodel, moderator) 

  set.seed(row$seed)
  
  # ---- SIMULATION GOES HERE ----
  # Example placeholder:
  value <- runif(1)
  

    # generate
   # params <- gen_matrixA(popmodel=popmodel, moderator=moderator, m1=TRUE, m2=0)
  #  sim <- gen_dataA(N=N, params=params)
  #  df <- sim$data
    
    # ensure columns required by analysis exist
   # if (!"m2" %in% names(df)) df$m2 <- rep(params$m2, N)
    
    # analyze
  #  res <- mi_analysis(df)
    
    # attach design info
   # cbind(
    #  data.frame(seed=seed, N=N, popmodel=popmodel, moderator=moderator),
    #  res
    #)
  }
  
  
  # ----------------------------------
  
  tibble(
    job_id    = row$job_id,
    popmodel    = row$popmodel,
    N           = row$N,
    reliability = row$reliability,
    intercepts  = row$intercepts,
    delta_lambda= row$delta_lambda,
    delta_nu    = row$delta_nu,
    moderator   = row$moderator,
    rep         = row$rep,
    
    mnlfa_model = row$mnlfa_model,
    mnlfa_det   = row$mnlfa_det,
    assumed_mod_form = row$assumed_mod_form,
    tree_splits_m1 = row$tree_splits_m1,
    tree_splits_m2 = row$tree_splits_m2,
    
    delta_bias        = delta_bias,
    delta_rmse        = delta_rmse,
    delta_coverage    = delta_coverage,
    value       = value#placeholder
  )
}

t1 <- Sys.time()
out <- run_one(DESIGN[1, ])
t2 <- Sys.time()

elapsed_one <- as.numeric(difftime(t2, t1, units = "secs"))
elapsed_one
out