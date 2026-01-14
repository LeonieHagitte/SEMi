
# This script is for the simulation itself, the evaluation will be done in a second script.

# Setup ------------------------------------------------------------------------

# clear workspace to avoid mistakes
#rm(list=ls())

# load packages
library(here) # to source functions
library(tidyverse) # for handling data
library(lavaan) # for mgsem
library(mxsem) # for mnlfa
library(glue)

#install.packages("devtools")
#library(devtools)
library(doParallel) # for parallelisation
library(doRNG) # for parallelisation and random seeds

# source functions
source(here("Simulation_1" , "Functions_1", "sim_functions_1.R"))

# define folder to save results into
results_folder <- here("Results")

# save session info
writeLines(capture.output(sessionInfo()), here(results_folder, "sim_sessioninfo.txt"))

# Define a function to generate different types of moderator
generate_moderator <- function(N, type = c("linear", "quadratic", "sigmoid")) {
  type <- match.arg(type)
  base <- runif(N, -2, 2)
  
  moderator <- switch(type,
                      linear = base,
                      quadratic = base^2,
                      sigmoid = 1 / (1 + exp(-2 * (base - 0)))  # centered sigmoid
  )
  
  return(moderator)
}


# Data-generating Models -----------------------------------------------------------------------

function(typ des models, loading_offset)
  
  for i in persons:
  
  if (typ des models = 1)
    loading1 <- 0.7
loading2 <- 0.7
else if (typ des models == 2)
  loading <- loading_offset + moderator_value * effect_size
loading <- 0.7
else if (typ des models == 3)
  loading <- loading_offset + moderator_value * effect_size

model1_nonmod <- '
  # first item fixed to 1
  latent1 =~ ',loading1,'*item1 + ',loading2,'*item2 + 0.7*item3 + 1*item4

  # variances
  item1 ~~ item1
  item2 ~~ item2
  item3 ~~ item3
  item4 ~~ item4
  
  # latent variances 
  latent1 ~~ latent1


  # intercepts
  item1 ~ 0*1
  item2 ~ 0*1
  item3 ~ 0*1
  item4 ~ 0*1
  
  # latent means
  latent1 ~ 0*1
  '
lavaan::simulateData(model1_named, sample.nobs = 1)


build_moderated_metricmodel <- function(num_moderated) {
  stopifnot(num_moderated %in% 0:2)
  
  all_items <- 1:4
  candidate_items <- 1:2
  moderated_items <- candidate_items[1:num_moderated]
  
  # Loadings
  loadings <- purrr::map_chr(all_items, function(i) {
    if (i %in% moderated_items) {
      glue::glue("{{lambda_item{i} := lambda_item{i}_0 + lambda_item{i}_1*data.moderator}}*item{i}")
    } else if (i == 4) {
      "1*item4"
    } else {
      glue::glue("lambda_item{i}*item{i}")
    }
  })
  
  # Construct model string
  model_string <- glue::glue('
    latent1 =~ {paste(loadings, collapse = " + ")}

    # variances
    item1 ~~ item1
    item2 ~~ item2
    item3 ~~ item3
    item4 ~~ item4

    # latent variances 
    latent1 ~~ latent1

    # intercepts
    item1 ~ 0*1
    item2 ~ 0*1
    item3 ~ 0*1
    item4 ~ 0*1

    # latent means
    latent1 ~ 0*1
  ')
  
  return(model_string)
}


model_0 <- build_moderated_metricmodel(0)  # no moderation
model_1 <- build_moderated_metricmodel(1)  # one moderated loading
model_2 <- build_moderated_metricmodel(2)  # two moderated loadings


build_moderated_scalarmodel <- function(num_moderated_items) {
  stopifnot(num_moderated_items %in% 0:2)
  
  all_items <- 1:4
  candidate_items <- 1:2  # Only item1 and item2 can be moderated
  moderated_items <- candidate_items[1:num_moderated_items]  # Moderated subset
  
  # Loadings
  loadings <- purrr::map_chr(all_items, function(i) {
    if (i %in% moderated_items) {
      glue::glue("{{lambda_item{i} := lambda_item{i}_0 + lambda_item{i}_1*data.moderator}}*item{i}")
    } else if (i == 4) {
      "1*item4"  # fixed for identification
    } else {
      glue::glue("lambda_item{i}*item{i}")
    }
  })
  
  # Intercepts
  intercepts <- purrr::map_chr(all_items, function(i) {
    if (i %in% moderated_items) {
      glue::glue("item{i} ~ {{nu_item{i} := nu_item{i}_0 + nu_item{i}_1*data.moderator}}*1")
    } else {
      glue::glue("item{i} ~ nu_item{i}*1")
    }
  })
  
  # Construct full model string
  model_string <- glue::glue('
    latent1 =~ {paste(loadings, collapse = " + ")}

    # variances
    item1 ~~ item1
    item2 ~~ item2
    item3 ~~ item3
    item4 ~~ item4

    # latent variances 
    latent1 ~~ latent1

    # intercepts
    {paste(intercepts, collapse = "\\n    ")}

    # latent means
    latent1 ~ 0*1
  ')
  
  return(model_string)
}

model_0 <- build_moderated_scalarmodel(0)  # No moderation
model_1 <- build_moderated_scalarmodel(1)  # item1 moderated
model_2 <- build_moderated_scalarmodel(2)  # item1 and item2 moderated


build_moderated_strictmodel <- function(num_moderated_items) {
  stopifnot(num_moderated_items %in% 0:2)
  
  all_items <- 1:4
  candidate_items <- 1:2  # Only item1 and item2 can be moderated
  moderated_items <- candidate_items[1:num_moderated_items]
  
  # Loadings
  loadings <- purrr::map_chr(all_items, function(i) {
    if (i %in% moderated_items) {
      glue::glue("{{lambda_item{i} := lambda_item{i}_0 + lambda_item{i}_1*data.moderator}}*item{i}")
    } else if (i == 4) {
      "1*item4"  # fixed for identification
    } else {
      glue::glue("lambda_item{i}*item{i}")
    }
  })
  
  # Intercepts
  intercepts <- purrr::map_chr(all_items, function(i) {
    if (i %in% moderated_items) {
      glue::glue("item{i} ~ {{nu_item{i} := nu_item{i}_0 + nu_item{i}_1*data.moderator}}*1")
    } else {
      glue::glue("item{i} ~ nu_item{i}*1")
    }
  })
  
  # Residual variances
  residuals <- purrr::map_chr(all_items, function(i) {
    if (i %in% moderated_items) {
      glue::glue("item{i} ~~ {{resid_item{i} := resid_item{i}_0 + resid_item{i}_1*data.moderator}}*item{i}")
    } else {
      glue::glue("item{i} ~~ item{i}")
    }
  })
  
  # Construct full model string
  model_string <- glue::glue('
    latent1 =~ {paste(loadings, collapse = " + ")}

    # residual variances
    {paste(residuals, collapse = "\\n    ")}

    # latent variances 
    latent1 ~~ latent1

    # intercepts
    {paste(intercepts, collapse = "\\n    ")}

    # latent means
    latent1 ~ 0*1
  ')
  
  return(model_string)
}

# Example usage:
model_strict_0 <- build_moderated_strictmodel(0)  # No moderation
model_strict_1 <- build_moderated_strictmodel(1)  # item1 moderated
model_strict_2 <- build_moderated_strictmodel(2)  # item1 and item2 moderated


# Conditions -------------------------------------------------------------------

# define conditions
conditions <- expand.grid(n = c(300, 500, 700, 1000), 
                          moderator_effect = c(0, 0.075, 0.15))

moderator_type <- sample(c("linear", "quadratic", "sigmoid"), 1)

# Generate moderator vector
moderator <- generate_moderator(N = N, type = moderator_type)

# Create data.frame to be passed to lavaan or other estimation function
data <- data.frame(moderator = moderator)
# Preparation ------------------------------------------------------------------

# set seed for this computer
set.seed(20251)

# define number of runs
n_runs <- 100

# get number of logical processors minus two so I can still use my laptop and it's an even number
n_cores <- parallel::detectCores() - 2 # 12 - 2 = 10

## Start simulation -----

# define how many runs each core/thread should do
n_runs_per_core <- floor(n_runs/n_cores) # rounded down to closest integer

# create cluster
cluster <- makeCluster(n_cores)

# make it parallel
registerDoParallel(cluster)


# Simulation -------------------------------------------------------------------

# start time simulation
start_time_simulation <- Sys.time()

# export all objects that the workers need
clusterExport(cluster, varlist = c("models1", "conditions", "n_runs_per_core", "results_folder"))

# run simulation
results_parallel <- foreach(core_nr = 1:n_cores,
                            .packages = c("lavaan", # necessary packages
                                          "mxsem",
                                          ""),
                            .export = c("generate_data", # individual functions
                                        "estimate_mnlfa",
                                        "estimate_semtrees",
                                        "get_results")) %dorng% { # use doRNG for reproducibility
                                          get_results(conditions, 
                                                      n_runs = n_runs_per_core,
                                                      models1,#############################
                                                      results_folder)
                                        }

# end cluster
stopCluster(cluster)

# time for simulation
# end time simulation
end_time_simulation <- Sys.time()
# total time for simulation
simtime_total <- as.numeric(difftime(end_time_simulation, start_time_simulation, units = "secs"))
# define path for simulation log
simulation_log <- file.path(results_folder, "sim_simlog.txt")
# get time broken down into days, hours, minutes, seconds
days <- simtime_total %/% (24 * 3600)
hours <- (simtime_total %% (24 * 3600)) %/% 3600
minutes <- (simtime_total %% 3600) %/% 60
seconds <- round(simtime_total %% 60)
# print into logfile
cat(paste0(
  "\nTotal time for simulation: ",
  days, " days, ",
  hours, " hours, ",
  minutes, " minutes, ",
  seconds, " seconds\n"),
  file = simulation_log, append = TRUE)


# Save results -----------------------------------------------------------------

# save results as R-object
saveRDS(results_parallel, here(results_folder, "sim_results.rds"))


