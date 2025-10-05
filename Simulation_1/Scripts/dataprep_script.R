# ---- dataprep_script.R ----
# Purpose: construct the condition grid and store it for simulation

suppressPackageStartupMessages({           # base::suppressPackageStartupMessages — hides "package loaded" messages
  library(here)                            # here — convenient project-rooted file paths (Rproj root or .here)
  library(tidyverse)                       # tidyverse — collection (dplyr, purrr, tibble, etc.) for data manipulation
})

# source functions
source(here("Simulation_1", "Functions", "dataprep_functions.R"))
# base::source — executes R files; here::here builds the OS-independent path from project root.
# The sourced file is expected to define helper functions used to build the condition grid,
# e.g., `make_conditions_grid()` (a user-defined function in your project).

source(here("Simulation_1", "Functions", "sim_functions.R"))  # for model builders
# Loads your simulation helpers (user-defined), e.g., build_model_moderated_* functions
# that return model syntax strings needed below for `alt_model_string`.

# ensure results dir
results_dir <- here("Simulation_1", "Results")
# character scalar; absolute-like path anchored at the project root (here::here)

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
# base::dir.exists checks folder presence.
# base::dir.create makes the directory tree; recursive=TRUE creates parents if missing.

# define a minimal illustrative grid; edit as needed
conditions <- make_conditions_grid(
  n_vec = c(300, 500, 1000),                 # numeric vector — sample sizes to simulate (per condition)
  moderator_types = c("linear", "linear", "sigmoid"), # character vector — type of moderator to generate
  mod_items_load_list = list(integer(0), c(1), c(1,2)), # list-columns — which item loadings are moderated in DGP
  mod_items_int_list  = list(integer(0), integer(0), c(1)), # list-columns — which item intercepts are moderated
  mod_items_res_list  = list(integer(0), integer(0), integer(0)), # list-columns — which residual variances are moderated
  slope_load_vec   = c(0, 0.10, 0.15),       # numeric vector — slope size for loading moderation (on m)
  slope_int_vec    = c(0, 0, 0.10),          # numeric vector — slope size for intercept moderation (on m)
  slope_logres_vec = c(0, 0, 0)              # numeric vector — slope size for log-residual moderation (on m)
)
# make_conditions_grid() is a user-defined function (in dataprep_functions.R).
# It should return a tibble/data.frame where each row is a simulation condition.
# Typical columns it creates (per our prior setup):
#   - n (int) — sample size
#   - moderator_type (chr) — "linear"/"quadratic"/"sigmoid" (how to generate moderator values)
#   - mod_items_load / mod_items_int / mod_items_res (list<int>) — which items are moderated (indices 1..4)
#   - slope_load / slope_int / slope_logres (dbl) — effect sizes in the DGP
# Optionally it may also include helper labels like alt_model_label, etc.

# attach one alternative model string per row (you may replace these)
conditions <- conditions %>%
  mutate(                                            # dplyr::mutate — add/modify columns
    alt_model_string = list(                         # create a single list-column with 3 elements (recycled row-wise)
      build_model_moderated_loadings(integer(0)),   # sim_functions.R (user-defined); returns mxsem model string
      build_model_moderated_loadings(c(1)),         # Alt model: moderation on loading of item 1
      build_model_moderated_loadings(c(1,2))        # Alt model: moderation on loadings of items 1 and 2
    )
  )
# `alt_model_string` is a list-column of character model syntax; each row will receive
# the corresponding element (vector recycling in mutate + list() yields row-wise assignment
# because the list has same length as nrow(conditions)). Each string is used later by
# run_condition() / estimate_mxsem() as the alternative analysis model.

# save conditions
saveRDS(conditions, file = file.path(results_dir, "conditions.rds"))
# base::saveRDS — writes a single R object to disk in RDS format (binary).
# The file contains the tibble `conditions` with all columns described above.

message("Saved conditions to: ", file.path(results_dir, "conditions.rds"))
# base::message — prints an informational message (to stderr) with the output path.
