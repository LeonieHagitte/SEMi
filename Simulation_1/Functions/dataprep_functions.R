# ---- dataprep_functions.R ----
# Purpose: helpers for preparing condition grids and moderator vectors

suppressPackageStartupMessages({
  library(tidyverse)
})

# Generate a moderator vector of length N according to the specified type.
generate_moderator <- function(N, type = c("linear", "quadratic", "sigmoid")) {
  type <- match.arg(type)
  base <- runif(N, -2, 2)
  out <- switch(
    type,
    linear    = base,
    quadratic = base^2,
    sigmoid   = 1 / (1 + exp(-2 * base))
  )
  out
}

# Build a tidy conditions tibble.
# Each row defines one scenario; alt model strings can be injected downstream.
make_conditions_grid <- function(n_vec,
                                 moderator_types,
                                 mod_items_load_list,
                                 mod_items_int_list,
                                 mod_items_res_list,
                                 slope_load_vec,
                                 slope_int_vec,
                                 slope_logres_vec) {
  tibble::tibble(
    n = n_vec,
    moderator_type = moderator_types,
    mod_items_load = mod_items_load_list,
    mod_items_int  = mod_items_int_list,
    mod_items_res  = mod_items_res_list,
    slope_load     = slope_load_vec,
    slope_int      = slope_int_vec,
    slope_logres   = slope_logres_vec
  )
}
