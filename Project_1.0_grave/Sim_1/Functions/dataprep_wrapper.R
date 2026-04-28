library(lavaan)
library(Matrix)


#source("dataprep_functions.R")

gen_pop_model_data <- function(model_type,
                               N,
                               reliability,
                               ...  # <- pass any extra knobs straight to gen_mat()
) {
  # Build matrices for the requested model_type.
  # Examples of extra knobs to pass through '...':
  #   moderator_1_value = 0.5,
  #   moderator_1_type  = "linear",
  #   delta_lambda_full = 0.2,
  #   delta_theta_full  = 0.2,
  #   delta_lambda_12   = 0.15,
  #   psi_var           = 1.0
  matrices <- gen_mat(model_type, reliability = reliability, ...)
  
  # Creates population syntax with embedded numeric values (required for simulateData()).
  syntax <- gen_pop_model_syntax(matrices, include.values = TRUE)
  
  # Draws the dataset of size N from the population model.
  data <- simulateData(model = syntax, sample.nobs = N)
  
  list(data = data, syntax = syntax, matrices = matrices)
}
