install.packages("mxsem")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("jhorzek/mxsem", 
                         ref = "main")

library(mxsem)
casp_modmxsem <- '
  # Latent variables
  control  =~ cC1 + cC2 + cC3 
  autonomy =~ cA1 + cA2 + cA3
  pleasure =~ cP1 + cP2 + cP3
  self_real =~ cS1 + cS2 + cS3

  # Regression paths from predictors to latent variables
  control ~ b1*gender + b2*age
  autonomy ~ b3*gender + b4*age
  pleasure ~ b5*gender + b6*age
  self_real ~ b7*gender + b8*age
'

mxsem(model = casp_modmxsem,
      data  = df5) |>
  mxTryHard() |>
  summary()
