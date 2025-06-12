install.packages("semTools")
library(semTools)
library(semtree)
library(lavaan)

# Remove rows with NA values
df5 <- na.omit(df5)

casp_modmxsem <- '
  # Latent variables
  control  =~ cC1 + cC2 + cC3 
  autonomy =~ cA1 + cA2 + cA3
  pleasure =~ cP1 + cP2 + cP3
  self_real =~ cS1 + cS2 + cS3

  # Regression paths from predictors to latent variables
  # control ~ b1*gender + b2*age
  # autonomy ~ b3*gender + b4*age
  # pleasure ~ b5*gender + b6*age
  # self_real ~ b7*gender + b8*age
'

fit <- lavaan::cfa(casp_modmxsem, data = df5)

tree1 <- semtree::semtree(fit, df5, # we need OpenMx for focus parameters!
                          #constraints = semtree.constraints(focus.parameters=fp),
                          control = semtree.control(method="score"))

plot(tree1)
