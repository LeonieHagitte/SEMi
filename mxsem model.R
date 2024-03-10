install.packages("mxsem")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("jhorzek/mxsem", 
                         ref = "main")

library(mxsem)

casp_modmxsem <- '
# loadings
control  =~ lC1*cC1 + lC2*cC2 + lC3*cC3 
autonomy =~ lA1*cA1 + lA2*cA2 + lA3*cA3
pleasure =~ lP1*cP1 + lP2*cP2 + lP3*cP3
self_real =~ lS1*cS1 + lS2*cS2 + lS3*cS3

# Explanation: Loadings specify the relationships between latent variables and observed indicators.
# For example, control is predicted by cC1, cC2, and cC3.

# latent variances and covariances
control ~~ 1*control
autonomy ~~ 1*autonomy
pleasure ~~ 1*pleasure
self_real ~~ 1*self_real + cov*control + cov*autonomy + cov*pleasure

# Explanation: Latent variances are set to 1 for scaling purposes.
# Covariances between latent variables (control, autonomy, pleasure, self_real) are specified.

# manifest variances
cC1  ~~ vC1*cC1
cC2  ~~ vC2*cC2
cC3  ~~ vC3*cC3

cA1  ~~ vA1*cA1
cA2  ~~ vA2*cA2
cA3  ~~ vA3*cA3

cP1  ~~ vP1*cP1
cP2  ~~ vP2*cP2
cP3  ~~ vP3*cP3

cS1  ~~ vS1*cS1
cS2  ~~ vS2*cS2
cS3  ~~ vS3*cS3

# Explanation: Manifest variances represent the variability of each observed indicator.

# intercepts
cC1  ~ iC1*1
cC2  ~ iC2*1
cC3  ~ iC3*1

cA1  ~ iA1*1
cA2  ~ iA2*1
cA3  ~ iA3*1

cP1  ~ iP1*1
cP2  ~ iP2*1
cP3  ~ iP3*1

cS1  ~ iS1*1
cS2  ~ iS2*1
cS3  ~ iS3*1

# Explanation: Intercepts represent the expected values of the observed indicators when latent variables are zero.
# They capture the intercept or constant term in the regression of each observed indicator on its latent variable.
'
#######
library(mxsem)

# MNLFA structure adapted to casp_modmxsem
casp_mnlfamxsem <- "
# Loadings
control  =~ {lC1 := lC0_1 + lC1_1*data.age + lC2_1*data.gender}*cC1 +
            {lC2 := lC0_2 + lC1_2*data.age + lC2_2*data.gender}*cC2 +
            {lC3 := lC0_3 + lC1_3*data.age + lC2_3*data.gender}*cC3

autonomy =~ {lA1 := lA0_1 + lA1_1*data.age + lA2_1*data.gender}*cA1 +
            {lA2 := lA0_2 + lA1_2*data.age + lA2_2*data.gender}*cA2 +
            {lA3 := lA0_3 + lA1_3*data.age + lA2_3*data.gender}*cA3

pleasure =~ {lP1 := lP0_1 + lP1_1*data.age + lP2_1*data.gender}*cP1 +
            {lP2 := lP0_2 + lP1_2*data.age + lP2_2*data.gender}*cP2 +
            {lP3 := lP0_3 + lP1_3*data.age + lP2_3*data.gender}*cP3

self_real =~ {lS1 := lS0_1 + lS1_1*data.age + lS2_1*data.gender}*cS1 +
             {lS2 := lS0_2 + lS1_2*data.age + lS2_2*data.gender}*cS2 +
             {lS3 := lS0_3 + lS1_3*data.age + lS2_3*data.gender}*cS3

# Explanation: Loadings specify the relationships between latent variables and observed indicators.
# Each loading is a function of intercepts, slopes (coefficients), and covariates (age and gender).

# Latent Variances and Covariances
control ~~ 1*control
autonomy ~~ 1*autonomy
pleasure ~~ 1*pleasure
self_real ~~ 1*self_real + {cov := cov0 + cov1*data.age + cov2*data.gender}*control + {cov := cov01 + cov01*data.age + cov02*data.gender}*autonomy + {cov := cov001 + cov001*data.age + cov002*data.gender}*pleasure

# Explanation: Latent variances are set to 1 for scaling purposes.
# Covariances between latent variables (control, autonomy, pleasure, self_real) are specified.

# Manifest Variances
cC1  ~~ {vC1 := exp(vC0_1 + vC1_1*data.age + vC2_1*data.gender)}*cC1
cC2  ~~ {vC2 := exp(vC0_2 + vC1_2*data.age + vC2_2*data.gender)}*cC2
cC3  ~~ {vC3 := exp(vC0_3 + vC1_3*data.age + vC2_3*data.gender)}*cC3

cA1  ~~ {vA1 := exp(vA0_1 + vA1_1*data.age + vA2_1*data.gender)}*cA1
cA2  ~~ {vA2 := exp(vA0_2 + vA1_2*data.age + vA2_2*data.gender)}*cA2
cA3  ~~ {vA3 := exp(vA0_3 + vA1_3*data.age + vA2_3*data.gender)}*cA3

cP1  ~~ {vP1 := exp(vP0_1 + vP1_1*data.age + vP2_1*data.gender)}*cP1
cP2  ~~ {vP2 := exp(vP0_2 + vP1_2*data.age + vP2_2*data.gender)}*cP2
cP3  ~~ {vP3 := exp(vP0_3 + vP1_3*data.age + vP2_3*data.gender)}*cP3

cS1  ~~ {vS1 := exp(vS0_1 + vS1_1*data.age + vS2_1*data.gender)}*cS1
cS2  ~~ {vS2 := exp(vS0_2 + vS1_2*data.age + vS2_2*data.gender)}*cS2
cS3  ~~ {vS3 := exp(vS0_3 + vS1_3*data.age + vS2_3*data.gender)}*cS3

# Explanation: Manifest variances are specified as exponential functions of intercepts, slopes, and covariates.

# Intercepts
cC1  ~ {iC1 := iC0_1 + iC1_1*data.age + iC2_1*data.gender}*1
cC2  ~ {iC2 := iC0_2 + iC1_2*data.age + iC2_2*data.gender}*1
cC3  ~ {iC3 := iC0_3 + iC1_3*data.age + iC2_3*data.gender}*1

cA1  ~ {iA1 := iA0_1 + iA1_1*data.age + iA2_1*data.gender}*1
cA2  ~ {iA2 := iA0_2 + iA1_2*data.age + iA2_2*data.gender}*1
cA3  ~ {iA3 := iA0_3 + iA1_3*data.age + iA2_3*data.gender}*1

cP1  ~ {iP1 := iP0_1 + iP1_1*data.age + iP2_1*data.gender}*1
cP2  ~ {iP2 := iP0_2 + iP1_2*data.age + iP2_2*data.gender}*1
cP3  ~ {iP3 := iP0_3 + iP1_3*data.age + iP2_3*data.gender}*1

cS1  ~ {iS1 := iS0_1 + iS1_1*data.age + iS2_1*data.gender}*1
cS2  ~ {iS2 := iS0_2 + iS1_2*data.age + iS2_2*data.gender}*1
cS3  ~ {iS3 := iS0_3 + iS1_3*data.age + iS2_3*data.gender}*1

# Explanation: Intercepts are functions of intercepts, slopes, and covariates, multiplied by 1 (constant term).
"

# Define your dataset 
# Make sure your dataset includes columns for all the variables used in the model.

# Pass the syntax to the mxsem() function
library(mxsem)
casp_mnlfa_model <- mxsem(model = casp_mnlfamxsem,
                          data = df5,
                          scale_loadings = FALSE,
                          scale_latent_variances = FALSE)

# Fit the model using mxRun() (default optimization)
casp_mnlfamxsem_result <- mxTryHard(casp_mnlfa_model)

# Alternative: Fit the model using mxTryHard() for more robust optimization
# Note: mxTryHard() employs various strategies to find the optimal solution,
# but it might take longer compared to mxRun().
# Use mxTryHard() if mxRun() struggles to converge.

# Display the summary
summary(casp_mnlfamxsem_result)

## Plotting individual parameters
# Get individual parameters for manifest variables (cC1, cC2, ..., cS3)
cC1_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cC1",
                                                 progress_bar = FALSE)

cC2_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cC2",
                                                 progress_bar = FALSE)

cC3_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cC3",
                                                 progress_bar = FALSE)

cA1_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cA1",
                                                 progress_bar = FALSE)

cA2_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cA2",
                                                 progress_bar = FALSE)

cA3_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cA3",
                                                 progress_bar = FALSE)

cP1_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cP1",
                                                 progress_bar = FALSE)

cP2_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cP2",
                                                 progress_bar = FALSE)

cP3_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cP3",
                                                 progress_bar = FALSE)

cS1_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cS1",
                                                 progress_bar = FALSE)

cS2_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cS2",
                                                 progress_bar = FALSE)

cS3_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem,
                                                 algebra_names = "cS3",
                                                 progress_bar = FALSE)

# Display the first few rows of individual parameters for manifest variables
head(cC1_individual$cC1)
head(cC2_individual$cC2)
head(cC3_individual$cC3)

head(cA1_individual$cA1)
head(cA2_individual$cA2)
head(cA3_individual$cA3)

head(cP1_individual$cP1)
head(cP2_individual$cP2)
head(cP3_individual$cP3)

head(cS1_individual$cS1)
head(cS2_individual$cS2)
head(cS3_individual$cS3)

## Plotting the Individual parameters
library(ggplot2)

# Plot individual parameters for cC1

ggplot(data = cC1_individual$cC1,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cC1") +
  geom_point()

# Plot individual parameters for cC2

ggplot(data = cC2_individual$cC2,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cC2") +
  geom_point()


# Plot individual parameters for cC3

ggplot(data = cC3_individual$cC3,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cC3") +
  geom_point()

# Plot individual parameters for cA1

ggplot(data = cA1_individual$cA1,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cA1") +
  geom_point()

# Plot individual parameters for cA2

ggplot(data = cA2_individual$cA2,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cA2") +
  geom_point()

# Plot individual parameters for cA3

ggplot(data = cA3_individual$cA3,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cA3") +
  geom_point()

# Plot individual parameters for cP1

ggplot(data = cP1_individual$cP1,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cP1") +
  geom_point()

# Plot individual parameters for cP2

ggplot(data = cP2_individual$cP2,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cP2") +
  geom_point()

# Plot individual parameters for cP3

ggplot(data = cP3_individual$cP3,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cP3") +
  geom_point()

# Repeat the above code for other manifest variables (cS1, cS2, cS3)

# Example for cS1

ggplot(data = cS1_individual$cS1,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cS1") +
  geom_point()

# Example for cS2

ggplot(data = cS2_individual$cS2,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cS2") +
  geom_point()

# Example for cS3

ggplot(data = cS3_individual$cS3,
       aes(x = Age,
           y = algebra_result,
           color = factor(Male))) +
  ylab("Individual Parameter Value for cS3") +
  geom_point()

