install.packages("mxsem")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("jhorzek/mxsem", 
                         ref = "main")

library(mxsem)

casp_modmxsem <- '
# loadings
control  =~ lC1*cC1_  + lC2*cC2_  + lC3*cC3_  
autonomy =~ lA1*cA1_  + lA2*cA2_  + lA3*cA3_ 
pleasure =~ lP1*cP1_  + lP2*cP2_  + lP3*cP3_ 
self_real =~ lS1*cS1_  + lS2*cS2_  + lS3*cS3_ 

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
cC1_  ~~ vC1*cC1
cC2_   ~~ vC2*cC2
cC3_   ~~ vC3*cC3

cA1_   ~~ vA1*cA1
cA2_   ~~ vA2*cA2
cA3_   ~~ vA3*cA3

cP1_   ~~ vP1*cP1
cP2_   ~~ vP2*cP2
cP3_   ~~ vP3*cP3

cS1_   ~~ vS1*cS1
cS2_   ~~ vS2*cS2
cS3_   ~~ vS3*cS3

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
control  =~ {lC1_ := lC0_1 + lC1_1*data.age + lC2_1*data.gender + lC3_1*data.wave}*cC1_ +
            {lC2_ := lC0_2 + lC1_2*data.age + lC2_2*data.gender + lC3_2*data.wave}*cC2_ +
            {lC3_ := lC0_3 + lC1_3*data.age + lC2_3*data.gender + lC3_3*data.wave}*cC3_

autonomy =~ {lA1_ := lA0_1 + lA1_1*data.age + lA2_1*data.gender + lA3_1*data.wave}*cA1_ +
            {lA2_ := lA0_2 + lA1_2*data.age + lA2_2*data.gender + lA3_2*data.wave}*cA2_ +
            {lA3_ := lA0_3 + lA1_3*data.age + lA2_3*data.gender + lA3_3*data.wave}*cA3_

pleasure =~ {lP1_ := lP0_1 + lP1_1*data.age + lP2_1*data.gender + lP3_1*data.wave}*cP1_ +
            {lP2_ := lP0_2 + lP1_2*data.age + lP2_2*data.gender + lP3_2*data.wave}*cP2_ +
            {lP3_ := lP0_3 + lP1_3*data.age + lP2_3*data.gender + lP3_3*data.wave}*cP3_

self_real =~ {lS1_ := lS0_1 + lS1_1*data.age + lS2_1*data.gender + lS3_1*data.wave}*cS1_ +
             {lS2_ := lS0_2 + lS1_2*data.age + lS2_2*data.gender + lS3_1*data.wave}*cS2_ +
             {lS3_ := lS0_3 + lS1_3*data.age + lS2_3*data.gender + lS3_1*data.wave}*cS3_

# Explanation: Loadings specify the relationships between latent variables and observed indicators.
# Each loading is a function of intercepts, slopes (coefficients), and covariates (age and gender).

# Latent Variances and Covariances
control ~~ 1*control
autonomy ~~ 1*autonomy
pleasure ~~ 1*pleasure
self_real ~~ 1*self_real 
# Explanation: Latent variances are set to 1 for scaling purposes.
# Covariances between latent variables (control, autonomy, pleasure, self_real) are specified.

# Manifest Variances
cC1_  ~~ {vC1_ := exp(vC0_1 + vC1_1*data.age + vC2_1*data.gender + vC3_1*data.wave)}*cC1_
cC2_  ~~ {vC2_ := exp(vC0_2 + vC1_2*data.age + vC2_2*data.gender + vC3_2*data.wave)}*cC2_
cC3_  ~~ {vC3_ := exp(vC0_3 + vC1_3*data.age + vC2_3*data.gender + vC3_3*data.wave)}*cC3_

cA1_  ~~ {vA1_ := exp(vA0_1 + vA1_1*data.age + vA2_1*data.gender + vA3_1*data.wave)}*cA1_
cA2_  ~~ {vA2_ := exp(vA0_2 + vA1_2*data.age + vA2_2*data.gender + vA3_2*data.wave)}*cA2_
cA3_  ~~ {vA3_ := exp(vA0_3 + vA1_3*data.age + vA2_3*data.gender + vA3_3*data.wave)}*cA3_

cP1_  ~~ {vP1_ := exp(vP0_1 + vP1_1*data.age + vP2_1*data.gender + vP3_1*data.wave)}*cP1_
cP2_  ~~ {vP2_ := exp(vP0_2 + vP1_2*data.age + vP2_2*data.gender + vP3_2*data.wave)}*cP2_
cP3_  ~~ {vP3_ := exp(vP0_3 + vP1_3*data.age + vP2_3*data.gender + vP3_3*data.wave)}*cP3_

cS1_  ~~ {vS1_ := exp(vS0_1 + vS1_1*data.age + vS2_1*data.gender + vS3_1*data.wave)}*cS1_
cS2_  ~~ {vS2_ := exp(vS0_2 + vS1_2*data.age + vS2_2*data.gender + vS3_2*data.wave)}*cS2_
cS3_  ~~ {vS3_ := exp(vS0_3 + vS1_3*data.age + vS2_3*data.gender + vS3_3*data.wave)}*cS3_

# Explanation: Manifest variances are specified as exponential functions of intercepts, slopes, and covariates.

# Intercepts
cC1_  ~ {iC1_ := iC0_1 + iC1_1*data.age + iC2_1*data.gender + iC3_1*data.wave}*1
cC2_  ~ {iC2_ := iC0_2 + iC1_2*data.age + iC2_2*data.gender + iC3_2*data.wave}*1
cC3_  ~ {iC3_ := iC0_3 + iC1_3*data.age + iC2_3*data.gender + iC3_3*data.wave}*1

cA1_  ~ {iA1_ := iA0_1 + iA1_1*data.age + iA2_1*data.gender + iA3_1*data.wave}*1
cA2_  ~ {iA2_ := iA0_2 + iA1_2*data.age + iA2_2*data.gender + iA3_2*data.wave}*1
cA3_  ~ {iA3_ := iA0_3 + iA1_3*data.age + iA2_3*data.gender + iA3_3*data.wave}*1

cP1_  ~ {iP1_ := iP0_1 + iP1_1*data.age + iP2_1*data.gender + iP3_1*data.wave}*1
cP2_  ~ {iP2_ := iP0_2 + iP1_2*data.age + iP2_2*data.gender + iP3_2*data.wave}*1
cP3_  ~ {iP3_ := iP0_3 + iP1_3*data.age + iP2_3*data.gender + iP3_3*data.wave}*1

cS1_  ~ {iS1_ := iS0_1 + iS1_1*data.age + iS2_1*data.gender + iS3_1*data.wave}*1
cS2_  ~ {iS2_ := iS0_2 + iS1_2*data.age + iS2_2*data.gender + iS3_2*data.wave}*1
cS3_  ~ {iS3_ := iS0_3 + iS1_3*data.age + iS2_3*data.gender + iS3_3*data.wave}*1

# Explanation: Intercepts are functions of intercepts, slopes, and covariates, multiplied by 1 (constant term).
"

# + {covRA := cov00 + cov01*data.age + cov02*data.gender}*autonomy + {covRP := cov000 + cov001*data.age + cov002*data.gender}*pleasure


# Define your dataset 
# Make sure your dataset includes columns for all the variables used in the model.
library(readr)
df4 <- read_delim("df4_long.csv", delim = ";", 
                  escape_double = FALSE, col_types = cols(...1 = col_skip()), 
                  trim_ws = TRUE)
df4 <- as.data.frame(df4)


# Pass the syntax to the mxsem() function
library(mxsem)
casp_mnlfa_model <- mxsem(model = casp_mnlfamxsem,
                          data = df4,
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
lC1_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lC1",
                                                 progress_bar = FALSE)

lC2_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lC2",
                                                 progress_bar = FALSE)

lC3_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lC3",
                                                 progress_bar = FALSE)

lA1_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lA1",
                                                 progress_bar = FALSE)

lA2_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lA2",
                                                 progress_bar = FALSE)

lA3_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lA3",
                                                 progress_bar = FALSE)

lP1_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lP1",
                                                 progress_bar = FALSE)

lP2_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lP2",
                                                 progress_bar = FALSE)

lP3_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lP3",
                                                 progress_bar = FALSE)

lS1_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lS1",
                                                 progress_bar = FALSE)

lS2_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lS2",
                                                 progress_bar = FALSE)

lS3_individual <- get_individual_algebra_results(mxModel = casp_mnlfamxsem_result,
                                                 algebra_names = "lS3",
                                                 progress_bar = FALSE)

# Display the first few rows of individual parameters for manifest variables
head(lC1_individual$lC1)
head(lC2_individual$lC2)
head(lC3_individual$lC3)

head(lA1_individual$lA1)
head(lA2_individual$lA2)
head(lA3_individual$lA3)

head(lP1_individual$lP1)
head(lP2_individual$lP2)
head(lP3_individual$lP3)

head(lS1_individual$lS1)
head(lS2_individual$lS2)
head(lS3_individual$lS3)

## Plotting the Individual parameters
library(ggplot2)

# Plot individual parameters for lC1

ggplot(data = lC1_individual$lC1,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lC1") +
  geom_point()

# Plot individual parameters for lC2

ggplot(data = lC2_individual$lC2,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lC2") +
  geom_point()


# Plot individual parameters for lC3

ggplot(data = lC3_individual$lC3,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lC3") +
  geom_point()

# Plot individual parameters for lA1

ggplot(data = lA1_individual$lA1,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lA1") +
  geom_point()

# Plot individual parameters for lA2

ggplot(data = lA2_individual$lA2,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lA2") +
  geom_point()

# Plot individual parameters for lA3

ggplot(data = lA3_individual$lA3,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lA3") +
  geom_point()

# Plot individual parameters for lP1

ggplot(data = lP1_individual$lP1,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lP1") +
  geom_point()

# Plot individual parameters for lP2

ggplot(data = lP2_individual$lP2,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lP2") +
  geom_point()

# Plot individual parameters for lP3

ggplot(data = lP3_individual$lP3,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lP3") +
  geom_point()

# Repeat the above code for other manifest variables (lS1, lS2, lS3)

# Example for lS1

ggplot(data = lS1_individual$lS1,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lS1") +
  geom_point()

# Example for lS2

ggplot(data = cS2_individual$lS2,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lS2") +
  geom_point()

# Example for lS3

ggplot(data = cS3_individual$lS3,
       aes(x = age,
           y = algebra_result,
           color = factor(gender))) +
  ylab("Individual Parameter Value for lS3") +
  geom_point()

#
# AB->LH here is a proposal how we can derive the
# model predictions much faster than going through
# the algebra computation
#

# for {lC1 := lC0_1 + lC1_1*data.age + lC2_1*data.gender}*cC1
params <- omxGetParameters(casp_mnlfamxsem_result)
grid <- expand.grid(seq(min(df5$age),max(df5$age),1), df5$gender)
names(grid) <- c("age","gender")
prediction <- apply(grid, 1, function(x) {
  params["lS0_1"]+x[1]*params["lS1_1"]+x[2]*params["lS2_1"]
})
grid <- cbind(grid, prediction)
grid %>% ggplot(aes(x=age,y=prediction, group=gender,color=gender))+geom_line()
