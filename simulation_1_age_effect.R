#
# This script simulates informative data
#
generate_data <- function(nobs=100, loading_matrix) {


casp_g_simulation <- paste0('control  =~ ',loading_matrix[1,1],'*cC1 + ',loading_matrix[2,1],'*cC2 + ',loading_matrix[3,1],'*cC3 
                       autonomy  =~ ',loading_matrix[1,2],'*cA1 + ',loading_matrix[2,2],'*cA2 + ',loading_matrix[3,2],'*cA3
                       pleasure  =~ ',loading_matrix[1,3],'*cP1 + ',loading_matrix[2,3],'*cP2 + ',loading_matrix[3,3],'*cP3
                       self_real  =~ ',loading_matrix[1,4],'*cS1 + ',loading_matrix[2,4],'*cS2 + ',loading_matrix[3,4],'*cS3
                     ')

# strange work-around; bug in lavaan code?
if (nobs==1) {
  dat = lavaan::simulateData(
    casp_g_simulation,sample.nobs = 2)
  dat <- dat[1,]
} else {
dat = lavaan::simulateData(
  casp_g_simulation,sample.nobs = nobs)
}

return(dat)
}

loading_matrix <- matrix(
  c(0.7,0.7,0.7,
    0.7,0.7,0.7,
    0.7,0.7,0.7,
    0.7,0.7,0.7),nrow=3,ncol=4)

df4 <- c()

for (i in 1:10000) {
  age <- sample(18:100, 1,replace=TRUE)
  age_std <- (age-18)/(100-18)
  gender <- sample(c(0,1), 1, replace=TRUE)
  # simulate an age effect on the loading 
  # structure of the control factor
  # the older, the weaker the loadings
  loading_matrix[1,1] <- 0.7-age_std*0.35
  loading_matrix[2,1] <- 0.7-age_std*0.35
  loading_matrix[3,1] <- 0.7-age_std*0.35
  new_row <- generate_data(n=1, loading_matrix)
  new_row$age = age
  new_row$gender = gender
  if (is.null(df4)) {
    df4 <- new_row
  } else {
    df4 <- rbind(df4, new_row)
  }
}

df4$gender <- factor(df4$gender)

#
# run a score-based SEM tree
#

casp_g <- 'control  =~ cC1 + cC2 + cC3 
                       autonomy  =~ cA1 + cA2 + cA3
                       pleasure  =~ cP1 + cP2 + cP3
                       self_real  =~ cS1 + cS2 + cS3
                     '

library(semtree)

fit <- lavaan::cfa(casp_g, df4)
tree <- semtree(model = fit, data=df4,
                # we need OpenMx for focus parameters!
                #constraints = semtree.constraints(focus.parameters=fp),
                control = semtree.control(method="score"))

plot(tree)

###############################################################################
# run a MNLFA
#
# Saving data frame as mx object:
mxdf1 <- mxData(observed = df4, type = "raw")
manVars <- colnames(df4[,-c(1,2)])
nv <- length(manVars)

#################################
# With Gender and Age
## Specify matrices for configural model
matT0 <- mxMatrix(type="Full", nrow=1, ncol=12,#baseline intercepts
                  free=TRUE,
                  values=1,
                  name="matT0")

matB1 <- mxMatrix(type = "Full", nrow = 1, ncol = nv,  # full matrix of background effects
                  free = TRUE,
                  values = 0,
                  name = "matB1")
matB2 <- mxMatrix(type = "Full", nrow = 1, ncol = nv,
                  free = TRUE,
                  values = 0,
                  name = "matB2")

matL0 <- mxMatrix(type = "Full", nrow = 12, ncol = 4,  # loadings ############
                  free = c(rep(c(TRUE, FALSE, FALSE, FALSE), 3),
                           rep(c(FALSE, TRUE, FALSE, FALSE), 3),
                           rep(c(FALSE, FALSE, TRUE, FALSE), 3),
                           rep(c(FALSE, FALSE, FALSE, TRUE), 3)
                  ),
                  values = c(rep(c(1, 0.7, 0.7, 0.7), 3),
                             rep(c(0.7, 1, 0.7, 0.7), 3),
                             rep(c(0.7, 0.7, 1, 0.7), 3),
                             rep(c(0.7, 0.7, 0.7, 1), 3)),
                  byrow = TRUE,
                  name = "matL0")
# direct effects of age and gender
matC1 <- mxMatrix(type = "Full", nrow = 12, ncol = 4,  
                  free = c(rep(c(TRUE, FALSE, FALSE, FALSE), 3),
                           rep(c(FALSE, TRUE, FALSE, FALSE), 3),
                           rep(c(FALSE, FALSE, TRUE, FALSE), 3),
                           rep(c(FALSE, FALSE, FALSE, TRUE), 3)),
                  values = c(rep(c(-0.35, 0, 0, 0), 3),
                             rep(c(0, 0, 0, 0), 3),
                             rep(c(0, 0, 0, 0), 3),
                             rep(c(0, 0, 0, 0), 3)),
                  byrow = TRUE,
                  name = "matC1")

matC2 <- mxMatrix(type = "Full", nrow = 12, ncol = 4,
                  free = c(rep(c(TRUE, FALSE, FALSE, FALSE), 3),
                           rep(c(FALSE, TRUE, FALSE, FALSE), 3),
                           rep(c(FALSE, FALSE, TRUE, FALSE), 3),
                           rep(c(FALSE, FALSE, FALSE, TRUE), 3)),
                  values = 0,
                  byrow = TRUE,
                  name = "matC2")
# matrix for residual covariance at baseline
matE0 <- mxMatrix(type = "Diag", nrow = nv, ncol = nv,  
                  free = TRUE,
                  values = 1,
                  name = "matE0")
# matrix for residual covariance of age
matD1 <- mxMatrix(type = "Diag", nrow = nv, ncol = nv,  
                  free = TRUE,
                  values = 0,
                  name = "matD1")
# matrix for residual covariance of gender
matD2 <- mxMatrix(type = "Diag", nrow = nv, ncol = nv,  
                  free = TRUE,
                  values = 0,
                  name = "matD2")

matP0 <- mxMatrix(type = "Symm", nrow = 4, ncol = 4,#variances
                  free = c(FALSE, TRUE, TRUE, TRUE,
                           TRUE, FALSE, TRUE, TRUE,
                           TRUE, TRUE, FALSE, TRUE,
                           TRUE, TRUE, TRUE, FALSE),
                  values = c(1, 0, 0, 0,
                             0, 1, 0, 0,
                             0, 0, 1, 0,
                             0, 0, 0, 1),
                  name = "matP0")
matH1 <- mxMatrix(type = "Symm", nrow = 4, ncol = 4,
                  free = c(FALSE, TRUE, TRUE, TRUE,
                           TRUE, FALSE, TRUE, TRUE,
                           TRUE, TRUE, FALSE, TRUE,
                           TRUE, TRUE, TRUE, FALSE),
                  values = 0,
                  name = "matH1")
matH2 <- mxMatrix(type = "Symm", nrow = 4, ncol = 4,
                  free = c(FALSE, TRUE, TRUE, TRUE,
                           TRUE, FALSE, TRUE, TRUE,
                           TRUE, TRUE, FALSE, TRUE,
                           TRUE, TRUE, TRUE, FALSE),
                  values = 0,
                  name = "matH2")

matA0 <- mxMatrix(type="Full", nrow=4, ncol=1,#factor means not estimated
                  free=FALSE,
                  values= c(0, 0, 0, 0),
                  name="matA0")
matG1 <- mxMatrix(type="Full", nrow=4, ncol=1,
                  free=FALSE, # to identify the model config to zero
                  values=0,
                  name="matG1")
matG2 <- mxMatrix(type="Full", nrow=4, ncol=1,
                  free=FALSE, # to identify the model config to zero
                  values=0,
                  name="matG2")

matV1 <- mxMatrix(type="Full", nrow=1, ncol=1, 
                  free=FALSE, 
                  labels="data.age", 
                  name = "age")
matV2 <- mxMatrix(type="Full", nrow=1, ncol=1, 
                  free=FALSE, 
                  labels="data.gender", 
                  name = "gender")

matIa <- mxMatrix(type="Diag", nrow=4, ncol=4, 
                  free=FALSE,
                  values=1, 
                  name="matIa")
matIb <- mxMatrix(type="Full", nrow=4, ncol=4, 
                  free=FALSE, 
                  values=c(0,1,1,1,
                           1,0,1,1,
                           1,1,0,1,
                           1,1,1,0),
                  name="matIb")


## Specify algebra for the dependent parameters

matT <- mxAlgebra(expression=matT0+matB1*gender+matB2*age, 
                  name="matT")
matL <- mxAlgebra(expression=matL0+matC1*gender+matC2*age, 
                  name="matL")
matE <- mxAlgebra(expression=matE0*exp(matD1*gender+matD2*age), 
                  name="matE")
matA <- mxAlgebra(expression=matA0+matG1*gender+matG2*age, 
                  name="matA")

## Specify algebra for covariance matrix of factors (transformed to ensure positive definite matrices)
matVar <- mxAlgebra(expression=(matP0*exp(matH1*gender+matH2*age)), 
                    name="matVar")
matR <- mxAlgebra(expression=(exp(2*(matP0+matH1*gender+matH2*age))-1)/
                    (exp(2*(matP0+matH1*gender+matH2*age))+1), 
                  name="matR")
matCov <- mxAlgebra(expression=(matIa*sqrt(matVar))%*%matR%*%(matIa*sqrt(matVar)), 
                    name="matCov")
matP <- mxAlgebra(expression=matIa*matVar+matIb*matCov, 
                  name="matP")

## Specify model-implied matrices
matC <- mxAlgebra(expression=matL%*%matP%*%t(matL)+matE, 
                  name="matC") 
matM <- mxAlgebra(expression=matT+t(matL%*%matA), 
                  name="matM") 

## Specify expectation and fit function
expF <- mxExpectationNormal(covariance="matC", 
                            means="matM",
                            dimnames=manVars)
fitF <- mxFitFunctionML() 

## Make mxModel object and run the model

modConfig <- mxModel(model="Configural", 
                     matT, matT0, matB1, matB2,
                     matL, matL0, matC1, matC2, 
                     matE, matE0, matD1, matD2,
                     matP, matP0, matH1, matH2,
                     matA, matA0, matG1, matG2,  
                     matIa, matIb, matV1, matV2, 
                     matVar, matR, matCov, matM, matC, 
                     expF, fitF, mxdf1)

fitConfig <- mxRun(modConfig)
summary(fitConfig) 
