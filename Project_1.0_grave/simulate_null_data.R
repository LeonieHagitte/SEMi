#
# This script simulates "null" data that follows the form
# of the required data but has no meaningful predictors
#
casp_g_simulation <- 'control  =~ cC1 + cC2 + cC3 
                       autonomy  =~ cA1 + cA2 + cA3
                       pleasure  =~ cP1 + cP2 + cP3
                       self_real  =~ cS1 + cS2 + cS3
                     '
df4 <- lavaan::simulateData(casp_g_simulation,sample.nobs = 1000)


df4$age <- sample(18:100, nrow(df4),replace=TRUE)
df4$gender <- sample(c(0,1),nrow(df4), replace=TRUE)
