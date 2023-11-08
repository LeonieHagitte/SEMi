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
