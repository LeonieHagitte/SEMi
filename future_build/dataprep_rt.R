# dataprep_rt.R
# Two functions, one population parameter generator and one data generator
# mxdata1 is used data label in analysis_rt.R

mod_h <- function(M, type = c("linear","sigmoid","quadratic","noise"), k = 2.0){

   if (any(!is.finite(M))) stop("M must be finite.")
   if (any(M <= -1 | M >= 1)) stop("M must lie in the interval (-1, 1).")

   type <- match.arg(type)
   switch(type,
         linear    = {M},
         sigmoid   =  {if (!is.numeric(k) || length(k) != 1 || !is.finite(k) || k <= 0) {
           stop("k must be a single finite number > 0.")
         }
         -1 + 2 / (1 + exp(-k * M))},
         quadratic = { 2 * M^2 - 1},
         noise     = { rep(0, length(M))}
   )
}

gen_matrixA <- function(popmodel,
                       nfactors = 1,
                       p = 4,
                       lambda = 0.70,
                       nu = 1,
                       reliability = 0.8,
                       moderator = c("linear","sigmoid","quadratic","noise"),
                       m1 = TRUE,
                       m2 = TRUE,
                       k = 2,
                       delta_lambda = 0.2,
                       delta_nu = 0.5,
                       delta_eta = 0.5,
                       psi = 1,
                       mu_eta = 0.0) {
  
  moderator <- match.arg(moderator)
  # ---------------- Moderators (both are always returned) ----------------
  # moderator_1 may be informative; moderator_2 is fixed noise (non-informative).
  eps <- .Machine$double.eps
  m1 <- if (isTRUE(m1)) runif(1, -1 + eps, 1 - eps) else as.numeric(m1)
  m2 <- if (isTRUE(m2)) runif(1, -1 + eps, 1 - eps) else as.numeric(m2)
  m0 <-0
 
  hm1 <- mod_h(m1, type = moderator, k = k)
  hm2 <- mod_h(m2, type = moderator, k = k)
  hm12 <- hm1* hm2
  hm0 <- 0  # by design, moderator_0 is always "noise" - ##### currently not returned #####
  
  # ---------------- Base (unmoderated) model ---------------- 
  # Lambda Matrix
  LAMBDA <- matrix(0, nrow = p, ncol = nfactors)
  # set loadings
  LAMBDA[, 1] <- lambda
  
  # Latent Variance
  PSI <- matrix(psi, 1, 1)
  
  # Intercepts
  NU <- matrix(0, nrow = p, ncol = 1)
  NU[, 1] <- nu
  
  # Latent Mean 
  MU_eta <- mu_eta #MU_eta <- mu_eta + delta_eta * hm1 - not sure if latent mean moderation should be included
  
  # ---------------- Moderation conditions ------------------
  if (popmodel == "1.1") #1.1 = Full loading moderation, no intercepts
  {
    LAMBDA[,1] <- LAMBDA[,1] + delta_lambda * hm1
  }
  else if (popmodel == "1.11") #1.11 = no moderation of loadings, full moderation of intercepts
  {
    NU[,1] <- NU[,1] + delta_nu * hm1
  }
  else if (popmodel == "1.12") #1.12 = full loadings, and intercepts
  {
    LAMBDA[,1] <- LAMBDA[,1] + delta_lambda * hm1
    NU[,1] <- NU[,1] + delta_nu * hm1
  }
  else if (popmodel == "1.2") #1.2 = partial loading moderation, no intercepts
  {
    LAMBDA[1:2, 1] <- LAMBDA[1:2, 1] + delta_lambda * hm1
  } 
  else if (popmodel == "1.21") #1.21 = no loading moderation, partial intercept moderation
  {
    NU[1:2, 1] <- NU[1:2, 1] + delta_nu * hm1
  }
  else if (popmodel == "1.22") #1.22 = partial loading moderation, and intercepts
  {
    LAMBDA[1:2, 1] <- LAMBDA[1:2, 1] + delta_lambda * hm1
    NU[1:2, 1] <- NU[1:2, 1] + delta_nu * hm1
  }
  else if (popmodel == "1.3") #1.3 = full loadings moderated by hm1 and moderated by hm2
  {
    LAMBDA[,1] <- LAMBDA[,1] + delta_lambda * hm1 + delta_lambda * hm2 + (delta_lambda*delta_lambda) * hm12 # interaction version 
  }
  
  # Residual Variance
  lambda_i <- LAMBDA[, 1]
  theta_vec <- (lambda_i^2 * psi * (1 - reliability)) / reliability
  
  THETA <- diag(theta_vec)
  
  MLIST <- list(LAMBDA = LAMBDA, PSI = PSI, THETA = THETA, NU = NU, MU_eta = MU_eta, hm1 = hm1, hm2 = hm2, hm12 = hm12, m1 = m1, m2 = m2, m0 = m0)
  return(MLIST)
}

# ------------- Generating Data ----------------------------

gen_dataA <- function(N, params, return_latent = TRUE){
  stopifnot(is.list(params), all(c("LAMBDA","NU","PSI","THETA","MU_eta","m1","hm1", "m2","hm2","hm12","m0") %in% names(params)))
  
  LAMBDA <- params$LAMBDA
  NU     <- params$NU
  PSI    <- params$PSI
  THETA  <- params$THETA
  MU_eta <- params$MU_eta
  
  p <- nrow(LAMBDA)
  psi <- PSI[1,1]
  
  eta <- rnorm(N, mean = MU_eta, sd = sqrt(psi))
  
  E <- matrix(rnorm(N*p), nrow = N, ncol = p) *
    matrix(sqrt(diag(THETA)), nrow = N, ncol = p, byrow = TRUE)
  
  NU_mat <- matrix(NU[,1], nrow = p, ncol = N)
  X <- t(NU_mat + LAMBDA %*% t(eta)) + E
  
  colnames(X) <- paste0("x", 1:p)
  dat <- as.data.frame(X)
  
  # moderator columns (trial-level replicated)
  dat$m1   <- rep(params$m1,  N)
  dat$hm1  <- rep(params$hm1, N)
  dat$m2   <- rep(params$m2,  N)
  dat$hm2  <- rep(params$hm2, N)
  dat$hm12 <- rep(params$hm12, N)
  dat$m0   <- rep(params$m0,  N)
  
  out <- list(data = dat, params = params)
  if (return_latent) out$eta <- eta
  out
}

# ------------------------ Test ---------------------------
set.seed(1)
pa <- gen_matrixA(popmodel="1.12", moderator="linear", m1=TRUE)
sa <- gen_dataA(N=200, params=pa)
head(sa$data)
data <- sa$data
head(data)

# ------------------------ Person-level moderation version -------------------

gen_paramsC <- function(popmodel,
                        lambda = 0.70,
                        nu = 1,
                        reliability = 0.8,
                        moderator = c("linear","sigmoid","quadratic","noise"),
                        k = 2,
                        delta_lambda = 0.2,
                        delta_nu = 0.5,
                        delta_eta = 0.5,
                        psi = 1,
                        mu_eta = 0.0) {
  
  moderator <- match.arg(moderator)
  
  list(popmodel = popmodel,
       lambda = lambda,
       nu = nu,
       reliability = reliability,
       moderator = moderator,
       k = k,
       delta_lambda = delta_lambda,
       delta_nu = delta_nu,
       delta_eta = delta_eta,
       psi = psi,
       mu_eta = mu_eta)
}


gen_dataC <- function(N, params, return_latent=TRUE){
  
  p <- 4
  
  lambda0 <- params$lambda
  nu0     <- params$nu
  rel     <- params$reliability
  moderator <- params$moderator
  k       <- params$k
  dlam    <- params$delta_lambda
  dnu     <- params$delta_nu
  deta    <- params$delta_eta ########### currently not used #########
  psi     <- params$psi
  mu0     <- params$mu_eta
  popmodel<- params$popmodel
  
  # person-level moderator
  eps <- .Machine$double.eps
  m1 <- runif(N, -1 + eps, 1 - eps)
  m2 <- runif(N, -1+eps, 1-eps)
  hm1  <- mod_h(m1, type = moderator, k = k)
  hm2 <- mod_h(m2, type=moderator, k=k)
  hm12 <- hm1 * hm2
    
  # item moderation pattern
  dlam1_vec  <- rep(0,p);
  dnu1_vec <- rep(0,p)
  dlam2_vec  <- rep(0,p);
  dnu2_vec <- rep(0,p)
  dlam12_vec <- rep(0,p);
  dnu12_vec <- rep(0,p)
  
  if (popmodel == "1.1") {
    dlam1_vec[] <- dlam
  } else if (popmodel == "1.11") {
    dnu1_vec[] <- dnu
  } else if (popmodel == "1.12") {
    dlam1_vec[] <- dlam; dnu1_vec[] <- dnu
  } else if (popmodel == "1.2") {
    dlam1_vec[1:2] <- dlam
  } else if (popmodel == "1.21") {
    dnu1_vec[1:2] <- dnu
  } else if (popmodel == "1.22") {
    dlam1_vec[1:2] <- dlam; dnu1_vec[1:2] <- dnu
  } else if (popmodel == "1.3") {
    dlam1_vec[] <- dlam; dlam2_vec[] <- dlam; dlam12_vec[] <- dlam^2
  } else if (popmodel == "1.32") {
    dlam1_vec[] <- dlam; dnu2_vec[] <- dnu
  }
  else stop("Unknown popmodel: ", popmodel)
  
  # person-specific loadings/intercepts
  Lambda_Np <- matrix(lambda0, N, p) + hm1  %o% dlam1_vec + hm2  %o% dlam2_vec + hm12 %o% dlam12_vec
  
  Nu_Np <- matrix(nu0, N, p) + hm1  %o% dnu1_vec + hm2  %o% dnu2_vec + hm12 %o% dnu12_vec
 
   # latent
  MU_eta_i <- rep(mu0, N) #+ deta*hm1
  eta <- rnorm(N, mean = MU_eta_i, sd = sqrt(psi))
  
  # residual variances from reliability (person-specific)
  ThetaVar_Np <- (Lambda_Np^2 * psi * (1-rel)) / rel
  
  E <- matrix(rnorm(N*p), N, p) * sqrt(ThetaVar_Np)
  
  X <- Nu_Np + Lambda_Np * eta + E
  
  colnames(X) <- paste0("x",1:p)
  data <- as.data.frame(X)
  data$m1 <- m1
  data$hm1  <- hm1
  data$m2 <- m2
  data$hm2 <- hm2
  data$hm12 <- hm12
  
  out <- list(data=data, params=params)
  if (return_latent) out$eta <- eta
  out
}
# ------------- Test ---------------------------
#set.seed(1)
#pc <- gen_paramsC(popmodel="1.12", moderator="linear")
#sc <- gen_dataC(N=20, params=pc)
#head(sc$data)
#data <- sc$data
#head(data)


