# dataprep_rt.R
# Two functions, one population parameter generator and one data generator
# mxdata1 is used data label in analysis_rt.R

mod_h <- function(M, type = c("linear","sigmoid","quadratic","noise"), k = 10.0){ #what k to choose?

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

# ------------------------ Person-level moderation version -------------------

gen_paramsC <- function(popmodel,
                        lambda = 0.70,
                        nu = 1,
                        reliability = 0.8,
                        moderator = c("linear","sigmoid","quadratic","noise"),
                        k = 10,
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
  m0 <- runif(N, -1 + eps, 1 - eps)   # uninformed moderator
  
  hm1  <- mod_h(m1, type = moderator, k = k)
  hm2 <- mod_h(m2, type=moderator, k=k)
  hm12 <- hm1 * hm2
  hm0  <- rep(0, N)   # keep uninformed
    
  # item moderation pattern
  dlam1_vec  <- rep(0,p);
  dnu1_vec <- rep(0,p)
  dlam2_vec  <- rep(0,p);
  dnu2_vec <- rep(0,p)
  dlam12_vec <- rep(0,p);
  dnu12_vec <- rep(0,p)
  
  if (popmodel == "0"|| popmodel == "NULL"){
    #no moderation; keep all d-vectors at 0
  } else if (popmodel == "1.1") {
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
  #ThetaVar_Np <- (Lambda_Np^2 * psi * (1-rel)) / rel
  theta0 <- (lambda0^2 * psi * (1 - rel)) / rel
  ThetaVar_Np <- matrix(theta0, nrow = N, ncol = p, byrow = TRUE)
  E <- matrix(rnorm(N*p), N, p) * sqrt(ThetaVar_Np)
  
  X <- Nu_Np + Lambda_Np * eta + E
  
  colnames(X) <- paste0("x",1:p)
  data <- as.data.frame(X)
  data$m1 <- m1
  data$hm1  <- hm1
  data$m2 <- m2
  data$hm2 <- hm2
  data$hm12 <- hm12
  data$m0 <- m0
  data$hm0 <- hm0
  
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


