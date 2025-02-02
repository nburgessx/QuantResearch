#####################################################################################
# Jim Gatheral, January 2024
#
# Implementations of the Black formula, implied volatility computations, greeks,
# and the Romano-Touzi formula.
#
#####################################################################################

BlackCall <- function (S0, K, T, sigma) 
{
  k <- log(K/S0)
  sig <- sigma * sqrt(T)
  d1 <- -k/sig + sig/2
  d2 <- d1 - sig
  return(S0 * pnorm(d1) - K * pnorm(d2))
}

BlackPut <- function (S0, K, T, sigma) 
{
  k <- log(K/S0)
  sig <- sigma * sqrt(T)
  d1 <- -k/sig + sig/2
  d2 <- d1 - sig
  return( K * pnorm(-d2) - S0 * pnorm(-d1) )
}

BlackOTM <- function (S0, K, T, sigma) 
{
  k <- log(K/S0)
  eps <- ifelse(k>=0,+1,-1)
  sig <- sigma * sqrt(T)
  d1 <- -k/sig + sig/2
  d2 <- d1 - sig
  return(eps*S0*pnorm(eps*d1) - eps*K*pnorm(eps*d2))
}

# This function works with vectors of strikes and option values
ivCall <- function(S0, K, T, C)
{
nK <- length(K)
sigmaL <- rep(1e-10,nK)
CL <- BlackCall(S0, K, T, sigmaL)
sigmaH <- rep(10,nK)
CH <- BlackCall(S0, K, T, sigmaH)
  while (mean(sigmaH - sigmaL) > 1e-10)
  {
    sigma <- (sigmaL + sigmaH)/2
    CM <- BlackCall(S0, K, T, sigma)
    CL <- CL + (CM < C)*(CM-CL)
    sigmaL <- sigmaL + (CM < C)*(sigma-sigmaL)
    CH <- CH + (CM >= C)*(CM-CH)
    sigmaH <- sigmaH + (CM >= C)*(sigma-sigmaH)
  }
  return(sigma)
}

# This function also works with vectors of strikes and option values  
ivPut <- function (S0, K, T, P) 
{
  nK <- length(K)
  sigmaL <- 1e-10
  PL <- BlackPut(S0, K, T, sigmaL)
  sigmaH <- 10
  PH <- BlackPut(S0, K, T, sigmaH)
  while (mean(sigmaH - sigmaL) > 1e-10) {
    sigma <- (sigmaL + sigmaH)/2
    PM <- BlackPut(S0, K, T, sigma)
    PL <- PL + (PM < P) * (PM - PL)
    sigmaL <- sigmaL + (PM < P) * (sigma - sigmaL)
    PH <- PH + (PM >= P) * (PM - PH)
    sigmaH <- sigmaH + (PM >= P) * (sigma - sigmaH)
  }
  return(sigma)
}

ivOTM <- function(S0, K, T, V.OTM){
  nK <- length(K)
  sigmaL <- 1e-10
  PL <- BlackOTM(S0, K, T, sigmaL)
  sigmaH <- 10
  PH <- BlackOTM(S0, K, T, sigmaH)
  while (mean(sigmaH - sigmaL) > 1e-10) {
    sigma <- (sigmaL + sigmaH)/2
    PM <- BlackOTM(S0, K, T, sigma)
    PL <- PL + (PM < V.OTM) * (PM - PL)
    sigmaL <- sigmaL + (PM < V.OTM) * (sigma - sigmaL)
    PH <- PH + (PM >= V.OTM) * (PM - PH)
    sigmaH <- sigmaH + (PM >= V.OTM) * (sigma - sigmaH)
  }
  return(sigma)
}

# Function to compute option prices and implied vols given vector of final values of underlying
ivS <- function (Sf, T, AK) 
{
  nK <- length(AK)
  N <- length(Sf)
  Sfbar <- mean(Sf)
  V <- numeric(nK)
  stErr <- numeric(nK)
  V.OTM <- numeric(nK)
  ivBlack <- numeric(nK)
  ivErr <- numeric(nK)
  for (j in 1:nK) {
    payoff <- (Sf - AK[j]) * (Sf > AK[j])
    
    V <- mean(payoff)
    stErr[j] <- sd(payoff)/sqrt(N)
    V.OTM[j] <- ifelse(Sfbar<AK[j], V, V + AK[j] - Sfbar)
    ivBlack[j] <- ivOTM(Sfbar, AK[j], T, V.OTM[j])
    ivErr[j] <-  ivOTM(Sfbar, AK[j], T, V.OTM[j]+stErr[j])-ivOTM(Sfbar, AK[j], T, V.OTM[j]-stErr[j])
  }
  return(ivBlack)
}

# Function to compute option prices and implied vols given vector of final values of underlying
ivS.list <- function (Sf, T, AK) 
{
  nK <- length(AK)
  N <- length(Sf)
  Sfbar <- mean(Sf)
  V <- numeric(nK)
  stErr <- numeric(nK)
  V.OTM <- numeric(nK)
  ivBlack <- numeric(nK)
  ivErr <- numeric(nK)
  for (j in 1:nK) {
    payoff <- (Sf - AK[j]) * (Sf > AK[j])
    
    V <- mean(payoff)
    stErr[j] <- sd(payoff)/sqrt(N)
    V.OTM[j] <- ifelse(Sfbar<AK[j], V, V + AK[j] - Sfbar)
    ivBlack[j] <- ivOTM(Sfbar, AK[j], T, V.OTM[j])
    ivErr[j] <-  ivOTM(Sfbar, AK[j], T, V.OTM[j]+stErr[j])-ivOTM(Sfbar, AK[j], T, V.OTM[j]-stErr[j])
  }
  res <- list(V.OTM=V.OTM, stErr=stErr, ivBlack=ivBlack,ivErr=ivErr)
  return(res)
}

# Functions to compute Black-Scholes greeks
delta.BS <- function (S0, K, T, sigma) 
{
  k <- log(K/S0)
  eps <- ifelse(k >= 0, +1, -1)
  sig <- sigma * sqrt(T)
  d1 <- -k/sig + sig/2
  d2 <- d1 - sig
  return(pnorm(eps * d1))
}

kappa.BS <- function (S0, K, T, sigma) 
{
  k <- log(K/S0)
  eps <- ifelse(k >= 0, +1, -1)
  sig <- sigma * sqrt(T)
  d1 <- -k/sig + sig/2
  d2 <- d1 - sig
  return(-pnorm(eps * d2))
}

vega.BS <- function (S0, K, T, sigma) 
{
  k <- log(K/S0)
  eps <- ifelse(k >= 0, +1, -1)
  sig <- sigma * sqrt(T)
  d1 <- -k/sig + sig/2
  d2 <- d1 - sig
  return(S0*dnorm(d1)*sqrt(T))
}

# Functions to compute option prices and sensitivities given a vector $Y_T = \int_t^T\,\sqrt{V_s}\,dW_s$
# and a vector $W_T =  \int_t^T\,V_s\,ds$ using the Romano-Touzi formula.

RomanoTouzi.OTM <- function(rho, y, w)Vectorize(function(k){
  
  s0 <- rho*y - rho^2/2*w
  eps <- ifelse(k>=0,+1,-1)
  sig <- sqrt((1-rho^2)*w)
  d1 <- (s0-k)/sig + sig/2
  d2 <- d1 - sig
  res <- eps*exp(s0)*pnorm(eps*d1) - eps*exp(k)*pnorm(eps*d2)
  return(mean(res))  
})

# Sensitivity wrt spot
RomanoTouzi.Delta <- function (rho, y, w) Vectorize(function(k)
{
  s0 <- rho * y - rho^2/2 * w
  eps <- ifelse(k >= 0, +1, -1)
  sig <- sqrt((1 - rho^2) * w)
  d1 <- (s0 - k)/sig + sig/2
  d2 <- d1 - sig
  res <- eps * exp(s0) * pnorm(eps * d1)
  return(mean(res))
})


# Sensitivity wrt strike
RomanoTouzi.Kappa <- function (rho, y, w)Vectorize(function(k) 
{
  s0 <- rho * y - rho^2/2 * w
  eps <- ifelse(k >= 0, +1, -1)
  sig <- sqrt((1 - rho^2) * w)
  d1 <- (s0 - k)/sig + sig/2
  d2 <- d1 - sig
  res <- -eps * pnorm(eps * d2)
  return(mean(res))
})

