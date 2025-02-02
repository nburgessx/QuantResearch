BSFormula <- function(S0, K, T, r, sigma)
{
x <- log(S0/K)+r*T;
sig <- sigma*sqrt(T);
d1 <- x/sig+sig/2;
d2 <- d1 - sig;
pv <- exp(-r*T);
return( S0*pnorm(d1) - pv*K*pnorm(d2));
}

BSFormulaPut <- function (S0, K, T, r, sigma) 
{
  x <- log(S0/K) + r * T
  sig <- sigma * sqrt(T)
  d1 <- x/sig + sig/2
  d2 <- d1 - sig
  pv <- exp(-r * T)
  return( pv * K * pnorm(-d2) - S0 * pnorm(-d1))
}

# This function now works with vectors of strikes and option values
BSImpliedVolCall <- function(S0, K, T, r, C)
{
nK <- length(K);
sigmaL <- rep(1e-10,nK);
CL <- BSFormula(S0, K, T, r, sigmaL);
sigmaH <- rep(10,nK);
CH <- BSFormula(S0, K, T, r, sigmaH);
  while (mean(sigmaH - sigmaL) > 1e-10)
  {
    sigma <- (sigmaL + sigmaH)/2;
    CM <- BSFormula(S0, K, T, r, sigma);
    CL <- CL + (CM < C)*(CM-CL);
    sigmaL <- sigmaL + (CM < C)*(sigma-sigmaL);
    CH <- CH + (CM >= C)*(CM-CH);
    sigmaH <- sigmaH + (CM >= C)*(sigma-sigmaH);
  }
  return(sigma);
}

# This function also works with vectors of strikes and option values  
BSImpliedVolPut <- function (S0, K, T, r, P) 
{
  nK <- length(K)
  sigmaL <- 1e-10
  PL <- BSFormulaPut(S0, K, T, r, sigmaL)
  sigmaH <- 10
  PH <- BSFormulaPut(S0, K, T, r, sigmaH)
  while (mean(sigmaH - sigmaL) > 1e-10) {
    sigma <- (sigmaL + sigmaH)/2
    PM <- BSFormulaPut(S0, K, T, r, sigma)
    PL <- PL + (PM < P) * (PM - PL)
    sigmaL <- sigmaL + (PM < P) * (sigma - sigmaL)
    PH <- PH + (PM >= P) * (PM - PH)
    sigmaH <- sigmaH + (PM >= P) * (sigma - sigmaH)
  }
  return(sigma)
}

BSImpliedVol.OTM <- function(S0, K, T, r, V){
  f <- S0 * exp(r*T)
  res <- ifelse(K>f,BSImpliedVolCall(S0, K, T, r, V),
                BSImpliedVolPut(S0, K, T, r, V))
  return(res)
}

# Function to compute option prices and implied vols given list of final values of underlying
bsOut <- function (xf, T, AK) 
{
  nK <- length(AK)
  N <- length(xf)
  xfbar <- mean(xf)
  CAV <- numeric(nK)
  BSV <- numeric(nK)
  BSVL <- numeric(nK)
  BSVH <- numeric(nK)
  for (j in 1:nK) {
    payoff <- (xf - AK[j]) * (xf > AK[j])
    CAV[j] <- sum(payoff)/N
    err <- sqrt(var(payoff)/N)
    BSV[j] <- BSImpliedVolCall(xfbar, AK[j], T, 0, CAV[j])
    BSVL[j] <- BSImpliedVolCall(xfbar, AK[j], T, 0, CAV[j] - 
                                  err)
    BSVH[j] <- BSImpliedVolCall(xfbar, AK[j], T, 0, CAV[j] + 
                                  err)
  }
  return(list(AK=AK, CAV=CAV, BSV=BSV, BSVL=BSVL, BSVH=BSVH, err=err))
}

# Function to return implied vols for a range of strikes
analyticOut <- function(callFormula,AK,T)
{
    nK <- length(AK);
    #callFormula is a function that computes the call price
    callPrice <- numeric(nK); BSV <- numeric(nK);
    for (j in 1:nK){
                    callPrice[j] <- callFormula(AK[j]);
                    BSV[j] <- BSImpliedVolCall(1, AK[j], T,0, callPrice[j]);
                    }
    return(data.frame(AK,callPrice,BSV));
}
