####################################################################################
# Jim Gatheral, 2023
#
# The following functions use the Lewis representation of the option price
# in terms of the characteristic function, equation (5.6) of 
# The Volatility Surface.
# 
# Integration and vectorization improved in 2023.
####################################################################################

source("BlackFormula.R")

option.OTM.raw <- function (phi, k, tau) 
{
  integrand <- function(u) {
    Re(exp(-(0+1i) * u * k) * phi(u - (0+1i)/2, tau)/(u^2 + 
                                                      1/4))
  }
  k.minus <- (k < 0) * k
  res <- exp(k.minus) - exp(k/2)/pi * 
    integrate(integrand,lower=0,upper=Inf,rel.tol=1e-10,subdivisions=1000)$value 
  return(ifelse(res<0,NA,res))
}

option.OTM <- Vectorize(option.OTM.raw,vectorize.args=c("k","tau"))

impvol.phi <- function(phi)function(k, t) {
  optVal <- option.OTM(phi, k, t)
  res <- ifelse(is.na(optVal),NA,
                ivOTM(1, exp(k), t, optVal))
  return(res)
}

####################################################################################
# Function to compute the ATM skew, equation (5.8) of 
# The Volatility Surface.
####################################################################################
atmSkew <- function(phi)Vectorize(function(T){
  atmVol <- impvol.phi(phi)(0,T)
  integrand <-  function(u){Im(u*phi(u - 1i/2, T)/(u^2 + 1/4))}
  res <- -integrate(integrand,lower=0,upper=Inf,rel.tol=1e-10,subdivisions=1000)$value/
    sqrt(T)*sqrt(2/pi)*exp(atmVol^2*T/8)
  return(res)})