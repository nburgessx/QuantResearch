#####################################################################################
# Jim Gatheral, June 2021
# Negative H added May 2023
#####################################################################################

library(gsl)
library(MittagLeffleR)

#######################################################
# Gamma kernel functions
#######################################################

# Gamma kernel
gGamma <- function(params) function(tau){
  al <- params$al
  H <- al - 1/2
  lam <- params$lam
  return(sqrt(2*H)*tau^{al-1}*exp(-lam*tau))
}

# Gamma variance G00
G00 <- function(params) function(tau){
  al <- params$al
  H <- al - 1/2
  H2 <- 2*H
  lam <- params$lam
  
  prefactor <- H2/((2*lam)^H2)
  bkt <- gamma(H2)- gamma_inc(H2,2*lam*tau)
  res2 <- tau^(2*H)
  res <- ifelse(lam>0,prefactor * bkt,res2)
  return(res)
}

# Gamma variance K00
K00 <- function(params) function(tau){params$eta^2*Vectorize(G00(params))(tau)}

# Gamma average K0
K0 <- function(params) function(tau){params$eta*Vectorize(G0(params))(tau)}

# G11
G11 <- function(params) function(tau){
  al <- params$al
  H <- al - 1/2
  H2 <- 2*H
  lam <- params$lam
  
  prefactor <- H2/((2*lam)^H2)
  bkt <- gamma_inc(H2,2*lam*tau)- gamma_inc(H2,4*lam*tau)
  res2 <- tau^(2*H)*(2^H2-1)
  res <- ifelse(lam>0,prefactor * bkt,res2)
  return(res)
}

# Gkk
Gkk <- function(params)function(tau) Vectorize(
  function(k) {
    al <- params$al
    H <- al - 1/2
    lam <- params$lam
    gGamma2 <-function(s){(gGamma(params)(s+k*tau))^2}
    res <- integrate(gGamma2,lower=0,upper=tau)$value
    return(res)
  })


# G0
G0 <- function(params) function(dt){
  
  al <- params$al
  H <- al - 1/2
  lam <- params$lam
  
  prefactor <- sqrt(2*H)/(lam^al)
  bkt <- gamma(al)- gamma_inc(al,lam*dt)
  res2 <- sqrt(2*H)/al*dt^(al)
  res <- ifelse(lam>0,prefactor * bkt,res2)
  return(res)
  
}

G1 <- function(params) function(dt){
  
  al <- params$al
  H <- al - 1/2
  lam <- params$lam
  
  prefactor <- sqrt(2*H)/(lam^al)
  bkt <- gamma_inc(al,lam*dt)- gamma_inc(al,2*lam*dt)
  res2 <- sqrt(2*H)/al*dt^(al)*(2^al-1)
  res <- ifelse(lam>0,prefactor * bkt,res2)
  return(res)
  
}

# Gamma covariance
G0k <- function(params,k)function(t){
  
  gp <- gGamma(params)
  eps <- 0
  integr <- function(s){gp(s)*gp(s+k*t)}
  res <- integrate(integr, lower=0,upper=t)$value
  return(res)
  
}

# Gamma first order covariance
G01 <- function(params)function(t){
  
  gp <- gGamma(params)
  eps <- 0
  integr <- function(s){gp(s)*gp(s+t)}
  res <- integrate(integr, lower=0,upper=t)$value
  return(res)
  
}

# Resolvent kernel of gGamma^2
bigK <- function(params)function(tau){
  al <- params$al
  H <- al - 1/2
  H.2 <- 2*al-1
  lam <- params$lam
  eta <- params$eta
  etaHat2 <- eta^2*H.2*gamma(H.2)
  tau.2H <- tau^(H.2)
  res <- etaHat2*exp(-2*lam*tau) * tau^(H.2-1)*mlf(etaHat2*tau.2H,H.2,H.2)
  return(res)    
}

# K0
bigK0.raw <- function(params)function(tau){
  
  bigKp <- bigK(params)
  integ <- function(s){bigKp(s)}
  res <- ifelse(tau>0,integrate(integ,lower=0,upper=tau)$value,0)
  return(res)    
}

bigK0 <- function(params){Vectorize(bigK0.raw(params))}
  
########################################################################
# Discrete version of resolvemt kernel computation
K.matrix <- function(params,tau)function(n){
  res <- array(0,dim=c(n,n))
  del <- tau/n
  for (i in 1:n){
    for (j in 1:i){
      res[i,j] <- K00(params)(del*(i-j+1))-K00(params)(del*(i-j))
    }
  }
  return(res)
}

# Resolvent kernel matrix
Ktilde.matrix <- function(K.matrix){
  n <- dim(K.matrix)[1]
  res <- solve(diag(1,n)-K.matrix) %*% K.matrix
  return(res)
}

XdmM <- function(params,ey,xi)Vectorize(function(T){
  
  del <- 1/12
  ybar <- ey(T+del/2)
  yFactor <- ybar/(ybar^2+params$c)
  
  c <- params$c
  del <- 1/12
  ybar <- ey(T+del/2)
  yFactor <- ybar/(ybar^2+c)
  kappaP <- function(tau){params$eta*Vectorize(gGamma(params))(tau)}
  
  K0p <- K0(params)
  bigK0p <- bigK0(params)
  bigKMp <- function(tau){1+ bigK0p(tau)}
  bigA <- function(s){1/del*(K0p(T+del-s)-K0p(T-s))}
  bigB <- function(r)function(s){
    integ <- function(u){1/del*kappaP(u-r)*kappaP(u-s)}
    return(integrate(integ, lower=T, upper=T+del)$value)
  }
  
  DF <- function(s){-2*yFactor*bigA(s)}
  
  integ.s <- function(r)function(s){sqrt(xi(s))*DF(s)*(
    -ey(s)*kappaP(s-r)*bigKMp(s-r)*DF(s) +
      2*xi(r)/(ybar^2+c)*bigB(r)(s) -
      sqrt(xi(r))*sqrt(xi(s))*DF(r)*DF(s)
  )}
  
  int.s <- Vectorize(function(r){integrate(integ.s(r),lower=r,upper=T)$value})
  #     return(int.s(.03))         
  integ.r <- function(r){ sqrt(xi(r))*DF(r)*int.s(r)}
  int.r <- integrate(integ.r,lower=0,upper=T)$value
  
  res <- int.r/4
  return(res) 
})
