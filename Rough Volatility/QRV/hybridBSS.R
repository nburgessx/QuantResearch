############################################################################
# Riemann-Liouville process (A special case of the gamma process below)
############################################################################

# The RL kernel
gRL <- function(params) function(tau){
  al <- params$al
  H <- al - 1/2
  return(sqrt(2*H)*tau^{al-1})
}

# RL variance
wRL <- function(params) function(tau){
  al <- params$al
  H <- al - 1/2
  return(tau^{2*H})
}

# RL c^star
cStarRL <- function(params) function(dt){
  al <- params$al
  H <- al - 1/2
  return(sqrt(2*H)/(H+1/2)*dt^(al-1))
}

# RL covariance
cRL <- function(params)function(s)function(u){
  
  gp <- gRL(params)
  eps <- 0
  integr <- function(r){gp(s-r)*gp(u+eps-r)}
  res <- integrate(integr, lower=0,upper=min(s,u))$value
  return(res)
  
}

# RL first order covariance
cov1RL <- function (params)function(dt){
  
  gp <- gRL(params)
  integ <- function(s){gp(s) * gp(dt+s)}
  res <- integrate(integ,lower=0,upper=dt)$value
  return(res/dt)
}

############################################################################
# Riemann-Liouville simulation
############################################################################
WtildeRL.sim <- function (params, hybrid = T) 
  function(W, Wperp) {
    library(stats)
    steps <- dim(W)[1]
    N <- dim(W)[2]
    stopifnot(dim(Wperp) == c(steps, N))
    dt <- 1/steps
    wp <- Vectorize(wRL(params))
    sqrt.dt <- sqrt(dt)
    tj <- (1:steps) * dt
    wpj <- c(0, wp(tj))
    bstar <- sqrt(diff(wpj)/dt)
    cstar <- cov1RL(params)(dt)
    rhostar <- cstar/(bstar[1] * bstar[2])
    rhobarstar <- sqrt(1 - rhostar^2)
    
    # Note that the convolution is path-by-path
    f <- function(n) {
      Wr <- W[steps:1, n]
      Y.Euler <- convolve(bstar, Wr, type = "open")[1:steps]
      Y.Correct <- bstar[1] * ((rhostar - 1) * W[, n] + rhobarstar * 
                                 Wperp[, n])
      return((Y.Correct * isTRUE(hybrid) + Y.Euler) * sqrt(dt))
    }
    Wtilde <- sapply(1:N, f)
    return(Wtilde)
  }

############################################################################
# Riemann-Liouville schemes
############################################################################
# Scheme to return S_T
######################################################
hybridSchemeRL.S <- function (params,xi) 
  function(paths, steps, expiries) {
    eta <- params$eta
    H <- params$al-1/2
    rho <- params$rho
    N <- paths
    W <- matrix(rnorm(N * steps), nrow = steps, ncol = N)
    Wperp <- matrix(rnorm(N * steps), nrow = steps, ncol = N)
    Zperp <- matrix(rnorm(N * steps), nrow = steps, ncol = N)
    Z <- rho * W + sqrt(1 - rho * rho) * Zperp
    Wtilde <- WtildeRL.sim(params)(W, Wperp)
    sim <- function(expiry) {
      dt <- expiry/steps
      ti <- (1:steps) * dt
      xi.t <- xi(ti)
      # Note the use of scaling in the below step
      v1 <- xi.t * exp(eta * expiry^H*Wtilde - 1/2 * eta^2 * ti^(2*H))
      v0 <- rep(xi(0), N)
      v <- rbind(v0, v1[-steps, ])
      logs <- apply(sqrt(v * dt) * Z - v/2 * dt, 2, sum)
      s <- exp(logs)
      return(s)
    }
    st <- t(sapply(expiries, sim))
    return(st)
  }

######################################################
# Returns matrix of y = martingale part of log S and w
######################################################
hybridSchemeRL.yw <- function (params,xi)function(paths, steps, expiries) {
  N <- paths
  eta <- params$eta
  H <- params$H
  W <- matrix(rnorm(N * steps), nrow = steps, ncol = N)
  Wperp <- matrix(rnorm(N * steps), nrow = steps, ncol = N)
  Wtilde <- WtildeRL.sim(params)(W, Wperp)
  yw <- function(expiry) {
    dt <- expiry/steps
    ti <- (1:steps) * dt
    xi.t <- xi(ti)
    v1 <- xi.t * exp(eta * expiry^H * Wtilde - 1/2 * eta^2 * ti^(2 * H))
    v0 <- rep(xi(0), N)
    v <- rbind(v0, v1[-steps, ])
    y <- apply(sqrt(v * dt) * W, 2, sum)
    w <- apply((v + v1)/2 * dt, 2, sum)
    
    # Now do some variance reduction
    #w.exact <- integrate(xi,lower=0,upper=expiry)$value
    #w <- w + w.exact - mean(w) # Control variance E[w]
    #y <- y - mean(y) # E[Y]=0
    
    yw.res <- array(dim = c(2, N))
    yw.res[1, ] <- y
    yw.res[2, ] <- w
    
    return(yw.res)
  }
  n.exp <- length(expiries)
  ywMatrix <- matrix(NA, nrow = 2 * n.exp, ncol = N)
  for (i in 1:n.exp) {
    ywMatrix[(1:2) + 2 * (i - 1), ] <- yw(expiries[i])
  }
  return(ywMatrix)
}

############################################################################
# Gamma process
############################################################################

# Gamma kernel
gGamma <- function(params) function(tau){
  al <- params$al
  H <- al - 1/2
  lam <- params$lam
  
  return(sqrt(2*H)*tau^{al-1}*exp(-lam*tau))
}

# Gamma variance
wGamma <- function(params) function(tau){
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


# Gmma c^star
cStarGamma <- function(params) function(dt){
  
  al <- params$al
  H <- al - 1/2
  lam <- params$lam
  
  prefactor <- sqrt(2*H)/(lam^al*dt)
  bkt <- gamma(al)- gamma_inc(al,lam*dt)
  res2 <- sqrt(2*H)/al*dt^(al-1)
  res <- ifelse(lam>0,prefactor * bkt,res2)
  return(res)
  
}

# Gamma covariance
cGamma <- function(params)function(t)function(u){
  
  gp <- gGamma(params)
  eps <- 0
  integr <- function(s){gp(t-s)*gp(u+eps-s)}
  res <- integrate(integr, lower=0,upper=min(t,u))$value
  return(res)
  
}

# Gamma first order covariance
cov1Gamma <- function (params)function(dt){
  
  gp <- gGamma(params)
  integ <- function(s){gp(s) * gp(dt+s)}
  res <- integrate(integ,lower=0,upper=dt)$value
  return(res/dt)
}

############################################################################
# Gamma process simulation (note the extra expiry argument)
############################################################################
# Scheme to return S_T
######################################################
WtildeGamma.sim <- function (params, hybrid=T)function(W, Wperp)function(expiry) 
{
  library(stats)
  library(gsl)
  steps <- dim(W)[1]
  N <- dim(W)[2]
  stopifnot(dim(Wperp) == c(steps, N))
  
  dt <- expiry/steps
  wp <- Vectorize(wGamma(params))
  sqrt.dt <- sqrt(dt)
  tj <- (1:steps) * dt
  wpj <- c(0, wp(tj))
  bstar <- sqrt(diff(wpj)/dt)
  cstar <- cov1Gamma(params)(dt)
  rhostar <- cstar/(bstar[1] * bstar[2])
  rhobarstar <- sqrt(1 - rhostar^2)
  
  # Note that the convolution is path-by-path
  f <- function(n) {
    Wr <- W[steps:1, n]
    Y.Euler <- convolve(bstar, Wr, type = "open")[1:steps]
    Y.Correct <- bstar[1] * ((rhostar - 1) * W[, n] + rhobarstar * 
                               Wperp[, n])
    return((Y.Correct * isTRUE(hybrid) + Y.Euler) * sqrt(dt))
  }
  Wtilde <- sapply(1:N, f)
  return(Wtilde)
}

############################################################################
# Gamma hybrid schemes
############################################################################
hybridSchemeGamma.S <- function (params,xi) 
  function(paths, steps, expiries) {
    eta <- params$eta
    H <- params$al-1/2
    rho <- params$rho
    N <- paths
    W <- matrix(rnorm(N * steps), nrow = steps, ncol = N)
    Wperp <- matrix(rnorm(N * steps), nrow = steps, ncol = N)
    Zperp <- matrix(rnorm(N * steps), nrow = steps, ncol = N)
    Z <- rho * W + sqrt(1 - rho * rho) * Zperp
    Wtilde.e <- WtildeGamma.sim(params)(W, Wperp) # This is a function of expiry
    wGamma.p <- Vectorize(wGamma(params))
    sim <- function(expiry) {
      dt <- expiry/steps
      ti <- (1:steps) * dt
      xi.t <- xi(ti)
      v1 <- xi.t * exp(eta * Wtilde.e(expiry) - 1/2 * eta^2 * wGamma.p(ti))
      v0 <- rep(xi(0), N)
      v <- rbind(v0, v1[-steps, ])
      logs <- apply(sqrt(v * dt) * Z - v/2 * dt, 2, sum)
      s <- exp(logs)
      return(s)
    }
    st <- t(sapply(expiries, S))
    return(st)
  }

######################################################
# Returns matrix of y = martingale part of log S and w
######################################################
hybridSchemeGamma.yw <- function (params,xi)function(paths, steps, expiries) {
  N <- paths
  eta <- params$eta
  H <- params$H
  W <- matrix(rnorm(N * steps), nrow = steps, ncol = N)
  Wperp <- matrix(rnorm(N * steps), nrow = steps, ncol = N)
  yw.res <- array(dim = c(2, N))
  Wtilde.e <- WtildeGamma.sim(params)(W, Wperp) # This is a function of expiry
  wGamma.p <- Vectorize(wGamma(params))
  yw <- function(expiry) {
    dt <- expiry/steps
    ti <- (1:steps) * dt
    xi.t <- xi(ti)
    v1 <- xi.t * exp(eta * Wtilde.e(expiry) - 1/2 * eta^2 * wGamma.p(ti))
    #print(mean(v1))
    v0 <- rep(xi(0), N)
    v <- rbind(v0, v1[-steps, ])
    y <- apply(sqrt(v * dt) * W, 2, sum)
    w <- apply((v + v1)/2 * dt, 2, sum)
    
    # Now do some variance reduction
    #w.exact <- integrate(xi,lower=0,upper=expiry)$value
    #w <- w + w.exact - mean(w) # Control variance E[w]
    #y <- y - mean(y) # E[Y]=0
    
    yw.res[1, ] <- y
    yw.res[2, ] <- w
    
    return(yw.res)
  }
  n.exp <- length(expiries)
  ywMatrix <- matrix(NA, nrow = 2 * n.exp, ncol = N)
  for (i in 1:n.exp) {
    ywMatrix[(1:2) + 2 * (i - 1), ] <- yw(expiries[i])
  }
  return(ywMatrix)
}
