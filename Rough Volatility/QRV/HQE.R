###############################################################
#
# RSQE and HQE codes
# Jim Gatheral, September 2022.
#
###############################################################

psiM <- function(psi,ev,w){
  beta2 <- 2/psi-1+sqrt(2/psi)*sqrt(abs(2/psi-1)) # The abs fixes situations where psi > 2
  alpha <- ev/(1+beta2)
  vf <- alpha*(sqrt(abs(beta2))+w)^2
  return(vf)
}

psiP <- function(psi,ev,u){
  p <- 2/(1+psi)
  gam <- ev/2*(1+psi)
  vf <- -(u<p)*gam*log(u/p)
  return(vf)
}

RSQE.sim <- function(params, xi)function(T, paths, steps){
  library(gsl)
  eta <- params$eta
  lam <- params$lambda
  H <- params$al - 1/2
  rho <- params$rho
  rho2m1 <- sqrt(1 - rho * rho)
  eps.0 <- 1e-10
  
  W <- matrix(rnorm(steps * paths), nrow = steps, ncol = paths)
  Wperp <- matrix(rnorm(steps * paths), nrow = steps, ncol = paths)
  U <- matrix(runif(steps * paths), nrow = steps, ncol = paths)
  
  G00p <- Vectorize(G00(params))
  
  dt <- T/steps
  sqrt.dt <- sqrt(dt)
  tj <- (1:steps) * dt
  xij <- xi(tj)
  G00del <- G00(params)(dt)
  G00j <- c(0, G00p(tj))
  bstar <- sqrt(diff(G00j)/dt)
  bstar1 <- bstar[1]
  u <- array(0, dim = c(steps, paths))
  v <- rep(xi(0), paths)
  xihat <- rep(xij[1], paths)
  x <- numeric(paths)
  y <- numeric(paths)
  w <- numeric(paths)
  
  for (j in 1:steps) {
    varv <- eta^2 * (xihat + 2 * H * v)/(1 + 2 * H) * 
      G00del
    psi <- varv/xihat^2
    vf <- ifelse(psi < 3/2, psiM(psi, xihat, W[j, ]), 
                 psiP(psi, xihat, U[j, ]))
    u[j, ] <- vf - xihat
    dw <- (v + vf)/2 * dt
    w <- w + dw
    dy <- as.numeric(u[j, ])/(eta * bstar1)
    y <- y + dy
    x <- x - dw/2 + sqrt(dw) * as.numeric(rho2m1 * Wperp[j, 
    ]) + rho * dy
    btilde <- rev(bstar[2:(j + 1)])
    if (j < steps) {  # This step can be much faster if H=1/2
      xihat <- xij[j + 1] + as.numeric(btilde %*% u[1:j, 
      ])/bstar1  
    }
    xihat <- ifelse(xihat > eps.0, xihat, eps.0)
    v <- vf
  }
  
  res <- list(x=x,v=v,w=w)
  return(res)
}

HQE.sim <- function (params, xi)function(T, paths, steps) {
  
  library(gsl)
  nu <- params$eta
  lam <- params$lambda
  H <- params$al - 1/2
  rho <- params$rho
  rho2m1 <- sqrt(1 - rho * rho)
  eps.0 <- 1e-10
  
  W <- matrix(rnorm(steps * paths), nrow = steps, ncol = paths)
  Wperp <- matrix(rnorm(steps * paths), nrow = steps, ncol = paths)
  Z <- matrix(rnorm(steps * paths), nrow = steps, ncol = paths)
  U <- matrix(runif(steps * paths), nrow = steps, ncol = paths)
  Uperp <- matrix(runif(steps * paths), nrow = steps, ncol = paths)
  
  dt <- T/steps
  sqrt.dt <- sqrt(dt)
  tj <- (1:steps) * dt
  xij <- xi(tj)
  G0del <- nu * G0(params)(dt)
  G1del <- nu * G1(params)(dt)
  G01del <- nu^2 * G01(params)(dt)
  Gjj <- nu^2 * Gkk(params)(dt)((1:steps)-1)
  G00del <- Gjj[1]
  G11del <- Gjj[2]
  bstar <- sqrt(Gjj/dt)
  bstar1 <- bstar[1]
  rho.vchi <- G0del/sqrt(G00del * dt)
  beta.vchi <- G0del/dt
  u <- array(0, dim = c(steps, paths))
  chi <- array(0, dim = c(steps, paths))
  v <- rep(xi(0), paths)
  xihat <- rep(xij[1], paths)
  x <- numeric(paths)
  y <- numeric(paths)
  w <- numeric(paths)
  
  for (j in 1:steps) {
    xibar <- (xihat + 2 * H * v)/(1 + 2 * H)
    var.eps <- xibar * G00del * (1 - rho.vchi^2)
    psi.chi <- 4 * G00del * rho.vchi^2 * xibar/xihat^2
    psi.eps <- 4 * G00del * (1 - rho.vchi^2) * xibar/xihat^2
    z.chi <- ifelse(psi.chi < 3/2, psiM(psi.chi, xihat/2, 
                                        W[j, ]), psiP(psi.chi, xihat/2, U[j, ]))
    z.eps <- ifelse(psi.eps < 3/2, psiM(psi.eps, xihat/2, 
                                        Wperp[j, ]), psiP(psi.eps, xihat/2, Uperp[j, 
                                        ]))
    chi[j, ] <- (z.chi - xihat/2)/beta.vchi
    eps <- z.eps - xihat/2
    u[j, ] <- beta.vchi * chi[j, ] + eps
    vf <- xihat + u[j, ]
    vf <- ifelse(vf > eps.0, vf, eps.0)
    dw <- (v + vf)/2 * dt
    w <- w + dw
    y <- y + chi[j, ]
    x <- x - dw/2 + sqrt(dw) * as.numeric(rho2m1 * Z[j, 
    ]) + rho * chi[j, ]
    btilde <- rev(bstar[2:(j + 1)])
    if (j < steps) {
      xihat <- xij[j + 1] + as.numeric(btilde %*% chi[1:j, 
      ])
    }
    v <- vf
  }
  res <- list(x=x,v=v,w=w)
  return(res)
}


