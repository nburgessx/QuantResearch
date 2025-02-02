# We construct a piecewise constant forward variance curve from variance swap data

# xi.curve <- function(expiries,w) function (t) 
# {
#   n <- length(w)
#   xi.vec <- c(w[1]/expiries[1], diff(w)/diff(expiries))
#   exp.vec <- c(0,expiries)
#   exp.vec[n+1] <- Inf # Long expiries are all in the last bucket
#   pos <- findInterval(t, vec = exp.vec)
#   res <- xi.vec[pos]
#   return(res)
# }

obj.w <- function(expiries,w.in)function(err.vec){
  
  w.in.1 <- w.in + 2*sqrt(w.in)*err.vec*sqrt(expiries)
  
  xi.vec <- c(w.in.1[1]/expiries[1],diff(w.in.1)/diff(expiries)) 
  
  dxi.dt <- diff(xi.vec)/diff(expiries)
  w.out <- c(0,cumsum(xi.vec[-1]*diff(expiries)))+xi.vec[1]*expiries[1]
  
  res <- sum((w.in-w.out)^2) + sum(dxi.dt^2)
  return(res*1e3)
}

xi.curve <- function(expiries,w.in,eps=0){
  
  n <- length(w.in)
  
  if(eps>0){
    res.optim <- optim(rep(0,n),obj.w(expiries,w.in),method="L-BFGS-B",
                     lower=rep(-eps,n),upper=rep(eps,n))
    err.vec <- res.optim$par
    w.in.1 <- w.in  + 2*sqrt(w.in)*err.vec*sqrt(expiries)}
  else {w.in.1 <- w.in}
  xi.vec.out <- c(w.in.1[1]/expiries[1],diff(w.in.1)/diff(expiries))
  xi.curve.raw <- function(t){
    res <- ifelse(t <= expiries[n], xi.vec.out[sum(expiries < t)+1],xi.vec.out[n])
    return(res)
  }
  xi.curve.out <- Vectorize(xi.curve.raw)
  fit.errs <- sqrt(w.in.1/expiries) - sqrt(w.in/expiries)
  
  return(list(xi.vec=xi.vec.out,xi.curve=xi.curve.out,fit.errs=fit.errs,w.out=w.in.1))
  
}

############################################################################
# Code due to Rick Cao September 2019
############################################################################
xi.curve.smooth <- function(expiries, w.in, xi=TRUE,eps=0){
  
  # Auxiliary functions
  phi <- function(tau)function(x) {
    min = min(x, tau)   
    return (1 - min^3 / 6 + x * tau * (2 + min) / 2)
  }
  phi.deri <- function(tau)function(x) {
    min = min(x, tau)   
    return (tau - min^2 / 2 + tau * min)
  }
  
  n <- length(expiries)
  c <- diag(n)
  
  A <- sapply(expiries, phi(expiries[1]))
  for (i in seq(2, n)) {
    A <- rbind(A, sapply(expiries, phi(expiries[i])))
  }
  
  obj.1 <- function(err.vec) {
    v <- w.in + 2*sqrt(w.in)*err.vec*sqrt(expiries)
    return (t(v) %*% solve(c %*% A %*% t(c)) %*% v)
  }
  
  res.optim <- optim(rep(0,n),obj.1,method="L-BFGS-B",
                     lower=rep(-eps,n),upper=rep(eps,n))
  err.vec <- res.optim$par
  
  w.in.1 <- w.in  + 2*sqrt(w.in)*err.vec*sqrt(expiries)
  
  Z <- t(c) %*% solve(c %*% A %*% t(c)) %*% w.in.1
  
  curve.raw <- function(x){ 
    
    sum.curve <- 0
    sum.curve.deri <- 0
    
    for (i in seq(1, n)) {
      sum.curve <- sum.curve + Z[i] * phi(expiries[i])(x)
      sum.curve.deri <- sum.curve.deri + Z[i] * phi.deri(expiries[i])(x)
    }
    if (xi) {
      return (sum.curve.deri)
    } else {
      return (sum.curve)
    }
  }
  
  xi.curve.out <- Vectorize(curve.raw)
  
  fit.errs <- sqrt(w.in.1/expiries) - sqrt(w.in/expiries)
  
  return(list(xi.curve=xi.curve.out, fit.errs=fit.errs, w.out=w.in.1))
  
}

CurveSmoothBuilder <- function(expiries, w.in, eps=0)
{

      ## Get the environment for this instance of the function.
      thisEnv <- environment()

      expiries <- expiries
      w.in <- w.in
      w.out <- w.in
      eps <- eps
      fit.errs <- 0
      
#       phi <- function(x)function(tau) {
#             min = min(x, tau)   
#             return (1 - min^3 / 6 + x * tau * (2 + min) / 2)
#       }
          
#       phi.deri <- function(x)function(tau) {
#             min = min(x, tau)   
#             return (tau - min^2 / 2 + tau * min)
#       }

      ## Create the list used to represent an object for this class
      me <- list(

              ## Define the environment where this list is defined so
              ## that I can refer to it later.
              thisEnv = thisEnv,

              fitCurve = function()
              {
                      expiries <- get("expiries",thisEnv)
                      w.in <- get("w.in",thisEnv)
                      eps <- get("eps",thisEnv)
                      
                      # Auxiliary functions
                      phi <- function(tau)function(x) {
                            min = min(x, tau)   
                            return (1 - min^3 / 6 + x * tau * (2 + min) / 2)
                      }

                      n <- length(expiries)
                      c <- diag(n)

                      A <- sapply(expiries, function(t){sapply(expiries, phi(t))})

                      obj <- function(err.vec) {
                        v <- w.in + 2*sqrt(w.in)*err.vec*sqrt(expiries)
                        return (t(v) %*% solve(c %*% A %*% t(c)) %*% v)
                      }

                      res.optim <- optim(rep(0,n),obj,method="L-BFGS-B",
                                         lower=rep(-eps,n),upper=rep(eps,n))
                      err.vec <- res.optim$par

                      w.out <- w.in  + 2*sqrt(w.in)*err.vec*sqrt(expiries)

                      Z.raw <- t(c) %*% solve(c %*% A %*% t(c)) %*% w.out
                      
                      fit.errs <- sqrt(w.out/expiries) - sqrt(w.in/expiries)

                      Z <- c(0, Z.raw, 0)
                      ts <- c(0, expiries, 0)
                      phi.list.raw <- sapply(expiries, function(tau){1 - tau^3 / 6})
                      phi.deri.list.raw <- sapply(expiries, function(tau){tau + tau^2 / 2})
                      phi.list <- c(0, phi.list.raw, 0)
                      phi.deri.list <- c(0, phi.deri.list.raw, 0)
                      cumsum.Z.phi.list <- cumsum(Z*phi.list)
                      cumsum.Z.phi.deri.list <- cumsum(Z*phi.deri.list)
                      revcumsum.Z.t.list <- cumsum(rev(Z*ts))
                      revcumsum.Z.list <- cumsum(rev(Z))
                      
                      assign("Z",Z.raw,thisEnv)
                      assign("w.out", w.out, thisEnv)
                      assign("fit.errs", fit.errs, thisEnv)
                      assign("cumsum.Z.phi.list",cumsum.Z.phi.list,thisEnv)
                      assign("cumsum.Z.phi.deri.list",cumsum.Z.phi.deri.list,thisEnv)
                      assign("revcumsum.Z.t.list",revcumsum.Z.t.list,thisEnv)
                      assign("revcumsum.Z.list",revcumsum.Z.list,thisEnv)
                      
              },

              getTotalVarCurve = function(x)
              {
                      expiries <- get("expiries",thisEnv)
                      cumsum.Z.phi.list <- get("cumsum.Z.phi.list",thisEnv)
                      cumsum.Z.phi.deri.list <- get("cumsum.Z.phi.deri.list",thisEnv)
                      revcumsum.Z.t.list <- get("revcumsum.Z.t.list",thisEnv)
                      revcumsum.Z.list <- get("revcumsum.Z.list",thisEnv)
                      curve <- function(x) {
                            a <- sum(expiries <= x) + 1
                            b <- sum(expiries > x) + 1
                            parta <- cumsum.Z.phi.list[a] + x * cumsum.Z.phi.deri.list[a]
                            partb <- revcumsum.Z.list[b] * (1 - (x^3) / 6) + revcumsum.Z.t.list[b] * (x + (x^2) / 2)
                            curve <- parta + partb
                            return(curve)
                      }
                      res <- Vectorize(curve)
                      return(res)
              },
          
              getForwardVarCurve = function()
              {
                      expiries <- get("expiries",thisEnv)
                      cumsum.Z.phi.deri.list <- get("cumsum.Z.phi.deri.list",thisEnv)
                      revcumsum.Z.t.list <- get("revcumsum.Z.t.list",thisEnv)
                      revcumsum.Z.list <- get("revcumsum.Z.list",thisEnv)
                      curve <- function(x) {
                            a <- sum(expiries <= x) + 1
                            b <- sum(expiries > x) + 1
                            curve <- cumsum.Z.phi.deri.list[a] + revcumsum.Z.t.list[b] * (1 + x) - (x^2 / 2) * revcumsum.Z.list[b]
                            return(curve)
                      }
                      res <- Vectorize(curve)
                      return(res)
              },
          
              getFitErrors = function() 
              {
                  return(get("fit.errs",thisEnv))
              },
              
              getWOut = function()
              {
                  return(get("w.out",thisEnv))
              }
              
        )

      ## Define the value of the list within the current environment.
      assign('this',me,envir=thisEnv)

      ## Set the name for the class
      class(me) <- append(class(me),"CurveSmoothBuilder")
      return(me)
}