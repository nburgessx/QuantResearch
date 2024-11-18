# We rewrite plotIvols to take an array st of stock prices for each expiry.
# In this version, each row of the array has to correspond to an option expiry.

library(stinepack);
source("BlackScholes.R")


plotIvolsMC <- function (ivolData, sviMatrix = NULL, slices = NULL, mcMatrix = NULL, 
                         plot = TRUE, colnum = NULL) 
{
  bidVols <- as.numeric(ivolData$Bid)
  askVols <- as.numeric(ivolData$Ask)
  expDates <- unique(ivolData$Texp)
  nSlices <- length(expDates)
  if (!is.null(slices)) {
    nSlices <- length(slices)
  }
  else {
    slices <- 1:nSlices
  }
  if(is.null(colnum)) {colnum <- sqrt(nSlices * 2)}
  rows <- round(colnum/2, 0)
  columns <- round(colnum, 0)
  while (rows * columns < nSlices) {
    rows <- rows + 1
  }
  atmVols <- numeric(nSlices)
  atmSkew <- numeric(nSlices)
  atmCurv <- numeric(nSlices)
  atmVolsMC <- numeric(nSlices)
  atmSkewMC <- numeric(nSlices)
  atmCurvMC <- numeric(nSlices)
  atmErrMC <- numeric(nSlices)
  par(mfrow = c(rows, columns), mex = 0.5)
  for (slice in slices) {
    t <- expDates[slice]
    texp <- ivolData$Texp
    bidVol <- bidVols[texp == t]
    askVol <- askVols[texp == t]
    midVol <- (bidVol + askVol)/2
    f <- (ivolData$Fwd[texp == t])[1]
    k <- log(ivolData$Strike[texp == t]/f)
    include <- !is.na(bidVol) & (bidVol > 0)
    kmin <- min(k[include])
    kmax <- max(k[include])
    ybottom <- 0.6 * min(bidVol[include])
    ytop <- 1.2 * max(askVol[include], na.rm = T)
    xrange <- c(kmin, kmax)
    yrange <- c(ybottom, ytop)
    if (plot == TRUE) { plot(k, bidVol, col = "red2", pch = 24, cex = 0.5, 
                             xlim = xrange, ylim = yrange, main = paste("T =", 
                                                                        format(t, digits = 2, nsmall = 2)), xlab = "Log-Strike", 
                             ylab = "Implied Vol.")
                        points(k, askVol, col = "royalblue", pch = 25, cex = 0.5)
                        abline(v=0, lty=2)
                        abline(h=0, lty=2) }
    if ((!is.null(sviMatrix))) {
      vol <- function(k) {
        sqrt(svi(sviMatrix[slice, ], k)/t)
      }
      if (plot == TRUE) {  curve(vol(x), from = kmin, to = kmax, col = "green4", lwd = 2, add = T)}
    }
    if ((!is.null(mcMatrix))) {
      spots <- mcMatrix[slice, ]
      s0 <- mean(spots)
      volMC <- function(k) {
        bsOut(spots, t, s0 * exp(k))$BSV
      }
      atmVolsMC[slice] <- volMC(0)
      sigMC <- volMC(0) * sqrt(t) 
      atmSkewMC[slice] <- (volMC(sigMC/10) - volMC(-sigMC/10))/(2*sigMC/10)
      atmCurv[slice] <- (volMC(sigMC/10)+volMC(-sigMC/10)-2*volMC(0))/(2*(sigMC/10)^2)
      errMC <- function(k) {
        bsOut(spots, t, s0 * exp(k))$err
      }
      atmErrMC[slice] <- errMC(0)
      if (plot == TRUE) {curve(volMC(x), from = kmin, to = kmax, col = "green4", lwd=2, add=T)}
    }
  kIn <- k[!is.na(midVol)]
  volIn <- midVol[!is.na(midVol)]
  volInterp <- function(xout) {
    stinterp(x = kIn, y = volIn, xout)$y
  }
  atmVols[slice] <- volInterp(0)
  sig <- atmVols[slice] * sqrt(t) 
  atmSkew[slice] <- (volInterp(sig/10)-volInterp(-sig/10))/(2*sig/10)
  atmCurv[slice] <- (volInterp(sig/10)+volInterp(-sig/10)-2*atmVols[slice])/(2*(sig/10)^2)
  
}
  par(mfrow = c(1, 1), mex = 1)
  par(new = F)
  return(list(expiries = expDates, atmVols = atmVols, atmSkew = atmSkew,atmCurv=atmCurv, 
              atmVolsMC = atmVolsMC, atmSkewMC = atmSkewMC, atmErrMC = atmErrMC))
}

#######################################################################
# Version that takes two matrices
plotIvolsMC2 <- function (ivolData, slices = NULL, mcMatrix = NULL, mcMatrix2 = NULL, 
                          plot = TRUE, rowcol = NULL) 
{
  bidVols <- as.numeric(ivolData$Bid)
  askVols <- as.numeric(ivolData$Ask)
  expDates <- unique(ivolData$Texp)
  nSlices <- length(expDates)
  gr <- (1 + sqrt(5))/2
  if (!is.null(slices)) {
    nSlices <- length(slices)
  }
  else {
    slices <- 1:nSlices
  }
  gr <- (1 + sqrt(5))/2
  colnum <- sqrt(nSlices * gr)
  rows <- round(colnum/2, 0)
  columns <- round(colnum, 0)
  
  while (rows * columns < nSlices) {
    rows <- rows + 1
  }
  if (!is.null(rowcol)){
    rows <- rowcol[1]
    columns <- rowcol[2]
  }
  atmVols <- numeric(nSlices)
  atmSkew <- numeric(nSlices)
  atmCurv <- numeric(nSlices)
  atmVolsMC <- numeric(nSlices)
  atmSkewMC <- numeric(nSlices)
  atmCurvMC <- numeric(nSlices)
  atmErrMC <- numeric(nSlices)
  atmVolsMC2 <- numeric(nSlices)
  atmSkewMC2 <- numeric(nSlices)
  atmCurvMC2 <- numeric(nSlices)
  atmErrMC2 <- numeric(nSlices)
  par(mfrow = c(rows, columns), mex = 0.5)
  for (slice in slices) {
    t <- expDates[slice]
    texp <- ivolData$Texp
    bidVol <- bidVols[texp == t]
    askVol <- askVols[texp == t]
    midVol <- (bidVol + askVol)/2
    f <- (ivolData$Fwd[texp == t])[1]
    k <- log(ivolData$Strike[texp == t]/f)
    include <- !is.na(bidVol) & (bidVol > 0)
    kmin <- min(k[include])
    kmax <- max(k[include])
    ybottom <- 0.6 * min(bidVol[include])
    ytop <- 1.2 * max(askVol[include], na.rm = T)
    xrange <- c(kmin, kmax)
    yrange <- c(ybottom, ytop)
    if (plot == TRUE) {
      plot(k, bidVol, col = "red1", pch = 24, cex = 0.5, 
           xlim = xrange, ylim = yrange, main = paste("T =", 
                                                      format(t, digits = 2, nsmall = 2)), xlab = "Log-Strike", 
           ylab = "Implied Vol.")
      points(k, askVol, col = "royalblue", pch = 25, cex = 0.5)
      abline(v = 0, lty = 2)
      abline(h = 0, lty = 2)
    }
    if ((!is.null(mcMatrix))) {
      spots <- mcMatrix[slice, ]
      s0 <- mean(spots)
      volMC <- function(k) {
        bsOut(spots, t, s0 * exp(k))$BSV
      }
      atmVolsMC[slice] <- volMC(0)
      sigMC <- volMC(0) * sqrt(t)
      atmSkewMC[slice] <- (volMC(sigMC/10) - volMC(-sigMC/10))/(2 * 
                                                                  sigMC/10)
      atmCurv[slice] <- (volMC(sigMC/10) + volMC(-sigMC/10) - 
                           2 * volMC(0))/(2 * (sigMC/10)^2)
      errMC <- function(k) {
        bsOut(spots, t, s0 * exp(k))$err
      }
      atmErrMC[slice] <- errMC(0)
      if (plot == TRUE) {
        curve(volMC(x), from = kmin, to = kmax, col = "green4", 
              lwd = 2, add = T)
      }
    }
    if ((!is.null(mcMatrix2))) {
      spots <- mcMatrix2[slice, ]
      s0 <- mean(spots)
      volMC <- function(k) {
        bsOut(spots, t, s0 * exp(k))$BSV
      }
      atmVolsMC2[slice] <- volMC(0)
      sigMC2 <- volMC(0) * sqrt(t)
      atmSkewMC2[slice] <- (volMC(sigMC2/10) - volMC(-sigMC2/10))/(2 * 
                                                                     sigMC2/10)
      atmCurvMC2[slice] <- (volMC(sigMC2/10) + volMC(-sigMC2/10) - 
                              2 * volMC(0))/(2 * (sigMC2/10)^2)
      errMC2 <- function(k) {
        bsOut(spots, t, s0 * exp(k))$err
      }
      atmErrMC2[slice] <- errMC2(0)
      if (plot == TRUE) {
        curve(volMC(x), from = kmin, to = kmax, col = "brown", 
              lwd = 2, lty = 1, add = T)
      }
    }
    kIn <- k[!is.na(midVol)]
    volIn <- midVol[!is.na(midVol)]
    volInterp <- function(xout) {
      stinterp(x = kIn, y = volIn, xout)$y
    }
    atmVols[slice] <- volInterp(0)
    sig <- atmVols[slice] * sqrt(t)
    atmSkew[slice] <- (volInterp(sig/10) - volInterp(-sig/10))/(2 * 
                                                                  sig/10)
    atmCurv[slice] <- (volInterp(sig/10) + volInterp(-sig/10) - 
                         2 * atmVols[slice])/(2 * (sig/10)^2)
  }
  par(mfrow = c(1, 1), mex = 1)
  par(new = F)
  return(list(expiries = expDates, atmVols = atmVols, atmSkew = atmSkew, 
              atmCurv = atmCurv, atmVolsMC = atmVolsMC, atmSkewMC = atmSkewMC, 
              atmErrMC = atmErrMC, atmVolsMC2 = atmVolsMC2, atmSkewMC2 = atmSkewMC2, 
              atmErrMC2 = atmErrMC2))
}

