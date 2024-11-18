library(stinepack);


plotIvols <- function (ivolData, sviMatrix = NULL, slices = NULL, rowcol=NULL, modelVol = FALSE, 
                       plot = T) 
{
  bidVols <- as.numeric(ivolData$Bid)
  include <- !is.na(bidVols)
  ivolData <- ivolData[include, ]
  bidVols <- as.numeric(ivolData$Bid)
  askVols <- as.numeric(ivolData$Ask)
  expDates <- unique(ivolData$Texp)
  k <- log(ivolData$Strike/ivolData$Fwd)
  texp <- ivolData$Texp
  midVols <- (bidVols + askVols)/2
  if (modelVol) {
    modelVols <- as.numeric(ivolData$modelVol)
  }
  nSlices <- length(expDates)
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
  atmVol <- numeric(nSlices)
  atmSkew <- numeric(nSlices)
  atmCurv <- numeric(nSlices)
  par(mfrow = c(rows, columns), mex = 0.5)
  for (slice in slices) {
    t <- expDates[slice]
    bidVols.t <- bidVols[texp == t]
    askVols.t <- askVols[texp == t]
    midVols.t <- midVols[texp == t]
    if (modelVol) {
      modelVols.t <- modelVols[texp == t]
    }
    k.t <- k[texp == t]
    kmin <- min(k.t)
    kmax <- max(k.t)
    ybottom <- 0.8 * min(bidVols.t, na.rm = T)
    ytop <- 1.2 * max(askVols.t, na.rm = T)
    xrange <- c(kmin, kmax)
    yrange <- c(ybottom, ytop)
    if (plot) {
      plot(k.t, bidVols.t, col = "red2", pch = 24, cex = 0.5, 
           xlim = xrange, ylim = yrange, main = paste("T =", 
                                                      format(t, digits = 2, nsmall = 2)), xlab = "Log-Strike", 
           ylab = "Implied Vol.")
      points(k.t, askVols.t, col = "royalblue", pch = 25, cex = 0.5)
      if (modelVol) {
        lines(k.t, modelVols.t, col = "green4", lwd = 2)
      }
      if ((!is.null(sviMatrix))) {
        vol <- function(k) {
          sqrt(svi(sviMatrix[slice, ], k)/t)
        }
        curve(vol(x), from = kmin, to = kmax, col = "brown", 
              lwd = 2, add = T)
      }
    }
    volInterp <- function(xout) {
      stinterp(x = k.t, y = midVols.t, xout)$y
    }
    atmVol[slice] <- volInterp(0)
    sig <- atmVol[slice] * sqrt(t)
    atmSkew[slice] <- (volInterp(sig/10) - volInterp(-sig/10))/(2 * 
                                                                  sig/10)
    atmCurv[slice] <- (volInterp(sig/10) + volInterp(-sig/10) - 
                         2 * atmVol[slice])/(2 * (sig/10)^2)
  }
  par(mfrow = c(1, 1), mex = 1)
  par(new = F)
  return(list(expiries = expDates, atmVol = atmVol, atmSkew = atmSkew, 
              atmCurv = atmCurv))
}