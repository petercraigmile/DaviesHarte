## ======================================================================
## Copyright 2002--2014, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================

## ======================================================================
## File     : FD_fGn_acvs.R
## Contains : R code to calculate the autocovariance sequence (ACVS)
##            for a fractionally differenced (FD) and fractional
##            Gaussian noise (fGn).
## Version  : 0.3
## Updated  : pfc@stat.osu.edu, April 2014.
##
## References :
##
## 1. J. Beran, Statistics for Long Memory Processes. New York:
##    Chapman and Hall, 1994.
##
## 2. P. F. Craigmile (2003). Simulating a class of stationary
##    Gaussian processes using the Davies-Harte algorithm, with
##    application to long memory processes. Journal of Time Series
##    Analysis, 24, 505-511. (http://dx.doi.org/10.1111/1467-9892.00318)
## ======================================================================




## ======================================================================
## Purpose : Calculate the autocovariance sequence for the non-negative
##           lags of a FD('d','sigma2') process, up to a maximum lag of
##           'lag.max'.
## Assumes : lag.max >= 0
## Created : pfc@stat.osu.edu, Feb 2003.
## ======================================================================

fd.acvs <- function (lag.max, d, sigma2=1)
{
  acvs0 <- sigma2 * exp(lgamma(1-2*d)-2*lgamma(1-d))
  if (lag.max>0)
  {
    ks <- 1:lag.max
    cumprod(c(acvs0, (ks-1+d)/(ks-d)))
  }
  else acvs0
}




## ======================================================================
## Purpose : Calculate the autocovariance sequence for the non-negative
##           lags of a fGn('H','sigma2') process, up to a maximum lag of
##           'lag.max'.
## Assumes : lag.max >= 0
## ======================================================================

fGn.acvs <- function (lag.max, H, sigma2=1) {

  lags <- 0:lag.max
  two.H <- 2 * H
  
  0.5 * sigma2 * (abs(lags+1)^two.H - 2*abs(lags)^two.H + abs(lags-1)^two.H)
}




## ======================================================================
## Purpose : Calculates the integer and fractional differencing
##           components for a difference parameter 'd'.
## Returns : a list of 'd', the fractional part 'fract', and the
##           integer part 'int'.
## Created : pfc@stat.osu.edu, Feb 2003. 
## ======================================================================

fd.fract.int.d <- function (d)
{
  int   <- floor(d+0.5)
  fract <- d-int
  list(d=d, int=int, fract=fract)
}

