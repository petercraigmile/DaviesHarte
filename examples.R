
## ======================================================================
## Copyright 2002--2014, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================

## ======================================================================
## File     : examples.R
## Contains : R code to demonstrate the simulation of fractionally
##            differenced (FD), fractional Gaussian noise (fGn), and
##            autoregressive fractionally integrated moving average
##            (ARFIMA) processes.
## Version  : 0.3
## Updated  : pfc@stat.osu.edu, April 2014.
## ======================================================================


source("DaviesHarte.R")

source("FD_fGn_acvs.R")


## Simulate a Gaussian fractionally differenced process with differencing
## parameter d=0.25 of length 256.
## (this is a stationary long memory process)
## If d < 0 we obtain a antipersistent process, d=0 is a white noise
## process, and d > 0 a long memory process.
x <- Davies.Harte.sim(256, fd.acvs, d=0.25)


## Produce a time series plot and estimate the sample autocorrelation
par(mfrow=c(2,2), cex=0.75, mar=c(4,4,1,1), mgp=c(2,0.5,0), bty="L")
plot(ts(x), ylab="Simulated FD(0.25) process")
acf(x, main="", ylab="Sample ACF")
pacf(x, main="", ylab="Sample PACF")




## When d>0.5 we need to decompose d into fractional (stationary)
## and integer (nonstationary) parts.
## We can use the fd.fract.int.d function to do this.
## See the example below.

## Simulate a Gaussian fractionally differenced process with d=0.85 of length 256
## This is a nonstationary process, so we simulate assuming that X_0 = 0.

d.parts <- fd.fract.int.d(0.85)

y <- Davies.Harte.sim(256, fd.acvs, d=d.parts$fract, csum=d.parts$int)

## Produce a time series plot and estimate the sample autocorrelation.
## Note the nonstationary nature of the time series realization.
par(mfrow=c(2,2), cex=0.75, mar=c(4,4,1,1), mgp=c(2,0.5,0), bty="L")
plot(ts(y), ylab="Simulated FD(0.85) process")
acf(y, main="", ylab="Sample ACF")
pacf(y, main="", ylab="Sample PACF")





## Simulate fractional Gaussian noise (fGn) with Hurst parameter H=0.75
## of length 200.
## (this is a stationary long memory process)
z <- Davies.Harte.sim(200, fGn.acvs, H=0.75)

## Produce a time series plot and estimate the sample autocorrelation.
par(mfrow=c(2,2), cex=0.75, mar=c(4,4,1,1), mgp=c(2,0.5,0), bty="L")
plot(ts(z), ylab="Simulated fGn process")
acf(z, main="", ylab="Sample ACF")
pacf(z, main="", ylab="Sample PACF")





## To simulate a Gaussian ARFIMA(p,d,q) process, we feed the output from
## simulating an FD(d) process to the arima.sim function.

## Here we simulate an ARFIMA(1,d,0) process of length 500
## with AR parameter 0.6, d=0.25, and sigma^2=2.

w <- arima.sim(model=list(ar=0.6), n=500,
               innov=Davies.Harte.sim(500, fd.acvs, d=0.25, sigma2=2))

## Produce a time series plot and estimate the sample autocorrelation.
par(mfrow=c(2,2), cex=0.75, mar=c(4,4,1,1), mgp=c(2,0.5,0), bty="L")
plot(ts(w), ylab="Simulated ARFIMA process")
acf(w, main="", ylab="Sample ACF")
pacf(w, main="", ylab="Sample PACF")
