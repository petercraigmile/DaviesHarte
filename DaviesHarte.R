## ======================================================================
## Copyright 2002--2010, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================

## ======================================================================
## File     : DaviesHarte.R
## Contains : R code to simulate stationary Gaussian processes
##            (that satisfy a certain positivity condition) using the
##            Davies-Harte algorithm.
## Version  : 0.3
## Updated  : pfc@stat.osu.edu, March 2007.
##
## References:
##
## 1. R. B. Davies and and Harte. (1987) , "Tests for Hurst Effect",
##    Biometrika, 74, 96-101.
## 
## 2. A.T.A. Wood and G. Chan (1994), "Simulation of Stationary
##    Gaussian Processes in $[0,1]^d$", Journal of Computational and
##    Graphical Statistics, 3, 409-432.
## 
## 3. T. Gneiting, T. (2000) "Power-law correlations, related models
##    for long-range dependence, and their simulation". Journal of
##    Applied Probability, 37, 1104-1109.
## 
## 4. P. F. Craigmile (2003). Simulating a class of stationary
##    Gaussian processes using the Davies-Harte algorithm, with
##    application to long memory processes. Journal of Time Series
##    Analysis, 24, 505-511. (http://dx.doi.org/10.1111/1467-9892.00318)
## ======================================================================



## =======================================================================
## Purpose  : Sets up variables for quick simulation of stationary
##            Gaussian process using the Davie Harte algorithm.
##            'N' denotes the required length of the time series.
##            The function acvs.fun(max.lag, ...) calculates the
##            autocovariance sequence at lags 0:'max.lag'.
## Assumes  : 'N' > 0;
##            Process satisfies a certain positivity condition
##            (see Craigmile, 2003).
## =======================================================================

Davies.Harte.sim.setup <- function (N, acvs.fun, ...)
{
  # 'Npower2' is the next power of two greater than or equal to 'N'.
  Npower2 <- 2^ceiling(log2(N))

  # Calculate the autocovariance sequence
  acvs <- acvs.fun(Npower2, ...)

  # Calculate 'Sk' for use in the Davies-Harte algorithm.
  Sk <- Re(fft(c(acvs, acvs[Npower2:2]), inverse=F))

  ## Check the positivity condition.
  ##  We do not need to actually check this for certain processes!
  ##  (see the references above).
  ##  I leave this in for generality.
  if (any(Sk<=0)) stop("Some values in Davies.Harte.sim.setup are <= 0!!")  
  
  ## return a list of the useful values for use in the Davies-Harte
  ## simulation function.
  list(N=N, Npower2=Npower2, M=2*Npower2, sqrt.Sk=sqrt(Sk), ks=2:Npower2)
}





## =======================================================================
## Purpose : Simulate a stationary Gaussian process of length 'N' with
##           acvs defined in the function 'acvs.fun', using the
##           Davies-Harte algorithm.
##           The function acvs.fun(max.lag, ...) calculates the
##           autocovariance sequence at lags 0:'max.lag'.
##           The resulting process is cumulatively summed 'csum' times.
## Returns : A vector of the simulated process.
## Assumes : 'N' > 0;
##           Process satisfies a certain positivity condition
##           (see Craigmile, 2003).
## Note    : 'DH.obj' can be used to store the useful information from
##           'Davies.Harte.sim.setup' to save calculation time on
##           repeated simulations.
## ======================================================================== 

Davies.Harte.sim <- function (N, acvs.fun,
                              DH.obj=Davies.Harte.sim.setup(N, acvs.fun, ...),
                              csum=0, ...)
{
  ## simulate the random normals
  zs <- rnorm(DH.obj$M)

  ## calculate 'xs' and 'ys'.
  ## 'xs' is used in the calculation for 'ys'.
  ## 'ys' is the circular process defined in the Fourier domain.
  xs <- sqrt(0.5) * DH.obj$sqrt.Sk[DH.obj$ks] *
        complex(real=zs[2*DH.obj$ks-2], imag=zs[2*DH.obj$ks-1])
  ys <- c(DH.obj$sqrt.Sk[1] * zs[1], xs,
          DH.obj$sqrt.Sk[DH.obj$Npower2+1] * zs[DH.obj$M], Conj(rev(xs)))

  ## Generate the data 'x' by taking the inverse FFT of the 'ys'.
  ## We need to take the real part, and scale appropriatly.
  ## Only return the first 'N' values.
  x  <- Re(fft(ys, inverse=T))[1:DH.obj$N]/sqrt(DH.obj$M)

  ## Cumulatively sum the series 'csum' times.
  if (csum>0)
    for (k in 1:csum) x <- cumsum(x)

  ## return the simulated time series.
  x
}

