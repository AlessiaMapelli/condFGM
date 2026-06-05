# Helper functions to compute basis function values at given time points.
library(fda)

# @param t      numeric vector of time points in [0, 1]
# @param M      integer; total number of Fourier basis functions
# @param l      integer; index of the basis function to evaluate
# @param omega  numeric; frequency (default 1)
# @return       numeric vector of length(t) with basis function values
fourier.basis <- function(t, M, l, omega=1){
  n <- length(t)
  b <- numeric(n)
  for (i in 1:n){ # basis function defined on [0,1]
    if (t[i]>=0 & t[i]<=1){
      if(l == 1) b[i] <- 1
      else if(l %% 2 == 0) b[i] <- sqrt(2) * cos(2 * pi * omega * (l/2) * t[i])
      else b[i] <- sqrt(2) * sin(2 * pi * omega * ((l-1)/2) * t[i])
    } else {
      b[i] <- 0
    }
  }
  return(b)
}

# @param obs.time  numeric vector; observation grid of length tau
# @param M         integer; number of Fourier basis functions
# @param omega     numeric; frequency (default 1)
# @return          tau x M matrix of observed basis function values
obs.fourier.bases <- function(obs.time, M, omega=1){
  tau <- length(obs.time)
  ObservMat <- matrix(0, nrow=n, ncol=M)
  for (k in 1:tau){
    for (l in 1:M){
      ObservMat[k, l] <- fourier.basis(obs.time[k], M, l, omega)
    }
  }
  return(ObservMat)
}

# @param obs.time  numeric vector; observation grid
# @param M         integer; number of Fourier basis functions
# @param period    numeric; period of the basis (must be <= 1; default 1)
# @return          tau x M matrix of Fourier basis values evaluated at obs.time
fda.fourier.mat <- function(obs.time, M, period=1){
  basis.fourier <- create.fourier.basis(rangeval = c(0,period), nbasis=M)
  ObservMat <- eval.basis(evalarg=obs.time, basisobj=basis.fourier)
  return(ObservMat)
}

# @param obs.time  numeric vector in [0, 1]; observation grid
# @param M         integer; number of B-spline basis functions
# @return          tau x M matrix of B-spline basis values evaluated at obs.time
fda.bspline.mat <- function(obs.time, M){
  basis.bspline <- create.bspline.basis(rangeval = c(0,1), nbasis=M)
  ObservMat <- eval.basis(evalarg=obs.time, basisobj=basis.bspline)
  return(ObservMat)
}

# @param obs.time  numeric vector; observation grid
# @param M         integer; number of exponential basis functions
# @return          tau x M matrix of exponential basis values evaluated at obs.time
fda.exp.mat <- function(obs.time, M){
  basis.exp <- create.exponential.basis(nbasis=M)
  ObservMat <- eval.basis(evalarg=obs.time, basisobj=basis.exp)
  return(ObservMat)
}

# @param obs.time  numeric vector; observation grid
# @param M         integer; number of monomial basis functions
# @return          tau x M matrix of monomial basis values evaluated at obs.time
fda.monomial.mat <- function(obs.time, M){
  basis.monomial <- create.monomial.basis(nbasis = M)
  ObservMat <- eval.basis(evalarg=obs.time, basisobj = basis.monomial)
  return(ObservMat)
}

# @param obs.time  numeric vector; observation grid
# @param M         integer; number of power basis functions
# @return          tau x M matrix of power basis values evaluated at obs.time
fda.power.mat <- function(obs.time, M){
  basis.power <- create.power.basis(nbasis = M)
  ObservMat <- eval.basis(evalarg=obs.time, basisobj = basis.power)
  return(ObservMat)
}
