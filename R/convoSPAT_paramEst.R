#======================================================================================
# Local likelihood estimation for covariance functions with spatially-varying
# parameters: the convoSPAT() package for R
# Mark D. Risser / The Ohio State University / 2014-2015
#======================================================================================

#======================================================================================
# Local/Global Parameter Estimation
#======================================================================================

#======================================================================================
# Global parameter estimation
#======================================================================================
# For a given correlation matrix, the following functions calculate the global
# parameters based on what is NOT specified to be spatially-varying.
#======================================================================================

#===================================
# First, models without kappa
# Estimates: nugget, variance

#ROxygen comments ----
#' Global parameter estimation for the nonstationary model.
#'
#' This function calculates the maximum likelihood estimates of
#' global variance parameters tausq, sigmasq with a fixed correlation
#' matrix (smoothness is fixed).
#'
#' @param params A vector of inputs: tausq and sigmasq.
#'
#' @return This function returns the loglikelihood for a particular
#' combination of (tausq, sigmasq).
#'
#' @examples
#' \dontrun{
#' overall_lik1( c(0.5,1) )
#' }
#'
#' @export

overall_lik1 <- function(params){

  # Parameters
  tausq <- params[1]
  sigmasq <- params[2]

  N <- dim(fit.data)[1]
  p <- dim(fit.data)[2]

  Cov.inv <- (1/tausq)*( diag(rep(1,N)) - (sigmasq/tausq)*Vmat
                         %*% diag( 1/(1/diag(Dmat) + rep(sigmasq/tausq,N)) )%*%t(Vmat) )

  XCX <- t(Xmat) %*% Cov.inv %*% Xmat
  XCX.inv <- chol2inv( chol(XCX) )

  Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

  # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
  log.det.Cov <- sum( log( 1/diag(Dmat) + rep(sigmasq/tausq,N) ) ) + sum( log(diag(Dmat)) ) + N*log(tausq)
  loglikelihood <- 0.5*( p*log.det.Cov + p*log(det(XCX)) + sum(diag(t(fit.data)%*%Ptemp%*%fit.data)) )

  if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
  return(loglikelihood)

}

# Estimates: variance

#ROxygen comments ----
#' Global parameter estimation for the nonstationary model.
#'
#' This function calculates the maximum likelihood estimates of
#' global variance parameter sigmasq with a fixed correlation
#' matrix (smoothness is fixed). The nugget variance is taken
#' to be spatially-varing.
#'
#' @param params Scalar input for sigmasq.
#'
#' @return This function returns the loglikelihood for a particular
#' value of sigmasq.
#'
#' @examples
#' \dontrun{
#' overall_lik2( 1 )
#' }
#'
#' @export

overall_lik2 <- function(params){

  # Parameters
  sigmasq <- params[1]

  N <- dim(fit.data)[1]
  p <- dim(fit.data)[2]

  Edecomp <- eigen(sigmasq*Corr + diag(fit.nuggets))
  Dmat.temp <- diag(Edecomp$values)
  Vmat.temp <- Edecomp$vectors

  Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

  XCX <- t(Xmat) %*% Cov.inv %*% Xmat
  XCX.inv <- chol2inv( chol(XCX) )

  Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

  # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
  loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(fit.data)%*%Ptemp%*%fit.data)) )

  if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
  return(loglikelihood)

}


# Estimates: nugget

#ROxygen comments ----
#' Global parameter estimation for the nonstationary model.
#'
#' This function calculates the maximum likelihood estimates of
#' global variance parameter tausq with a fixed correlation
#' matrix (smoothness is fixed). The process variance is taken
#' to be spatially-varing.
#'
#' @param params Scalar input for sigmasq.
#'
#' @return This function returns the loglikelihood for a particular
#' value of sigmasq.
#'
#' @examples
#' \dontrun{
#' overall_lik2( 1 )
#' }
#'
#' @export

overall_lik3 <- function(params){

  # Parameters
  tausq <- params[1]

  N <- dim(fit.data)[1]
  p <- dim(fit.data)[2]

  Edecomp <- eigen(diag(sqrt(fit.variance)) %*% Corr %*% diag(sqrt(fit.variance)) + diag(rep(tausq,N)))
  Dmat.temp <- diag(Edecomp$values)
  Vmat.temp <- Edecomp$vectors

  Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

  XCX <- t(Xmat) %*% Cov.inv %*% Xmat
  XCX.inv <- chol2inv( chol(XCX) )

  Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

  # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
  loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(fit.data)%*%Ptemp%*%fit.data)) )

  if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
  return(loglikelihood)

}

#===================================

#===================================
# Next, models with kappa
# Estimates: nugget, variance, kappa

#ROxygen comments ----
#' Global parameter estimation for the nonstationary model.
#'
#' This function calculates the maximum likelihood estimates of
#' global variance parameters tausq, sigmasq, and nu.
#'
#' @param params A vector of inputs: tausq, sigmasq, and nu.
#'
#' @return This function returns the loglikelihood for a particular
#' combination of (tausq, sigmasq, nu).
#'
#' @examples
#' \dontrun{
#' overall_lik1_kappa( c(0.5,1,1) )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

overall_lik1_kappa <- function(params){

  # Parameters
  tausq <- params[1]
  sigmasq <- params[2]
  kapp <- params[3]

  Unscl.corr <- cov.spatial( Distmat, cov.model = covmodel,
                             cov.pars = c(1,1), kappa = kapp )
  Corr <- Scalemat*Unscl.corr

  N <- dim(fit.data)[1]
  p <- dim(fit.data)[2]

  Edecomp <- eigen(sigmasq*Corr + diag(rep(tausq,N)))
  Dmat.temp <- diag(Edecomp$values)
  Vmat.temp <- Edecomp$vectors

  Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

  XCX <- t(Xmat) %*% Cov.inv %*% Xmat
  XCX.inv <- chol2inv( chol(XCX) )

  Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

  # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
  loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(fit.data)%*%Ptemp%*%fit.data)) )

  if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
  return(loglikelihood)

}


# Estimates: variance, kappa
#ROxygen comments ----
#' Global parameter estimation for the nonstationary model.
#'
#' This function calculates the maximum likelihood estimates of
#' global variance parameters sigmasq and nu. The nugget variance is
#' taken to be spatially-varying.
#'
#' @param params A vector of inputs: sigmasq and nu.
#'
#' @return This function returns the loglikelihood for a particular
#' combination of (sigmasq, nu).
#'
#' @examples
#' \dontrun{
#' overall_lik2_kappa( c(1,1) )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

overall_lik2_kappa <- function(params){

  # Parameters
  sigmasq <- params[1]
  kapp <- params[2]

  Unscl.corr <- cov.spatial( Distmat, cov.model = covmodel,
                             cov.pars = c(1,1), kappa = kapp )
  Corr <- Scalemat*Unscl.corr

  N <- dim(fit.data)[1]
  p <- dim(fit.data)[2]

  Edecomp <- eigen(sigmasq * Corr + diag(fit.nuggets))
  Dmat.temp <- diag(Edecomp$values)
  Vmat.temp <- Edecomp$vectors

  Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

  XCX <- t(Xmat) %*% Cov.inv %*% Xmat
  XCX.inv <- chol2inv( chol(XCX) )

  Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

  # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
  loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(fit.data)%*%Ptemp%*%fit.data)) )

  if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
  return(loglikelihood)

}

# Estimates: nugget, kappa
#ROxygen comments ----
#' Global parameter estimation for the nonstationary model.
#'
#' This function calculates the maximum likelihood estimates of
#' global variance parameters tausq and nu. The process variance is
#' taken to be spatially-varying.
#'
#' @param params A vector of inputs: tausq and nu.
#'
#' @return This function returns the loglikelihood for a particular
#' combination of (tausq, nu).
#'
#' @examples
#' \dontrun{
#' overall_lik3_kappa( c(0.5,1) )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

overall_lik3_kappa <- function(params){

  # Parameters
  tausq <- params[1]
  kapp <- params[2]

  Unscl.corr <- cov.spatial( Distmat, cov.model = covmodel,
                             cov.pars = c(1,1), kappa = kapp )
  Corr <- Scalemat*Unscl.corr

  N <- dim(fit.data)[1]
  p <- dim(fit.data)[2]

  Edecomp <- eigen(diag( sqrt(fit.variance) ) %*% Corr %*% diag( sqrt(fit.variance) ) + tausq*diag(rep(1,N)))
  Dmat.temp <- diag(Edecomp$values)
  Vmat.temp <- Edecomp$vectors

  Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

  XCX <- t(Xmat) %*% Cov.inv %*% Xmat
  XCX.inv <- chol2inv( chol(XCX) )

  Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

  # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
  loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(fit.data)%*%Ptemp%*%fit.data)) )

  if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
  return(loglikelihood)

}

# Estimates: kappa
#ROxygen comments ----
#' Global parameter estimation for the nonstationary model.
#'
#' This function calculates the maximum likelihood estimates of
#' global variance parameters nu. The process variance
#' and nugget variance are taken to be spatially-varying.
#'
#' @param params A vector of inputs: nu.
#'
#' @return This function returns the loglikelihood for a particular
#' combination of (nu).
#'
#' @examples
#' \dontrun{
#' overall_lik4_kappa( c(1) )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

overall_lik4_kappa <- function(params){

  # Parameters
  kapp <- params[1]

  Unscl.corr <- cov.spatial( Distmat, cov.model = covmodel,
                             cov.pars = c(1,1), kappa = kapp )
  Corr <- Scalemat*Unscl.corr

  N <- dim(fit.data)[1]
  p <- dim(fit.data)[2]

  Edecomp <- eigen(diag( sqrt(fit.variance) ) %*% Corr %*% diag( sqrt(fit.variance) ) + diag(fit.nuggets))
  Dmat.temp <- diag(Edecomp$values)
  Vmat.temp <- Edecomp$vectors

  Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

  XCX <- t(Xmat) %*% Cov.inv %*% Xmat
  XCX.inv <- chol2inv( chol(XCX) )

  Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

  # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
  loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(fit.data)%*%Ptemp%*%fit.data)) )

  if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
  return(loglikelihood)

}
#======================================================================================

#======================================================================================
# Function to calculate the locally anisotropic model
#======================================================================================
# Using a subset of the data, calculate MLEs of covariance and mean
# parameters. Options for estimating models with and without kappa.
#======================================================================================

# Likelihood for Gaussian data with anisotropic covariance,
# without kappa:

#ROxygen comments ----
#' Local parameter estimation for the nonstationary model.
#'
#' This function calculates the maximum likelihood estimates of
#' local variance parameters lam1, lam2, eta, tausq, and sigmasq,
#' assuming the smoothness is fixed, using a Gaussian likelihood with
#' an anisotropic covariance structure.
#'
#' @param params A vector of inputs: lam1, lam2, eta, tausq, and sigmasq
#'
#' @return This function returns the loglikelihood for a particular
#' combination of (lam1, lam2, eta, tausq, sigmasq).
#'
#' @examples
#' \dontrun{
#' anis_model( c(1,1,0,0.1,1) )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom StatMatch mahalanobis.dist

anis_model <- function(params){

  # Kernel parameters
  lam1 <- params[1]
  lam2 <- params[2]
  eta <- params[3]

  # Nugget
  tausq <- params[4]

  # Process variance
  sigmasq <- params[5]

  #================================
  N <- dim(temp.locations)[1]
  p <- dim(temp.data)[2]

  Pmat <- matrix(c(cos(eta),-sin(eta),sin(eta),cos(eta)),nrow=2,byrow=T)
  Lmat <- diag(c(lam1,lam2))

  Sigma <- Pmat %*% Lmat %*% t(Pmat)

  distances <- mahalanobis.dist( data.x = temp.locations, vc = Sigma )
  NS.cov <- sigmasq*cov.spatial(distances, cov.model = covmodel,
                         cov.pars = c(1,1), kappa = KAPPA)

  Edecomp <- eigen(NS.cov + diag(rep(tausq,N)))
  Dmat.temp <- diag(Edecomp$values)
  Vmat.temp <- Edecomp$vectors

  Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

  XCX <- t(Xtemp) %*% Cov.inv %*% Xtemp
  XCX.inv <- chol2inv( chol(XCX) )

  Ptemp <- Cov.inv - Cov.inv %*% Xtemp %*% XCX.inv %*% t(Xtemp) %*% Cov.inv

  # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
  loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(temp.data)%*%Ptemp%*%temp.data)) )

  if(abs(loglikelihood) == Inf){ loglikelihood <- 1000000 }

  return(loglikelihood)

}

# with kappa:
#ROxygen comments ----
#' Local parameter estimation for the nonstationary model.
#'
#' This function calculates the maximum likelihood estimates of
#' local variance parameters lam1, lam2, eta, tausq, sigmasq, and
#' nu (smoothness) using a Gaussian likelihood with
#' an anisotropic covariance structure.
#'
#' @param params A vector of inputs: lam1, lam2, eta, tausq, sigmasq, and nu
#'
#' @return This function returns the loglikelihood for a particular
#' combination of (lam1, lam2, eta, tausq, sigmasq, nu).
#'
#' @examples
#' \dontrun{
#' anis_model_kappa( c(1,1,0,0.1,1,1) )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom StatMatch mahalanobis.dist

anis_model_kappa <- function(params){

  # Kernel parameters
  lam1 <- params[1]
  lam2 <- params[2]
  eta <- params[3]

  # Nugget
  tausq <- params[4]

  # Process variance
  sigmasq <- params[5]

  # Smoothness
  KAPPA <- params[6]

  #================================
  N <- dim(temp.locations)[1]
  p <- dim(temp.data)[2]

  Pmat <- matrix(c(cos(eta),-sin(eta),sin(eta),cos(eta)),nrow=2,byrow=T)
  Dmat <- diag(c(lam1,lam2))

  Sigma <- Pmat %*% Dmat %*% t(Pmat)

  distances <- mahalanobis.dist( data.x = temp.locations, vc = Sigma )
  NS.cov <- sigmasq*cov.spatial(distances, cov.model = covmodel,
                                cov.pars = c(1,1), kappa = KAPPA)

  Edecomp <- eigen(NS.cov + diag(rep(tausq,N)))
  Dmat.temp <- diag(Edecomp$values)
  Vmat.temp <- Edecomp$vectors

  Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

  XCX <- t(Xtemp) %*% Cov.inv %*% Xtemp
  XCX.inv <- chol2inv( chol(XCX) )

  Ptemp <- Cov.inv - Cov.inv %*% Xtemp %*% XCX.inv %*% t(Xtemp) %*% Cov.inv

  # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
  loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(temp.data)%*%Ptemp%*%temp.data)) )

  if(abs(loglikelihood) == Inf){ loglikelihood <- 1000000 }

  return(loglikelihood)

}

