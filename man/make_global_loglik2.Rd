% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convoSPAT_paramEst.R
\name{make_global_loglik2}
\alias{make_global_loglik2}
\title{Constructor functions for global parameter estimation.}
\usage{
make_global_loglik2(data, Xmat, Corr, obs.nuggets, nugg2.var)
}
\arguments{
\item{data}{A vector or matrix of data to use in the likelihood
calculation.}

\item{Xmat}{The design matrix for the mean model.}

\item{Corr}{The correlation matrix.}

\item{obs.nuggets}{A vector containing the spatially-varying nuggets
corresponding to each data location.}

\item{nugg2.var}{Fixed values for the covariance of the second nugget term.}
}
\value{
This function returns another function for use in \code{optim}.
}
\description{
This function generates another function to be used within \code{optim} to
obtain maximum likelihood estimates of
global variance parameter sigmasq with a fixed correlation
matrix (smoothness is fixed). The nugget variance is taken
to be spatially-varing.
}
\examples{
\dontrun{
make_global_loglik2( data, Xmat, Corr, obs.nuggets, nugg2.var )
}

}

