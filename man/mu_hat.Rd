% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPADE_estimate.R
\name{mu_hat}
\alias{mu_hat}
\title{The estimate of mu}
\usage{
mu_hat(delta, UTy, UT1, S)
}
\arguments{
\item{delta}{Delta value in the variance function.}

\item{UTy}{The expression UTy will need to be re-computed for each gene,}

\item{UT1}{The expression UT1 only depends on the coordinates X and can be pre-computed and reused for each gene.}

\item{S}{A diagonal matrix used in the spectral decomposition of the covariance.}
}
\value{
This function returns the estimate of mean expression level for each gene.
\item{mu_h}{The estimate of mean expression level for each gene.}
}
\description{
The estimate of mean expression level for each gene.
}
