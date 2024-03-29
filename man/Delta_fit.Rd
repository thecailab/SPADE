% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPADE_DE.R
\name{Delta_fit}
\alias{Delta_fit}
\title{Estimation of the optimal delta}
\usage{
Delta_fit(location, y, L)
}
\arguments{
\item{location}{Spatial coordinates for all spots.}

\item{y}{Read counts data for the gene.}

\item{L}{The length-scale hyperparameter in the kernel function.}
}
\value{
This function returns the optimal log transformed delta value.
\item{delta}{Delta in the variance function.}
\item{Tao_hat}{Estimated optimal gamma value in the variance function.}
}
\description{
Estimation of the optimal delta in the variance function between groups.
}
