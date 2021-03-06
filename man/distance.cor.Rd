% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures.R
\name{distance.cor}
\alias{distance.cor}
\title{Measure an impact of the covariates on the response using the distance correlation}
\usage{
distance.cor(x, y, subsamples, index = 1, ...)
}
\arguments{
\item{x}{Matrix with \code{n} observations of \code{p} covariates in each row.}

\item{y}{Response vector with \code{n} observations.}

\item{subsamples}{Matrix with \code{m} indices of \code{N} subsamples in each column.}

\item{index}{Positive scalar.}

\item{...}{Not in use.}
}
\value{
Numeric \code{p} by \code{N} matrix with distance correlations evaluated for each subsample.
}
\description{
This function evaluates the distance correlation between the response \code{y} and each column in the design matrix \code{x} over subsamples in \code{subsamples}.
}
\references{
Maria L. Rizzo and Gabor J. Szekely (2014). energy: E-statistics
(energy statistics). R package version 1.6.1 (\url{http://CRAN.R-project.org/package=energy}).
}
