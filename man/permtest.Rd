% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{permtest}
\alias{permtest}
\title{Randomization test}
\usage{
permtest(f, trtname, runit, strat = NULL, data, nperm = 1000,
  ncores = 1, seed, quietly = T)
}
\arguments{
\item{f}{fitted model object or function. If \code{f} is a fitted model
object, then the coefficient corresponding to \code{trtname} is used as the
test statistic. If \code{f} is a function, it must be defined such that when
applied to \code{data}, it returns the observed univariate test statistic.
This function can be as simple or complex as desired as long as its one input
argument is a data frame structured the same as \code{data}.}

\item{trtname}{character string specifying the name of randomized treatment
variable in \code{data} (variable to permute)}

\item{runit}{character string specifying the name of unit of randomization
in \code{data}}

\item{strat}{an optional character string specifying the name of the variable
in \code{data} upon which randomization was stratified}

\item{data}{a data frame containing the variables necessary for the function
\code{f}. This argument is passed to \code{f}.}

\item{nperm}{number of permutations for randomization test}

\item{ncores}{number of cores to use for computation. If ncores > 1, permtest
runs in parallel.}

\item{seed}{a numerical seed to use, passed to \code{\link[base]{set.seed}}
(if \code{ncores == 1}) or \code{\link[doRNG]{registerDoRNG}} (if
\code{ncores > 1}).}

\item{quietly}{logical; if TRUE (and if ncores == 1), status updates will be
printed to Console otherwise, suppress updates.}
}
\description{
Carry out a randomization test for a treatment effect using a fitted model
object or user defined test statistic.
}
\seealso{
\code{\link[permuter]{permci}} for a randomization-based CI
}
