% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{update_rm}
\alias{update_rm}
\title{Perform one update using Robbins-Monro search process}
\usage{
update_rm(method, init, thetahat, t, tstar, alpha, i, m, k, Ps, bound)
}
\arguments{
\item{method}{if method = 'G' (default), then search is carried out as
described in \href{http://doi.org/10.2307/2532852}{Garthwaite (1996)}. For
longer searches (nperm >= 200,000), method = 'GJ' is recommended and carried
out as outlined in
\href{https://doi.org/10.1198/jcgs.2009.0011}{Garthwaite and Jones (2009)}.}

\item{init}{initial (or most recent) estimate (e.g. L_i or U_i)}

\item{thetahat}{point estimate of theta (parameter of interest) based on
original data}

\item{t}{test statistic for new permutation}

\item{tstar}{test statistic for original permutation}

\item{alpha}{corresponds to 100(1 - alpha)\% CI}

\item{i}{iteration of the search process}

\item{m}{an optional initial magnitude of the steps; if left unspecified,
m defaults to recommended value proposed in Garthwaite and Buckland (1992)}

\item{k}{step length multiplier}

\item{Ps}{if method = 'GJ', vector of search lengths for each phase (if
unspecified, defaults to recommended values in
\href{https://doi.org/10.1198/jcgs.2009.0011}{Garthwaite and Jones (2009)})}

\item{bound}{"lower" or "upper"; i.e. which confidence bound you want}
}
\description{
\code{update_rm} performs one update using Robbins-Monro search process
}
\details{
This function inputs values necessary to perform one update using the
Robbins-Monro search procedure specific to confidence intervals proposed by
\href{http://doi.org/10.2307/2532852}{Garthwaite (1996)} (for method = 'G')
and
\href{https://doi.org/10.1198/jcgs.2009.0011}{Garthwaite and Jones (2009)}
(for method = 'GJ'). This is a general utility function called within
the confidence interval functions (e.g permci_glm()).
}
