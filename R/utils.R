#' Permute randomized treatment assignments
#'
#' \code{permute} returns data with treatment assigments permuted
#'
#' This function inputs a data frame and outputs the same data frame with the
#' \code{trtname} column permuted according to the randomization scheme
#' provided. This is a general utility function used within test and confidence
#' interval functions.
#'
#' @param data a data frame
#' @param trtname character string specifying the name of randomized treatment
#' variable in \code{data} (variable to permute)
#' @param runit character string specifying the name of unit of randomization
#' in \code{data}
#' @param strat an optional character string specifying the name of the variable
#' in \code{data} upon which randomization was stratified
permute <- function(data, trtname, runit, strat = NULL) {
  if (is.null(strat)) {
    # permute
    design <- unique(data[, c(runit, trtname)])
    pdesign <- design
    pdesign[, trtname] <- sample(design[, trtname])
  } else {
    # stratified permute
    design <- unique(data[, c(runit, strat, trtname)])
    pdesign <- design
    for (s in unique(design[, strat])) {
      ptrt <- sample(design[design[, strat] == s, trtname])
      pdesign[design[, strat] == s, trtname] <- ptrt
    }
    pdesign <- pdesign[, - which(names(pdesign) == strat)]
  }

  data <- merge(data[, -which(names(data) == trtname)], pdesign, by = runit)
  return(data) # return data frame with permuted trtname
}

#' Perform one update using Robbins-Monro search process
#'
#' \code{update_rm} performs one update using Robbins-Monro search process
#'
#' This function inputs values necessary to perform one update using the
#' Robbins-Monro search procedure specific to confidence intervals proposed by
#' \href{http://doi.org/10.2307/2532852}{Garthwaite (1996)} (for method = 'G')
#' and
#' \href{https://doi.org/10.1198/jcgs.2009.0011}{Garthwaite and Jones (2009)}
#' (for method = 'GJ'). This is a general utility function called within
#' the confidence interval functions (e.g permci_glm()).
#'
#' @param method if method = 'G' (default), then search is carried out as
#' described in \href{http://doi.org/10.2307/2532852}{Garthwaite (1996)}. For
#' longer searches (nperm >= 200,000), method = 'GJ' is recommended and carried
#' out as outlined in
#' \href{https://doi.org/10.1198/jcgs.2009.0011}{Garthwaite and Jones (2009)}.
#' @param init initial (or most recent) estimate (e.g. L_i or U_i)
#' @param thetahat point estimate of theta (parameter of interest) based on
#' original data
#' @param t test statistic for new permutation
#' @param tstar test statistic for original permutation
#' @param alpha corresponds to 100(1 - alpha)\% CI
#' @param i iteration of the search process
#' @param m an optional initial magnitude of the steps; if left unspecified,
#' m defaults to recommended value proposed in Garthwaite and Buckland (1992)
#' @param k step length multiplier
#' @param Ps if method = 'GJ', vector of search lengths for each phase (if
#' unspecified, defaults to recommended values in
#' \href{https://doi.org/10.1198/jcgs.2009.0011}{Garthwaite and Jones (2009)})
#' @param bound "lower" or "upper"; i.e. which confidence bound you want
update_rm <- function(method, init, thetahat, t, tstar, alpha, i, m, k, Ps,
                      bound) {
  if (!(bound %in% c("lower", "upper")))
    stop("please choose 'lower' or 'upper' for bound")

  # step length constant
  if (bound == "lower") {
    c <- k * (thetahat - init)
  } else {
    c <- k * (init - thetahat)
  }

  # step size depends on method/phase/step
  if (method == 'G') {
    a <- c / (m + i)
  } else if (method == 'GJ') {
    if (i < Ps[1]) {
      # Phase 1
      a <- c / (m + i)
    } else if (i < Ps[1] + Ps[2]) {
      # Phase 2
      a <- c / (m + Ps[1])
    } else {
      # Phase 3
      a <- c / ((m + i) * (m + Ps[1]) / (m + Ps[1] + Ps[2]))
    }
  }

  if (bound == "lower") {
    if (t < tstar) {
      out <- init + a * (alpha / 2)
    } else {
      out <- init - a * (1 - (alpha / 2))
    }
  } else if (bound == "upper") {
    if (t > tstar) {
      out <- init - a * (alpha / 2)
    } else {
      out <- init + a * (1 - (alpha / 2))
    }
  }

  return(out)
}
