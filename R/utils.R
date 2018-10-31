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
#' \href{http://doi.org/10.2307/2532852}{Garthwaite (1996)}. Most of the
#' parameter names are consistent with the notation in
#' \href{http://doi.org/10.2307/2532852}{Garthwaite (1996)}. This is a general
#' utility function used within test and confidence interval functions.
#'
#' @param init initial (or most recent) estimate (e.g. L_i or U_i)
#' @param thetahat point estimate of theta (parameter of interest) based on
#' original data
#' @param t test statistic for new permutation
#' @param tstar test statistic for original permutation
#' @param alpha corresponds to 100(1 - alpha)\% CI
#' @param i iteration of the search process
#' @param m an optional initial magnitude of the steps; if left unspecified,
#' m defaults to recommended value proposed in Garthwaite and Buckland (1992)
#' @param bound "lower" or "upper"; i.e. which confidence bound you want
update_rm <- function(init, thetahat, t, tstar, alpha, i, m, bound) {
  if (!(bound %in% c("lower", "upper")))
    stop("please choose 'lower' or 'upper' for bound")

  g.alpha <- alpha / 2 # alpha as defined in Garthwaite paper
  # "m" defaults to suggestion in Garthwaite and Buckland 92
  if (missing(m))
    m <- min(c(50, ceiling(0.3 * (2 - g.alpha) / g.alpha)))

  # k for c length constant
  z <- qnorm(1 - g.alpha)
  k <- 2 * sqrt(2 * pi) * exp(z^2 / 2) / z

  if (bound == "lower") {
    c <- k * (thetahat - init)
    if (t < tstar) {
      out <- init + (c * g.alpha / (m + i - 1))
    } else {
      out <- init - (c * (1 - g.alpha) / (m + i - 1))
    }
  } else if (bound == "upper") {
    c <- k * (init - thetahat)
    if (t > tstar) {
      out <- init - (c * g.alpha / (m + i - 1))
    } else {
      out <- init + (c * (1 - g.alpha) / (m + i - 1))
    }
  }

  return(out)
}
