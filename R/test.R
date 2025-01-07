#' Randomization test for regression
#'
#' Carry out a randomization test for a treatment effect using a regression
#' parameter estimate as the test statistic.
#'
#' These functions are used to carry out randomization tests for a treatment
#' effect based on a parameter estimate from a regression model.
#' Different functions corrrespond to different regression models:
#' \itemize{
#'   \item \code{permtest_glm}: randomization test based on
#'   \code{\link[stats]{glm}}
#'   \item \code{permtest_survreg}: randomization test based on
#'   \code{\link[survival]{survreg}}
#'   \item \code{permtest_coxph}: randomization test based on
#'   \code{\link[survival]{coxph}}
#' }
#' To ensure correct specification of the parameters passed to the models above
#' (e.g. \code{formula} in \code{\link[survival]{survreg}}), please refer to
#' their documentation.
#'
#' @seealso \code{\link[permuter]{permci_glm}},
#' \code{\link[permuter]{permci_survreg}}, \code{\link[permuter]{permci_coxph}}
#' for corresponding randomization-based CIs
#'
#' @param formula an object of class "\code{\link[stats]{formula}}"
#' (or one that can be coerced to that class): a symbolic description of the
#' model to be fitted. This argument is passed to the corresponding regression
#' function, e.g. \code{\link[stats]{glm}} (see Details).
#' @inheritParams permute
#' @param family a description of the error distribution and link function to
#'     be used in the model. This can be a character string naming a family
#'     function, a family function or the result of a call to a family function.
#'     (See \code{\link[stats]{family}} for details.)
#' @param data a data frame containing the variables in the model. This argument
#' is passed to the corresponding regression function, e.g.
#' \code{\link[stats]{glm}} (see Details).
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two-sided" (default), "greater", or "less".
#' @param theta a number indicating the null treatment effect value of interest
#' for the randomization test
#' @param nperm number of permutations for randomization test
#' @param ncores number of cores to use for computation. If ncores > 1, permtest
#' runs in parallel.
#' @param seed a numerical seed to use, passed to \code{\link[base]{set.seed}}
#' (if \code{ncores == 1}) or \code{\link[doRNG]{registerDoRNG}} (if
#' \code{ncores > 1}).
#' @param quietly logical; if TRUE (and if ncores == 1), status updates will be
#' printed to Console otherwise, suppress updates.
#' @param dist assumed distribution for y variable. If the argument is a
#' character string, then it is assumed to name an element from
#' \code{\link[survival]{survreg.distributions}}. These include "weibull",
#' "exponential", "gaussian", "logistic","lognormal" and "loglogistic".
#' Otherwise, it is assumed to be a user defined list conforming to the format
#' described in \code{\link[survival]{survreg.distributions}}.
#' @import stats
#' @importFrom foreach %dopar%
#' @importFrom doRNG %dorng%
#' @importFrom survival Surv
#'
#' @examples
#' # Carry out a randomization test to determine whether there was a difference
#' # in the rate of bacterial pneumonia episodes between the two intervention
#' # groups in pneumovac data set. The test statistic we will use is the
#' # estimated log incidence rate ratio (IRR) from a Poisson GLM. (Note, it will
#' # take a few seconds to run 1,000 permutations)
#'
#' head(pneumovac) # visualize data
#' test <- permtest_glm(bpepisodes ~ spnvac, trtname = 'spnvac',
#'                      runit = 'randunit', family = poisson, data = pneumovac,
#'                      nperm = 1000, ncores = 2, seed = 445)
#' print(c(test$coef, test$pval))
#' # [1] -0.4466939  0.0560000
#' plot(test) # visualize Monte Carlo randomization distribution
#' @export
permtest_glm <- function(formula, trtname, runit, strat = NULL,
                         family = gaussian, data,
                         alternative = 'two-sided', theta = 0,
                         nperm = 1000, ncores = 1,
                         seed, quietly = T) {
  if (!(alternative %in% c('two-sided', 'greater', 'less')))
    stop("'alternative' must be one of c('two-sided', 'greater', 'less')")

  call <- match.call()
  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
  } else {
    foreach::registerDoSEQ()
  }
  if (!missing(seed)) {
    set.seed(seed)
  } else {
    seed <- NA
  }

  # add null value to data set
  data$null <- theta

  # fit glm
  m1 <- glm(formula = formula, family = family, data = data)
  obs1 <- as.numeric(coef(m1)[trtname])

  # permute based on runit
  perm.stat <- foreach::foreach(i = 2:nperm, .combine = c) %dorng% {
    data.tmp <- permute(data, trtname, runit, strat) # permuted data
    model.tmp <- update(m1, formula. = as.formula(paste0("~ . + offset(", trtname, " * null)")), data = data.tmp) # fit

    if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                            ceiling(nperm / 10)))
      cat(i, "of", nperm, "permutations complete\n")

    as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
  }

  comb.stat1 <- c(obs1, perm.stat)

  # prop more extreme than obs
  if (alternative == 'two-sided') {
    pval.perm <- mean(abs(comb.stat1) >= abs(obs1))
  } else if (alternative == 'less') {
    pval.perm <- mean(comb.stat1 < obs1)
  } else if (alternative == 'greater') {
    pval.perm <- mean(comb.stat1 > obs1)
  }

  out <- list(coef = obs1,
              pval = pval.perm,
              permCoefs = comb.stat1,
              call = call,
              args = list(
                trtname = trtname,
                runit = runit,
                strat = strat,
                nperm = nperm,
                ncores = ncores,
                seed = seed
              ))
  class(out) <- 'permtest'
  return(out)
}

#' @rdname permtest_glm
#' @export
permtest_survreg <- function(formula, trtname, runit, strat = NULL, data,
                             alternative = 'two-sided', theta = 0,
                             dist = "weibull", nperm = 1000, ncores = 1,
                             seed, quietly = T) {
  if (!(alternative %in% c('two-sided', 'greater', 'less')))
    stop("'alternative' must be one of c('two-sided', 'greater', 'less')")

  call <- match.call()
  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
  } else {
    foreach::registerDoSEQ()
  }
  if (!missing(seed)) {
    set.seed(seed)
  } else {
    seed <- NA
  }

  # add null value to data set
  data$null <- theta

  # fit
  m1 <- survival::survreg(formula = formula, data = data, dist = dist)
  obs1 <- as.numeric(coef(m1)[trtname])

  # permute based on runit
  perm.stat <- foreach::foreach(i = 2:nperm, .combine = c) %dorng% {
    data.tmp <- permute(data, trtname, runit, strat) # permuted data
    model.tmp <- update(m1, formula. = as.formula(paste0("~ . + offset(", trtname, " * null)")), data = data.tmp) # fit

    if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                            ceiling(nperm / 10)))
      cat(i, "of", nperm, "permutations complete\n")

    as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
  }

  comb.stat1 <- c(obs1, perm.stat)

  # prop more extreme than obs
  if (alternative == 'two-sided') {
    pval.perm <- mean(abs(comb.stat1) >= abs(obs1))
  } else if (alternative == 'less') {
    pval.perm <- mean(comb.stat1 < obs1)
  } else if (alternative == 'greater') {
    pval.perm <- mean(comb.stat1 > obs1)
  }

  out <- list(coef = obs1,
              pval = pval.perm,
              permCoefs = comb.stat1,
              call = call,
              args = list(
                trtname = trtname,
                runit = runit,
                strat = strat,
                nperm = nperm,
                ncores = ncores,
                seed = seed
              ))
  class(out) <- 'permtest'
  return(out)
}

#' @rdname permtest_glm
#' @export
permtest_coxph <- function(formula, trtname, runit, strat = NULL, data,
                           alternative = 'two-sided', theta = 0,
                           nperm = 1000, ncores = 1, seed, quietly = T) {
  if (!(alternative %in% c('two-sided', 'greater', 'less')))
    stop("'alternative' must be one of c('two-sided', 'greater', 'less')")

  call <- match.call()
  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
  } else {
    foreach::registerDoSEQ()
  }
  if (!missing(seed)) {
    set.seed(seed)
  } else {
    seed <- NA
  }

  # add null value to data set
  data$null <- theta

  # fit
  m1 <- survival::coxph(formula = formula, data = data)
  obs1 <- as.numeric(coef(m1)[trtname])

  # permute based on runit
  perm.stat <- foreach::foreach(i = 2:nperm, .combine = c) %dorng% {
    data.tmp <- permute(data, trtname, runit, strat) # permuted data
    model.tmp <- update(m1, formula. = as.formula(paste0("~ . + offset(", trtname, " * null)")), data = data.tmp) # fit

    if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                            ceiling(nperm / 10)))
      cat(i, "of", nperm, "permutations complete\n")

    as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
  }

  comb.stat1 <- c(obs1, perm.stat)

  # prop more extreme than obs
  if (alternative == 'two-sided') {
    pval.perm <- mean(abs(comb.stat1) >= abs(obs1))
  } else if (alternative == 'less') {
    pval.perm <- mean(comb.stat1 < obs1)
  } else if (alternative == 'greater') {
    pval.perm <- mean(comb.stat1 > obs1)
  }

  out <- list(coef = obs1,
              pval = pval.perm,
              permCoefs = comb.stat1,
              call = call,
              args = list(
                trtname = trtname,
                runit = runit,
                strat = strat,
                nperm = nperm,
                ncores = ncores,
                seed = seed
              ))
  class(out) <- 'permtest'
  return(out)
}

#' Randomization test
#'
#' Carry out a randomization test for a treatment effect using a fitted model
#' object or user defined test statistic.
#'
#' @seealso \code{\link[permuter]{permci}} for a randomization-based CI
#'
#' @param f fitted model object or function. If \code{f} is a fitted model
#' object, then the coefficient corresponding to \code{trtname} is used as the
#' test statistic. If \code{f} is a function, it must be defined such that when
#' applied to \code{data}, it returns the observed univariate test statistic.
#' This function can be as simple or complex as desired as long as its one input
#' argument is a data frame structured the same as \code{data}.
#' @inheritParams permute
#' @param data a data frame containing the variables necessary for the function
#' \code{f}. This argument is passed to \code{f}.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two-sided" (default), "greater", or "less".
#' @param theta a number indicating the null treatment effect value of interest
#' for the randomization test
#' @param nperm number of permutations for randomization test
#' @param ncores number of cores to use for computation. If ncores > 1, permtest
#' runs in parallel.
#' @param seed a numerical seed to use, passed to \code{\link[base]{set.seed}}
#' (if \code{ncores == 1}) or \code{\link[doRNG]{registerDoRNG}} (if
#' \code{ncores > 1}).
#' @param quietly logical; if TRUE (and if ncores == 1), status updates will be
#' printed to Console otherwise, suppress updates.
#' @export
permtest <- function(f, trtname, runit, strat = NULL, data,
                     alternative = 'two-sided', theta = 0,
                     nperm = 1000, ncores = 1, seed, quietly = T) {
  if (!(alternative %in% c('two-sided', 'greater', 'less')))
    stop("'alternative' must be one of c('two-sided', 'greater', 'less')")

  call <- match.call()
  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
  } else {
    foreach::registerDoSEQ()
  }
  if (!missing(seed)) {
    set.seed(seed)
  } else {
    seed <- NA
  }

  # add null value to data set
  data$null <- theta

  # calculate test statistic for obs data
  if (is.function(f)) {
    obs1 <- f(data)
  } else {
    obs1 <- as.numeric(coef(f)[trtname])
  }

  # permute based on runit
  perm.stat <- foreach::foreach(i = 2:nperm, .combine = c) %dorng% {
    data.tmp <- permute(data, trtname, runit, strat) # permuted data
    if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                            ceiling(nperm / 10)))
      cat(i, "of", nperm, "permutations complete\n")

    if (is.function(f)) {
      if (theta != 0)
        stop("Non-zero null only supported when 'f' is a fitted model object")
      f(data.tmp) # return test statistic for perm data
    } else {
      model.tmp <- update(f, formula. = as.formula(paste0("~ . + offset(", trtname, " * null)")), data = data.tmp) # fit
      as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
    }
  }

  comb.stat1 <- c(obs1, perm.stat)

  if (alternative == 'two-sided') {
    pval.perm <- mean(abs(comb.stat1) >= abs(obs1)) # prop more extreme than obs
  } else if (alternative == 'less') {
    pval.perm <- mean(comb.stat1 < obs1)
  } else if (alternative == 'greater') {
    pval.perm <- mean(comb.stat1 > obs1)
  }

  out <- list(coef = obs1,
              pval = pval.perm,
              permCoefs = comb.stat1,
              call = call,
              args = list(
                trtname = trtname,
                runit = runit,
                strat = strat,
                nperm = nperm,
                ncores = ncores,
                seed = seed
              ))
  class(out) <- 'permtest'
  return(out)
}
