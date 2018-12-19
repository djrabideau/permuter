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
#'   \item \code{permtest_ic_sp}: randomization test based on
#'   \code{\link[icenReg]{ic_sp}}
#'   \item \code{permtest_survreg}: randomization test based on
#'   \code{\link[survival]{survreg}}
#'   \item \code{permtest_coxph}: randomization test based on
#'   \code{\link[survival]{coxph}}
#' }
#' To ensure correct specification of the parameters passed to the models above
#' (e.g. \code{formula} in \code{\link[icenReg]{ic_sp}}), please refer to their
#' documentation.
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
#' @param nperm number of permutations for randomization test
#' @param ncores number of cores to use for computation. If ncores > 1, permtest
#' runs in parallel.
#' @param quietly logical; if TRUE (and if ncores == 1), status updates will be
#' printed to Console otherwise, suppress updates.
#' @param dist assumed distribution for y variable. If the argument is a
#' character string, then it is assumed to name an element from
#' \code{\link[survival]{survreg.distributions}}. These include "weibull",
#' "exponential", "gaussian", "logistic","lognormal" and "loglogistic".
#' Otherwise, it is assumed to be a user defined list conforming to the format
#' described in \code{\link[survival]{survreg.distributions}}.
#' @importFrom foreach %dopar%
#' @export
permtest_glm <- function(formula, trtname, runit, strat = NULL,
                         family = gaussian, data, nperm = 999, ncores = 1,
                         quietly = F) {
  # fit glm
  m1 <- glm(formula = formula, family = family, data = data)
  obs1 <- as.numeric(coef(m1)[trtname])

  perm.stat <- rep(0, nperm)

  # permute based on runit
  doMC::registerDoMC(ncores)
  perm.stat <- foreach::foreach(i = 1:nperm, .combine = c) %dopar% {
    data.tmp <- permute(data, trtname, runit, strat) # permuted data
    model.tmp <- glm(formula = formula, family = family, data = data.tmp) # fit
    if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                            ceiling(nperm / 10)))
      cat(i, "of", nperm, "permutations complete\n")

    as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
  }

  comb.stat1 <- c(obs1, perm.stat)

  pval.perm <- mean(abs(comb.stat1) >= abs(obs1)) # prop more extreme than obs

  out <- c(obs1, pval.perm)
  names(out) <- c(trtname, "p.perm")
  return(out)
}


#' @rdname permtest_glm
#' @export
permtest_ic_sp <- function(formula, trtname, runit, strat = NULL, data,
                           nperm = 999, ncores = 1, quietly = F) {
  # fit ic_sp
  m1 <- icenReg::ic_sp(formula = formula, data = data)
  obs1 <- as.numeric(coef(m1)[trtname])

  perm.stat <- rep(0, nperm)

  # permute based on runit
  doMC::registerDoMC(ncores)
  perm.stat <- foreach::foreach(i = 1:nperm, .combine = c) %dopar% {
    data.tmp <- permute(data, trtname, runit, strat) # permuted data
    model.tmp <- icenReg::ic_sp(formula = formula, data = data.tmp) # fit

    if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                            ceiling(nperm / 10)))
      cat(i, "of", nperm, "permutations complete\n")

    as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
  }

  comb.stat1 <- c(obs1, perm.stat)

  pval.perm <- mean(abs(comb.stat1) >= abs(obs1)) # prop more extreme than obs

  out <- c(obs1, pval.perm)
  names(out) <- c(trtname, "p.perm")
  return(out)
}

#' @rdname permtest_glm
#' @export
permtest_survreg <- function(formula, trtname, runit, strat = NULL, data,
                             dist = "weibull", nperm = 999, ncores = 1,
                             quietly = F) {
  # fit
  m1 <- survival::survreg(formula = formula, data = data, dist = dist)
  obs1 <- as.numeric(coef(m1)[trtname])

  perm.stat <- rep(0, nperm)

  # permute based on runit
  doMC::registerDoMC(ncores)
  perm.stat <- foreach::foreach(i = 1:nperm, .combine = c) %dopar% {
    data.tmp <- permute(data, trtname, runit, strat) # permuted data
    model.tmp <- survival::survreg(formula = formula, data = data.tmp,
                                   dist = dist) # fit

    if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                            ceiling(nperm / 10)))
      cat(i, "of", nperm, "permutations complete\n")

    as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
  }

  comb.stat1 <- c(obs1, perm.stat)

  pval.perm <- mean(abs(comb.stat1) >= abs(obs1)) # prop more extreme than obs

  out <- c(obs1, pval.perm)
  names(out) <- c(trtname, "p.perm")
  return(out)
}

#' @rdname permtest_glm
#' @export
permtest_coxph <- function(formula, trtname, runit, strat = NULL, data,
                             nperm = 999, ncores = 1, quietly = F) {
  # fit
  m1 <- survival::coxph(formula = formula, data = data)
  obs1 <- as.numeric(coef(m1)[trtname])

  perm.stat <- rep(0, nperm)

  # permute based on runit
  doMC::registerDoMC(ncores)
  perm.stat <- foreach::foreach(i = 1:nperm, .combine = c) %dopar% {
    data.tmp <- permute(data, trtname, runit, strat) # permuted data
    model.tmp <- survival::coxph(formula = formula, data = data.tmp) # fit

    if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                            ceiling(nperm / 10)))
      cat(i, "of", nperm, "permutations complete\n")

    as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
  }

  comb.stat1 <- c(obs1, perm.stat)

  pval.perm <- mean(abs(comb.stat1) >= abs(obs1)) # prop more extreme than obs

  out <- c(obs1, pval.perm)
  names(out) <- c(trtname, "p.perm")
  return(out)
}
