#' Randomization CI for regression
#'
#' Calculate a randomization confidence interval (CI) for a regression parameter.
#'
#' These functions are used to calculate randomization confidence intervals (CI)
#' for a regression parameter. These CIs correspond to inverting randomization
#' tests by using an offset to test non-zero "null" values. To invert the
#' randomization test, these functions use a computationally efficient CI
#' algorithm proposed by
#' \href{http://doi.org/10.2307/2532852}{Garthwaite (1996)}, which is based on
#' the Robbins-Monro search process.
#'
#' Different functions
#' corrrespond to different regression models:
#' \itemize{
#'   \item \code{permci_glm}: randomization CI based on
#'   \code{\link[stats]{glm}}
#'   \item \code{permci_ic_sp}: randomization CI based on
#'   \code{\link[icenReg]{ic_sp}}
#'   \item \code{permci_survival::survreg}: randomization CI based on
#'   \code{\link[survival]{survreg}}
#'   \item \code{permci_coxph}: randomization CI based on
#'   \code{\link[survival]{coxph}}
#' }
#' To ensure correct specification of the parameters passed to the models above
#' (e.g. \code{formula} in \code{\link[icenReg]{ic_sp}}), please refer to their
#' documentation.
#'
#' @inheritParams permtest_glm
#' @param level two-sided confidence level (e.g. level = 0.95 for 95\% CI)
#' @param init vector of initial values for CI, with lower bound as
#' first element and upper bound as second. If \code{init} not provided,
#' initial bounds are based on \code{initmethod}.
#' @param initmethod character; indicates the method to be used for initial
#' values for CI. If "asymp", initial bounds are based on asymptotic
#' approximation (e.g. Wald CI for GLM). If "perm" (default), initial bounds are based
#' on the permutation approach used in Garthwaite (1996) with \eqn{\hat{\theta}
#' \pm \{(t_2 - t_1)/2\}}, where \eqn{t_1} and \eqn{t_2} denote the second
#' smallest and second largest estimates from the permutation test.
#' @param nperm number of permutations for each randomization CI bound
#' @param nburn number of ``burn-in'' permutations. I.e. algorithm will start at
#' \code{init} bounds, run for \code{nburn} permutations, then restart algorithm
#' at latest estimates from ``burn-in'' phase and run for another \code{nperm}
#' permutations until the final CI estimates are reached. Increasing
#' \code{nburn} may help convergence if \code{init} CI bounds are poor.
#' @param ncores number of cores to use for computation. If \code{ncores} > 1,
#' lower and upper bound search procedures run in parallel across 2 cores.
#' @param quietly logical; if TRUE (and if ncores == 1), status updates will be
#' printed to Console otherwise, suppress updates.
#' @param ... optional arguments to \code{\link[permuter]{update_rm}}, e.g.
#' \code{m} controls the initial magnitude of the steps
#' @export
permci_glm <- function(formula, trtname, runit, strat = NULL,
                       family = gaussian, data, nperm = 1000, nburn = 0,
                       level = 0.95, init, initmethod = 'perm',
                       ncores = 1, seed, quietly = F, ...) {
  call <- match.call()
  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
    if (!missing(seed))
      doRNG::registerDoRNG(seed)
  } else {
    if (!missing(seed))
      set.seed(seed)
  }

  data[, paste0(trtname, ".obs")] <- data[, trtname] # obs trt for offset

  alpha <- 1 - level

  # get lower/upper with GLM norm approx
  m1 <- glm(formula = formula, family = family, data = data)
  Vcov <- vcov(m1, useScale = FALSE)
  obs1 <- as.numeric(coef(m1)[trtname])
  trt.se <- sqrt(Vcov[trtname, trtname])
  lower <- obs1 - qnorm(1 - alpha / 2) * trt.se
  upper <- obs1 + qnorm(1 - alpha / 2) * trt.se

  if (missing(init)) {
    if (initmethod == 'asymp') {
      # initialize at asymptotic lower/upper
      data$low <- low <- lower
      data$up <- up <- upper
    } else if (initmethod == 'perm') {
      # initialize using quick randomization test of H0: theta = obs1,
      # as recommended in Garthwaite (1996)
      g.alpha <- alpha / 2 # alpha as defined in Garthwaite paper
      nperm_init <- ceiling((2 - g.alpha) / g.alpha)
      data$obs1 <- obs1
      formula.tmp <- update(formula,
                  as.formula(paste0("~ . + offset(", trtname, ".obs * obs1)")))
      perm.stat <- foreach::foreach(i = 1:nperm_init, .combine = c) %dopar% {
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- glm(formula = formula.tmp, family = family,
                         data = data.tmp) # fit
        as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
      }
      t1 <- sort(perm.stat)[2] # 2nd to smallest
      t2 <- sort(perm.stat)[nperm_init - 1] # 2nd to largest
      data$low <- low <- obs1 - ((t2 - t1) / 2)
      data$up <- up <- obs1 + ((t2 - t1) / 2)
    } else {
      print(initmethod)
      stop("'initmethod' not recognized")
    }
  } else {
    # initialize at given initial values
    data$low <- low <- init[1]
    data$up <- up <- init[2]
  }

  inits <- c(low, up)

  # if more than 1 core, run lower/upper in parallel
  # invert test for CI using offset approach
  trace <- foreach::foreach(j = 1:min(ncores, 2), .combine = cbind) %dopar% {
    if (j == 1 | ncores == 1) {
      # search for lower
      low.vec <- rep(NA, nperm + nburn)
      formula.tmp <- update(formula,
                  as.formula(paste0("~ . + offset(", trtname, ".obs * low)")))
      for (i in 1:(nperm + nburn)) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- glm(formula = formula.tmp, family = family,
                         data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - low) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        ii <- i - as.numeric(i > nburn) * nburn # reset i <- 1 after nburn perms
        low <- update_rm(init = low, thetahat = obs1, t, tstar, alpha, ii,
                         bound = "lower", ...)
        data$low <- low
        low.vec[i] <- low

        if (ncores == 1 & !quietly & i %in% seq(ceiling((nperm + nburn) / 10), (nperm + nburn),
                                  ceiling((nperm + nburn) / 10)))
          cat(i, "of", (nperm + nburn), "permutations complete\n")
      }
      if (ncores == 1 & !quietly) cat("lower bound complete\n")
      low.vec
    }

    if (j == 2 | ncores == 1) {
      # search for upper
      up.vec <- rep(NA, nperm + nburn)
      formula.tmp <- update(formula,
                    as.formula(paste0("~ . + offset(", trtname, ".obs * up)")))
      for (i in 1:(nperm + nburn)) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- glm(formula = formula.tmp, family = family,
                         data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - up) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        ii <- i - as.numeric(i > nburn) * nburn # reset i <- 1 after nburn perms
        up <- update_rm(init = up, thetahat = obs1, t, tstar, alpha, ii,
                        bound = "upper", ...)
        data$up <- up
        up.vec[i] <- up

        if (ncores == 1 & !quietly & i %in% seq(ceiling((nperm + nburn) / 10), (nperm + nburn),
                                  ceiling((nperm + nburn) / 10)))
          cat(i, "of", (nperm + nburn), "permutations complete\n")
      }
      if (ncores == 1 & !quietly) cat("upper bound complete\n")
      up.vec
    }

    # return these values
    if (ncores == 1) {
      cbind(low.vec, up.vec)
    } else if (j == 1) {
      low.vec
    } else {
      up.vec
    }
  } # end foreach

  if (missing(seed)) seed <- NA
  dimnames(trace)[[2]] <- c("lower", "upper")
  out <- list(ci = c(trace[nperm + nburn, 1], trace[nperm + nburn, 2]),
              trace = trace,
              init = inits,
              call = call,
              args = list(
                trtname = trtname,
                runit = runit,
                strat = strat,
                nperm = nperm,
                nburn = nburn,
                level = level,
                initmethod = initmethod,
                ncores = ncores,
                seed = seed
              ))
  class(out) <- 'permci'
  return(out)
}


#' @rdname permci_glm
#' @export
permci_ic_sp <- function(formula, trtname, runit, strat = NULL, data,
                         nperm = 1000, nburn = 0,
                         level = 0.95, init, initmethod = 'perm',
                         ncores = 1, seed, quietly = F, ...) {
  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
    if (!missing(seed))
      doRNG::registerDoRNG(seed)
  } else {
    if (!missing(seed))
      set.seed(seed)
  }

  data[, paste0(trtname, ".obs")] <- data[, trtname] # obs trt for offset

  alpha <- 1 - level

  # get lower/upper with weibull model for interval censored
  m1 <- survival::survreg(formula = formula, data = data)
  Vcov <- vcov(m1, useScale = FALSE)
  obs1 <- as.numeric(coef(m1)[trtname])
  trt.se <- sqrt(Vcov[trtname, trtname])
  lower.sr <- obs1 - qnorm(1 - alpha / 2) * trt.se
  upper.sr <- obs1 + qnorm(1 - alpha / 2) * trt.se
  # re-parameterize from survival::survreg to ic_sp
  lower <- - upper.sr / m1$scale
  upper <- - lower.sr / m1$scale

  # reset m1 and obs1 corresponding to ic_sp
  m1 <- icenReg::ic_sp(formula = formula, data = data)
  obs1 <- as.numeric(coef(m1)[trtname])

  if (missing(init)) {
    if (initmethod == 'asymp') {
      # initialize at asymptotic lower/upper
      data$low <- low <- lower
      data$up <- up <- upper
    } else if (initmethod == 'perm') {
      # initialize using quick randomization test of H0: theta = obs1,
      # as recommended in Garthwaite (1996)
      g.alpha <- alpha / 2 # alpha as defined in Garthwaite paper
      nperm_init <- ceiling((2 - g.alpha) / g.alpha)
      data$obs1 <- obs1
      formula.tmp <- update(formula,
                  as.formula(paste0("~ . + offset(", trtname, ".obs * obs1)")))
      perm.stat <- foreach::foreach(i = 1:nperm_init, .combine = c) %dopar% {
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- icenReg::ic_sp(formula = formula.tmp, data = data.tmp) # fit
        as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
      }
      t1 <- sort(perm.stat)[2] # 2nd to smallest
      t2 <- sort(perm.stat)[nperm_init - 1] # 2nd to largest
      data$low <- low <- obs1 - ((t2 - t1) / 2)
      data$up <- up <- obs1 + ((t2 - t1) / 2)
    } else {
      print(initmethod)
      stop("'initmethod' not recognized")
    }
  } else {
    # initialize at given initial values
    data$low <- low <- init[1]
    data$up <- up <- init[2]
  }

  inits <- c(low, up)

  # if more than 1 core, run lower/upper in parallel
  # search for lower
  trace <- foreach::foreach(j = 1:min(ncores, 2), .combine = cbind) %dopar% {
    if (j == 1 | ncores == 1) {
      low.vec <- rep(NA, nperm + nburn)
      formula.tmp <- update(formula,
                  as.formula(paste0("~ . + offset(", trtname, ".obs * low)")))
      for (i in 1:(nperm + nburn)) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- icenReg::ic_sp(formula = formula.tmp, data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - low) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        ii <- i - as.numeric(i > nburn) * nburn # reset i <- 1 after nburn perms
        low <- update_rm(init = low, thetahat = obs1, t, tstar, alpha, ii,
                         bound = "lower", ...)
        data$low <- low
        low.vec[i] <- low

        if (ncores == 1 & !quietly & i %in% seq(ceiling((nperm + nburn) / 10), (nperm + nburn),
                                              ceiling((nperm + nburn) / 10)))
          cat(i, "of", (nperm + nburn), "permutations complete\n")
      }
      if (ncores == 1 & !quietly) cat("lower bound complete\n")
      low.vec
    }

    if (j == 2 | ncores == 1) {

      # search for upper
      up.vec <- rep(NA, nperm + nburn)
      formula.tmp <- update(formula,
                    as.formula(paste0("~ . + offset(", trtname, ".obs * up)")))
      for (i in 1:(nperm + nburn)) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- icenReg::ic_sp(formula = formula.tmp, data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - up) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        ii <- i - as.numeric(i > nburn) * nburn # reset i <- 1 after nburn perms
        up <- update_rm(init = up, thetahat = obs1, t, tstar, alpha, ii,
                        bound = "upper", ...)
        data$up <- up
        up.vec[i] <- up

        if (ncores == 1 & !quietly & i %in% seq(ceiling((nperm + nburn) / 10), (nperm + nburn),
                                              ceiling((nperm + nburn) / 10)))
          cat(i, "of", (nperm + nburn), "permutations complete\n")
      }
      if (ncores == 1 & !quietly) cat("upper bound complete\n")
      up.vec
    }

    # return these values
    if (ncores == 1) {
      cbind(low.vec, up.vec)
    } else if (j == 1) {
      low.vec
    } else {
      up.vec
    }
  } # end foreach

  if (missing(seed)) seed <- NA
  dimnames(trace)[[2]] <- c("lower", "upper")
  out <- list(ci = c(trace[nperm + nburn, 1], trace[nperm + nburn, 2]),
              trace = trace,
              init = inits,
              call = call,
              args = list(
                trtname = trtname,
                runit = runit,
                strat = strat,
                nperm = nperm,
                nburn = nburn,
                level = level,
                initmethod = initmethod,
                ncores = ncores,
                seed = seed
              ))
  class(out) <- 'permci'
  return(out)
}


#' @rdname permci_glm
#' @export
permci_survreg <- function(formula, trtname, runit, strat = NULL, data,
                           dist = "weibull", nperm = 1000, nburn = 0,
                           level = 0.95, init, initmethod = 'perm',
                           ncores = 1, seed, quietly = F, ...) {
  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
    if (!missing(seed))
      doRNG::registerDoRNG(seed)
  } else {
    if (!missing(seed))
      set.seed(seed)
  }

  data[, paste0(trtname, ".obs")] <- data[, trtname] # obs trt for offset

  alpha <- 1 - level

  # get lower/upper with survival::survreg
  m1 <- survival::survreg(formula = formula, data = data, dist = dist)
  Vcov <- vcov(m1, useScale = FALSE)
  obs1 <- as.numeric(coef(m1)[trtname])
  trt.se <- sqrt(Vcov[trtname, trtname])
  lower <- obs1 - qnorm(1 - alpha / 2) * trt.se
  upper <- obs1 + qnorm(1 - alpha / 2) * trt.se

  if (missing(init)) {
    if (initmethod == 'asymp') {
      # initialize at asymptotic lower/upper
      data$low <- low <- lower
      data$up <- up <- upper
    } else if (initmethod == 'perm') {
      # initialize using quick randomization test of H0: theta = obs1,
      # as recommended in Garthwaite (1996)
      g.alpha <- alpha / 2 # alpha as defined in Garthwaite paper
      nperm_init <- ceiling((2 - g.alpha) / g.alpha)
      data$obs1 <- obs1
      formula.tmp <- update(formula,
                  as.formula(paste0("~ . + offset(", trtname, ".obs * obs1)")))
      perm.stat <- foreach::foreach(i = 1:nperm_init, .combine = c) %dopar% {
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- survival::survreg(formula = formula.tmp, data = data.tmp, dist = dist) # fit
        as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
      }
      t1 <- sort(perm.stat)[2] # 2nd to smallest
      t2 <- sort(perm.stat)[nperm_init - 1] # 2nd to largest
      data$low <- low <- obs1 - ((t2 - t1) / 2)
      data$up <- up <- obs1 + ((t2 - t1) / 2)
    } else {
      print(initmethod)
      stop("'initmethod' not recognized")
    }
  } else {
    # initialize at given initial values
    data$low <- low <- init[1]
    data$up <- up <- init[2]
  }

  inits <- c(low, up)

  # if more than 1 core, run lower/upper in parallel
  # search for lower
  trace <- foreach::foreach(j = 1:min(ncores, 2), .combine = cbind) %dopar% {
    if (j == 1 | ncores == 1) {
      low.vec <- rep(NA, nperm + nburn)
      formula.tmp <- update(formula,
                  as.formula(paste0("~ . + offset(", trtname, ".obs * low)")))
      for (i in 1:(nperm + nburn)) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- survival::survreg(formula = formula.tmp, data = data.tmp,
                             dist = dist) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - low) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        ii <- i - as.numeric(i > nburn) * nburn # reset i <- 1 after nburn perms
        low <- update_rm(init = low, thetahat = obs1, t, tstar, alpha, ii,
                         bound = "lower", ...)
        data$low <- low
        low.vec[i] <- low

        if (ncores == 1 & !quietly & i %in% seq(ceiling((nperm + nburn) / 10), (nperm + nburn),
                                              ceiling((nperm + nburn) / 10)))
          cat(i, "of", (nperm + nburn), "permutations complete\n")
      }
      if (ncores == 1 & !quietly) cat("lower bound complete\n")
      low.vec
    }

    if (j == 2 | ncores == 1) {

      # search for upper
      up.vec <- rep(NA, nperm + nburn)
      formula.tmp <- update(formula,
                    as.formula(paste0("~ . + offset(", trtname, ".obs * up)")))
      for (i in 1:(nperm + nburn)) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- survival::survreg(formula = formula.tmp, data = data.tmp,
                             dist = dist) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - up) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        ii <- i - as.numeric(i > nburn) * nburn # reset i <- 1 after nburn perms
        up <- update_rm(init = up, thetahat = obs1, t, tstar, alpha, ii,
                        bound = "upper", ...)
        data$up <- up
        up.vec[i] <- up

        if (ncores == 1 & !quietly & i %in% seq(ceiling((nperm + nburn) / 10), (nperm + nburn),
                                              ceiling((nperm + nburn) / 10)))
          cat(i, "of", (nperm + nburn), "permutations complete\n")
      }
      if (ncores == 1 & !quietly) cat("upper bound complete\n")
      up.vec
    }

    # return these values
    if (ncores == 1) {
      cbind(low.vec, up.vec)
    } else if (j == 1) {
      low.vec
    } else {
      up.vec
    }
  } # end foreach

  if (missing(seed)) seed <- NA
  dimnames(trace)[[2]] <- c("lower", "upper")
  out <- list(ci = c(trace[nperm + nburn, 1], trace[nperm + nburn, 2]),
              trace = trace,
              init = inits,
              call = call,
              args = list(
                trtname = trtname,
                runit = runit,
                strat = strat,
                nperm = nperm,
                nburn = nburn,
                level = level,
                initmethod = initmethod,
                ncores = ncores,
                seed = seed
              ))
  class(out) <- 'permci'
  return(out)
}


#' @rdname permci_glm
#' @export
permci_coxph <- function(formula, trtname, runit, strat = NULL, data,
                           nperm = 1000, nburn = 0,
                           level = 0.95, init, initmethod = 'perm',
                           ncores = 1, seed, quietly = F, ...) {
  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
    if (!missing(seed))
      doRNG::registerDoRNG(seed)
  } else {
    if (!missing(seed))
      set.seed(seed)
  }

  data[, paste0(trtname, ".obs")] <- data[, trtname] # obs trt for offset

  alpha <- 1 - level

  # get lower/upper with survival::coxph
  m1 <- survival::coxph(formula = formula, data = data)
  Vcov <- vcov(m1, useScale = FALSE)
  obs1 <- as.numeric(coef(m1)[trtname])
  trt.se <- sqrt(Vcov[trtname, trtname])
  lower <- obs1 - qnorm(1 - alpha / 2) * trt.se
  upper <- obs1 + qnorm(1 - alpha / 2) * trt.se

  if (missing(init)) {
    if (initmethod == 'asymp') {
      # initialize at asymptotic lower/upper
      data$low <- low <- lower
      data$up <- up <- upper
    } else if (initmethod == 'perm') {
      # initialize using quick randomization test of H0: theta = obs1,
      # as recommended in Garthwaite (1996)
      g.alpha <- alpha / 2 # alpha as defined in Garthwaite paper
      nperm_init <- ceiling((2 - g.alpha) / g.alpha)
      data$obs1 <- obs1
      formula.tmp <- update(formula,
                            as.formula(paste0("~ . + offset(", trtname, ".obs * obs1)")))
      perm.stat <- foreach::foreach(i = 1:nperm_init, .combine = c) %dopar% {
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- survival::coxph(formula = formula.tmp, data = data.tmp) # fit
        as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
      }
      t1 <- sort(perm.stat)[2] # 2nd to smallest
      t2 <- sort(perm.stat)[nperm_init - 1] # 2nd to largest
      data$low <- low <- obs1 - ((t2 - t1) / 2)
      data$up <- up <- obs1 + ((t2 - t1) / 2)
    } else {
      print(initmethod)
      stop("'initmethod' not recognized")
    }
  } else {
    # initialize at given initial values
    data$low <- low <- init[1]
    data$up <- up <- init[2]
  }

  inits <- c(low, up)

  # if more than 1 core, run lower/upper in parallel
  # search for lower
  trace <- foreach::foreach(j = 1:min(ncores, 2), .combine = cbind) %dopar% {
    if (j == 1 | ncores == 1) {
      low.vec <- rep(NA, nperm + nburn)
      formula.tmp <- update(formula,
                            as.formula(paste0("~ . + offset(", trtname, ".obs * low)")))
      for (i in 1:(nperm + nburn)) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- survival::coxph(formula = formula.tmp, data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - low) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        ii <- i - as.numeric(i > nburn) * nburn # reset i <- 1 after nburn perms
        low <- update_rm(init = low, thetahat = obs1, t, tstar, alpha, ii,
                         bound = "lower", ...)
        data$low <- low
        low.vec[i] <- low

        if (ncores == 1 & !quietly & i %in% seq(ceiling((nperm + nburn) / 10), (nperm + nburn),
                                                ceiling((nperm + nburn) / 10)))
          cat(i, "of", (nperm + nburn), "permutations complete\n")
      }
      if (ncores == 1 & !quietly) cat("lower bound complete\n")
      low.vec
    }

    if (j == 2 | ncores == 1) {

      # search for upper
      up.vec <- rep(NA, nperm + nburn)
      formula.tmp <- update(formula,
                            as.formula(paste0("~ . + offset(", trtname, ".obs * up)")))
      for (i in 1:(nperm + nburn)) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- survival::coxph(formula = formula.tmp, data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - up) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        ii <- i - as.numeric(i > nburn) * nburn # reset i <- 1 after nburn perms
        up <- update_rm(init = up, thetahat = obs1, t, tstar, alpha, ii,
                        bound = "upper", ...)
        data$up <- up
        up.vec[i] <- up

        if (ncores == 1 & !quietly & i %in% seq(ceiling((nperm + nburn) / 10), (nperm + nburn),
                                                ceiling((nperm + nburn) / 10)))
          cat(i, "of", (nperm + nburn), "permutations complete\n")
      }
      if (ncores == 1 & !quietly) cat("upper bound complete\n")
      up.vec
    }

    # return these values
    if (ncores == 1) {
      cbind(low.vec, up.vec)
    } else if (j == 1) {
      low.vec
    } else {
      up.vec
    }
  } # end foreach

  if (missing(seed)) seed <- NA
  dimnames(trace)[[2]] <- c("lower", "upper")
  out <- list(ci = c(trace[nperm + nburn, 1], trace[nperm + nburn, 2]),
              trace = trace,
              init = inits,
              call = call,
              args = list(
                trtname = trtname,
                runit = runit,
                strat = strat,
                nperm = nperm,
                nburn = nburn,
                level = level,
                initmethod = initmethod,
                ncores = ncores,
                seed = seed
              ))
  class(out) <- 'permci'
  return(out)
}
