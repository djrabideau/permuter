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
#'   \item \code{permtest_ic_sp}: randomization CI based on
#'   \code{\link[icenReg]{ic_sp}}
#' }
#' To ensure correct specification of the parameters passed to the models above
#' (e.g. \code{formula} in \code{\link[icenReg]{ic_sp}}), please refer to their
#' documentation.
#'
#' @inheritParams permtest_glm
#' @param level two-sided confidence level (e.g. level = 0.95 for 95\% CI)
#' @param ncores number of cores to use for computation. If ncores > 1, lower
#' and upper bound search procedures run in parallel.
#' @param quietly logical; if TRUE (and if ncores == 1), status updates will be
#' printed to Console otherwise, suppress updates.
#' @importFrom foreach %dopar%
#' @export
permci_glm <- function(formula, trtname, runit, strat = NULL,
                       family = gaussian, data, nperm = 999, level = 0.95,
                       quietly = F, ncores = 1) {
  data[, paste0(trtname, ".obs")] <- data[, trtname] # obs trt for offset

  alpha <- 1 - level

  # get lower/upper with GLM norm approx
  m1 <- glm(formula = formula, family = family, data = data)
  Vcov <- vcov(m1, useScale = FALSE)
  obs1 <- as.numeric(coef(m1)[trtname])
  trt.se <- sqrt(Vcov[trtname, trtname])
  lower <- obs1 - qnorm(1 - alpha / 2) * trt.se
  upper <- obs1 + qnorm(1 - alpha / 2) * trt.se

  # initialize at lower/upper
  data$low <- low <- lower
  data$up <- up <- upper

  # if more than 1 core, run lower/upper in parallel
  doMC::registerDoMC(ncores)

  # invert test for CI using offset approach
  trace <- foreach::foreach(j = 1:min(ncores, 2), .combine = cbind) %dopar% {
    if (j == 1 | ncores == 1) {
      # search for lower
      low.vec <- c(low, rep(NA, nperm - 1))
      formula.tmp <- update(formula,
                  as.formula(paste0("~ . + offset(", trtname, ".obs * low)")))
      for (i in 1:nperm) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- glm(formula = formula.tmp, family = family,
                         data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - low) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        low <- update_rm(init = low, thetahat = obs1, t, tstar, alpha, i,
                         bound = "lower")
        data.tmp$low <- low
        low.vec[i] <- low

        if (!quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                  ceiling(nperm / 10)))
          cat(i, "of", nperm, "permutations complete\n")
      }
      if (ncores == 1 & !quietly) cat("lower bound complete\n")
      low.vec
    }

    if (j == 2 | ncores == 1) {
      # search for upper
      up.vec <- c(up, rep(NA, nperm - 1))
      formula.tmp <- update(formula,
                    as.formula(paste0("~ . + offset(", trtname, ".obs * up)")))
      for (i in 1:nperm) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- glm(formula = formula.tmp, family = family,
                         data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - up) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        up <- update_rm(init = up, thetahat = obs1, t, tstar, alpha, i,
                        bound = "upper")
        data.tmp$up <- up
        up.vec[i] <- up

        if (!quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                  ceiling(nperm / 10)))
          cat(i, "of", nperm, "permutations complete\n")
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

  dimnames(trace)[[2]] <- c("lower", "upper")
  return(list(ci = c(trace[nperm, 1], trace[nperm, 2]), trace = trace))
}


#' @rdname permci_glm
#' @export
permci_ic_sp <- function(formula, trtname, runit, strat = NULL, data,
                         nperm = 999, level = 0.95, quietly = F, ncores = 1) {
  data[, paste0(trtname, ".obs")] <- data[, trtname] # obs trt for offset

  alpha <- 1 - level

  # get lower/upper with weibull model for interval censored
  m1 <- survreg(formula = formula, data = data)
  Vcov <- vcov(m1, useScale = FALSE)
  obs1 <- - (as.numeric(coef(m1)[trtname])) # weibull flips effect
  trt.se <- sqrt(Vcov[trtname, trtname])
  lower <- obs1 - qnorm(1 - alpha / 2) * trt.se
  upper <- obs1 + qnorm(1 - alpha / 2) * trt.se

  # reset m1 and obs1 corresponding to ic_sp
  m1 <- ic_sp(formula = formula, data = data)
  obs1 <- as.numeric(coef(m1)[trtname])

  # initialize at lower/upper
  data$low <- low <- lower
  data$up <- up <- upper

  # if more than 1 core, run lower/upper in parallel
  doMC::registerDoMC(ncores)

  # search for lower
  trace <- foreach::foreach(j = 1:min(ncores, 2), .combine = cbind) %dopar% {
    if (j == 1 | ncores == 1) {
      low.vec <- c(low, rep(NA, nperm - 1))
      formula.tmp <- update(formula,
                  as.formula(paste0("~ . + offset(", trtname, ".obs * low)")))
      for (i in 1:nperm) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- ic_sp(formula = formula.tmp, data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - low) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        low <- update_rm(init = low, thetahat = obs1, t, tstar, alpha, i,
                         bound = "lower")
        data.tmp$low <- low
        low.vec[i] <- low

        if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                              ceiling(nperm / 10)))
          cat(i, "of", nperm, "permutations complete\n")
      }
      if (ncores == 1 & !quietly) cat("lower bound complete\n")
      low.vec
    }

    if (j == 2 | ncores == 1) {

      # search for upper
      up.vec <- c(up, rep(NA, nperm - 1))
      formula.tmp <- update(formula,
                    as.formula(paste0("~ . + offset(", trtname, ".obs * up)")))
      for (i in 1:nperm) {

        # permute based on runit
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- ic_sp(formula = formula.tmp, data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
        tstar <- (obs1 - up) # tx effect estimate from original permutation

        # update using Robbins-Monro step
        up <- update_rm(init = up, thetahat = obs1, t, tstar, alpha, i,
                        bound = "upper")
        data.tmp$up <- up
        up.vec[i] <- up

        if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                              ceiling(nperm / 10)))
          cat(i, "of", nperm, "permutations complete\n")
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

  dimnames(trace)[[2]] <- c("lower", "upper")
  return(list(ci = c(trace[nperm, 1], trace[nperm, 2]), trace = trace))
}
