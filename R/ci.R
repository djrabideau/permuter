# function to calculate confidence intervals by inverting permutation test
# using an offset in GLM for non-zero "null" values.
# Algorithm: Garthwaite 1996 - Biometrics (uses Robbins-Monro search process)
# uses m instead of i as recommended in Garthwaite 1996

permci_glm <- function(formula, trt.name, runit, strat = NULL,
                       family = gaussian, data, nperm = 1000, level = 0.95,
                       quietly = F, ncores = 1) {
  # formula:          formula argument passed to glmer(),
  #                     e.g. y ~ trt + factor(period) + (1 | clusid)
  # trt.name:         character string specifying the name of randomized
  #                     treatment variable in data frame (variable to permute)
  # runit:            character string specifying the name of unit
  #                     of randomization in data frame
  # strat:            character string specifying the name of the variable
  #                     upon which randomization was stratified
  # family:           family argument passed to glmer(), e.g. gaussian, binomial
  # data:             data argument passed to glmer()
  # quietly:          T if you want some status updates, convergence notes, etc.
  # ncores:           if >1, then lower/upper run in parallel

  data[, paste0(trt.name, ".obs")] <- data[, trt.name] # obs trt for offset

  alpha <- 1 - level

  # get lower/upper with GLM norm approx
  m1 <- glm(formula = formula, family = family, data = data)
  Vcov <- vcov(m1, useScale = FALSE)
  obs1 <- as.numeric(coef(m1)[trt.name])
  trt.se <- sqrt(Vcov[trt.name, trt.name])
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
                  as.formula(paste0("~ . + offset(", trt.name, ".obs * low)")))
      for (i in 1:nperm) {

        # permute based on runit
        data.tmp <- permute(data, trt.name, runit, strat) # permuted data
        model.tmp <- glm(formula = formula.tmp, family = family,
                         data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trt.name]) # return tx effect estimate
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
                    as.formula(paste0("~ . + offset(", trt.name, ".obs * up)")))
      for (i in 1:nperm) {

        # permute based on runit
        data.tmp <- permute(data, trt.name, runit, strat) # permuted data
        model.tmp <- glm(formula = formula.tmp, family = family,
                         data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trt.name]) # return tx effect estimate
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


# function to calculate confidence intervals by inverting permutation test
# using an offset in Cox PH model with interval censoring for non-zero "null"
# values. Algorithm: Garthwaite 1996 - Biometrics (uses Robbins-Monro search
# process) uses m instead of i as recommended in Garthwaite 1996

permci_ic_sp <- function(formula, trt.name, runit, strat = NULL, data,
                         nperm = 999, level = 0.95, quietly = F, ncores = 1) {
  # formula:          formula argument passed to ic_sp(), e.g. cbind(l, u) ~ trt
  # trt.name:         character string specifying the name of randomized
  #                     treatment variable in data frame (variable to permute)
  # runit:            character string specifying the name of unit
  #                     of randomization in data frame
  # strat:            character string specifying the name of the variable
  #                     upon which randomization was stratified
  # data:             data argument passed to ic_sp()
  # quietly:          T if you want some status updates, convergence notes, etc.
  # ncores:           if >1, then lower/upper run in parallel

  data[, paste0(trt.name, ".obs")] <- data[, trt.name] # obs trt for offset

  alpha <- 1 - level

  # get lower/upper with weibull model for interval censored
  m1 <- survreg(formula = formula, data = data)
  Vcov <- vcov(m1, useScale = FALSE)
  obs1 <- - (as.numeric(coef(m1)[trt.name])) # weibull flips effect
  trt.se <- sqrt(Vcov[trt.name, trt.name])
  lower <- obs1 - qnorm(1 - alpha / 2) * trt.se
  upper <- obs1 + qnorm(1 - alpha / 2) * trt.se

  # reset m1 and obs1 corresponding to ic_sp
  m1 <- ic_sp(formula = formula, data = data)
  obs1 <- as.numeric(coef(m1)[trt.name])

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
                  as.formula(paste0("~ . + offset(", trt.name, ".obs * low)")))
      for (i in 1:nperm) {

        # permute based on runit
        data.tmp <- permute(data, trt.name, runit, strat) # permuted data
        model.tmp <- ic_sp(formula = formula.tmp, data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trt.name]) # return tx effect estimate
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
                    as.formula(paste0("~ . + offset(", trt.name, ".obs * up)")))
      for (i in 1:nperm) {

        # permute based on runit
        data.tmp <- permute(data, trt.name, runit, strat) # permuted data
        model.tmp <- ic_sp(formula = formula.tmp, data = data.tmp) # fit
        t <- as.numeric(coef(model.tmp)[trt.name]) # return tx effect estimate
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
