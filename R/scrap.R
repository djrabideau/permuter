#' SEQ version of permci_survreg (for debugging)
#' @export
permci_survreg_SEQ <- function(formula, trtname, runit, strat = NULL, data,
                               dist = "weibull", nperm = 1000, nburn = 0,
                               level = 0.95, init, initmethod = 'perm',
                               seed, quietly = F, ...) {
  set.seed(seed)

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
      perm.stat <- rep(NA, nperm_init)
      for (i in 1:nperm_init) {
        data.tmp <- permute(data, trtname, runit, strat) # permuted data
        model.tmp <- survival::survreg(formula = formula.tmp, data = data.tmp, dist = dist) # fit
        perm.stat[i] <- as.numeric(coef(model.tmp)[trtname])
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
  trace <- matrix(NA, nrow = nperm + nburn, ncol = 2)
  # search for lower
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

    if (!quietly & i %in% seq(ceiling((nperm + nburn) / 10), (nperm + nburn),
                              ceiling((nperm + nburn) / 10)))
      cat(i, "of", (nperm + nburn), "permutations complete\n")
  }
  if (!quietly) cat("lower bound complete\n")
  trace[, 1] <- low.vec

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

    if (!quietly & i %in% seq(ceiling((nperm + nburn) / 10), (nperm + nburn),
                              ceiling((nperm + nburn) / 10)))
      cat(i, "of", (nperm + nburn), "permutations complete\n")
  }
  if (!quietly) cat("upper bound complete\n")
  trace[, 2] <- up.vec

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
                seed = seed
              ))
  class(out) <- 'permci'
  return(out)
}
