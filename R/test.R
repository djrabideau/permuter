# This function carries out a permutation test for a GLM
permtest_glm <- function(formula, trt.name, runit, strat = NULL,
                         family = gaussian, data, nperm = 2000, ncores = 1,
                         quietly = T) {
  # formula:          formula argument passed to glm(), e.g. y ~ x1 + x2
  # trt.name:         character string specifying the name of randomized
  #                     treatment variable in data frame (variable to permute)
  # runit:            character string specifying the name of unit
  #                     of randomization in data frame
  # strat:            character string specifying the name of the variable
  #                     upon which randomization was stratified
  # family:           family argument passed to glm(), e.g. gaussian, binomial
  # data:             data argument passed to glm()
  # ncores:           if >1, then permtest runs in parallel
  # quietly:          T if you want some status updates (only if ncores == 1)

  # fit glm
  m1 <- glm(formula = formula, family = family, data = data)
  obs1 <- as.numeric(coef(m1)[trt.name])

  perm.stat <- rep(0, nperm)

  # permute based on runit
  doMC::registerDoMC(ncores)
  perm.stat <- foreach::foreach(i = 1:nperm, .combine = c) %dopar% {
    data.tmp <- permute(data, trt.name, runit, strat) # permuted data
    model.tmp <- glm(formula = formula, family = family, data = data.tmp) # fit
    if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                            ceiling(nperm / 10)))
      cat(i, "of", nperm, "permutations complete\n")

    as.numeric(coef(model.tmp)[trt.name]) # return tx effect estimate
  }

  comb.stat1 <- c(obs1, perm.stat)

  pval.perm <- mean(abs(comb.stat1) >= abs(obs1)) # prop more extreme than obs

  out <- c(obs1, pval.perm)
  names(out) <- c(trt.name, "p.perm")
  return(out)
}


# Carries out a permutation test for a Cox PH model with interval censoring
permtest_ic_sp <- function(formula, trt.name, runit, strat = NULL, data,
                           nperm = 999, ncores = 1, quietly = F) {
  # formula:          formula argument passed to ic_sp(), e.g. cbind(l, u) ~ trt
  # trt.name:         character string specifying the name of randomized
  #                     treatment variable in data frame (variable to permute)
  # runit:            character string specifying the name of unit
  #                     of randomization in data frame
  # strat:            character string specifying the name of the variable
  #                     upon which randomization was stratified
  # data:             data argument passed to ic_sp()
  # ncores:           if >1, then permtest runs in parallel
  # quietly:          T if you want some status updates (only if ncores == 1)

  # fit ic_sp
  m1 <- ic_sp(formula = formula, data = data)
  obs1 <- as.numeric(coef(m1)[trt.name])

  perm.stat <- rep(0, nperm)

  # permute based on runit
  doMC::registerDoMC(ncores)
  perm.stat <- foreach::foreach(i = 1:nperm, .combine = c) %dopar% {
    data.tmp <- permute(data, trt.name, runit, strat) # permuted data
    model.tmp <- ic_sp(formula = formula, data = data.tmp) # fit

    if (ncores == 1 & !quietly & i %in% seq(ceiling(nperm / 10), nperm,
                                            ceiling(nperm / 10)))
      cat(i, "of", nperm, "permutations complete\n")

    as.numeric(coef(model.tmp)[trt.name]) # return tx effect estimate
  }

  comb.stat1 <- c(obs1, perm.stat)

  pval.perm <- mean(abs(comb.stat1) >= abs(obs1)) # prop more extreme than obs

  out <- c(obs1, pval.perm)
  names(out) <- c(trt.name, "p.perm")
  return(out)
}
