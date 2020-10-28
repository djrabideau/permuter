#' Permute randomized treatment assignments
#'
#' \code{permute} returns data with treatment assigments permuted
#'
#' This function inputs a data frame and outputs the same data frame with the
#' \code{trtname} column permuted according to the randomization scheme
#' provided. This is a general utility function used within test and confidence
#' interval functions.
#'
#' For example, in a parallel cluster randomized trial (CRT), the unit of
#' randomization is the cluster. Cluster is also the unit corresponding to each
#' unique scalar treatment assigment value (i.e. every observation within a
#' cluster is assigned the same treatment). In this case, \code{runit =
#' trtunit = 'cluster'}.
#'
#' In a stepped wedge CRT, the unit of randomization is also the cluster, but
#' each cluster is randomized to crossover from control to treatment at a
#' particular time period (i.e. each cluster is randomized to a sequence of
#' treatment assigments). In this case, \code{runit = 'cluster'} and
#' \code{trtunit = c('cluster', 'period')}. This will tell \code{permute()} to
#' permute the entire treatment sequence at a cluster level.
#'
#' @param data a data frame
#' @param trtname character string specifying the name of randomized treatment
#' variable in \code{data} (i.e. specifies the variable to permute)
#' @param runit character string specifying the name of unit of randomization
#' in \code{data} (i.e. specifies the level at which trtname will be permuted)
#' @param trtunit character string specifying the name of the unit that
#' corresponds to a unique scalar treatment assignment value (i.e. specifies how
#' trtname will be permuted at runit level, see Details below). If trtunit =
#' NULL (default), trtunit is assumed equivalent to runit.
#' @param strat an optional character string specifying the name of the variable
#' in \code{data} upon which randomization was stratified
#' @param restrictedDf an optional matrix of restricted randomizations from
#' which to choose for randomization test/CI. Rows correspond to randomization
#' units (e.g. individuals, clusters), first column contains randomization unit
#' labels corresponding with data, other columns correspond to different possible
#' randomization patterns. Entries depend on "restrictedDfType", see below.
#' @param restrictedDfType type of restricted randomization matrix provided.
#' If "parallel", entries indicate control (0) or intervention (1). If "SWT",
#' entries indicate period at which clusters crossover.
#' @export
permute <- function(data, trtname, runit, trtunit = NULL, strat = NULL,
                    restrictedDf = NULL, restrictedDfType = NULL) {
  if (is.null(trtunit))
    trtunit <- runit
  if (!(runit %in% trtunit))
    stop('runit must be a subset of trtunit')
  if (is.null(strat) & is.null(restrictedDf)) {
    # permute
    design <- unique(data[, c(trtunit, trtname)])
    if (nrow(as.matrix(unique(data[, c(trtunit)]))) != nrow(design))
      stop('each unique trtunit must correspond to only one value of trtname')
    pdesign <- design
    ldesign <- split(pdesign, pdesign[, runit])
    shuffled_trt_list <- sample(lapply(ldesign, function(x) x[, trtname]))
    pdesign[, trtname] <- unlist(shuffled_trt_list)
  } else if (is.null(restrictedDf)) {
    # stratified permute
    design <- unique(data[, c(trtunit, strat, trtname)])
    if (nrow(as.matrix(unique(data[, c(trtunit)]))) != nrow(design))
      stop('each unique trtunit must correspond to only one value of trtname')
    pdesign <- design
    for (s in unique(design[, strat])) {
      pdesign_tmp <- pdesign[pdesign[, strat] == s, ]
      ldesign_tmp <- split(pdesign_tmp, pdesign_tmp[, runit])
      shuffled_trt_list_tmp <- sample(lapply(ldesign_tmp, function(x) x[, trtname]))
      pdesign[pdesign[, strat] == s, trtname] <- unlist(shuffled_trt_list_tmp)
    }
    pdesign <- pdesign[, - which(names(pdesign) == strat)]
  } else if (restrictedDfType == 'parallel') {
    stop('restrictedDfType not implemented yet')
  } else if (restrictedDfType == 'SWT') {
    # choose one of the restricted randomizations
    dimRestr <- dim(restrictedDf)
    if (dimRestr[1] != length(unique(ds[, runit])))
      stop('ncol(restrictedDf) != length(unique(runit))')
    # periodsDf <- unique(as.vector(as.matrix(restrictedDf[, -1])))
    # periodsDs <- unique(data[, trtunit[2]])
    # if (sum(periodsDs %in% periodsDf) < length(periodsDs))
    #   stop('restrictedDf has entries that do not align with data runit values')
    col <- sample(2:ncol(restrictedDf), 1)
    rand <- restrictedDf[, c(1, col)]
    names(rand) <- c(runit, 'crossoverPeriod')

    # design df
    design <- unique(data[, c(trtunit, trtname)])
    if (nrow(as.matrix(unique(data[, c(trtunit)]))) != nrow(design))
      stop('each unique trtunit must correspond to only one value of trtname')
    design2 <- merge(design, rand, all.x = T)

    # replace actual trt with chosen restricted values
    pdesign <- design
    pdesign[, trtname] <- as.numeric(design2[, trtunit[2]] >= design2[, 'crossoverPeriod'])
  } else {
    stop('restrictedDfType not accepted')
  }

  data <- merge(data[, -which(names(data) == trtname)], pdesign, by = trtunit)
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

#' Get initial values for CI search
#'
#' \code{getinits} returns initial values for CI search, \code{c(lower, upper)}
#'
#' If \code{initmethod = "asymp"}, initial bounds are based on asymptotic
#' approximation (e.g. Wald CI for GLM). If \code{initmethod = "perm"}, initial
#' bounds are based on the permutation approach used in Garthwaite (1996) with
#' \eqn{\hat{\theta}\pm \{(t_2 - t_1)/2\}}, where \eqn{t_1} and \eqn{t_2} denote
#' the second smallest and second largest estimates from the permutation test.
#'
#' This is a general utility function used within \code{permci}.
#'
getInits <- function(model, trtname, runit, trtunit, strat, data, initmethod, alpha, obs1, restrictedDf, restrictedDfType) {
  if (initmethod == 'asymp') {
    if (sum(c('glm', 'survreg', 'coxph') %in% class(model)) > 0) {
      Vcov <- vcov(model, useScale = FALSE)
      trt.se <- sqrt(Vcov[trtname, trtname])
      lower <- obs1 - qnorm(1 - alpha / 2) * trt.se
      upper <- obs1 + qnorm(1 - alpha / 2) * trt.se
      return(c(lower, upper))
    } else {
      stop(paste0("initmethod = '", initmethod,
                  "' not available for class(model)='", paste(class(model), collapse = ', '),
                  "'. Please provide starting values with 'init'"))
    }

  } else if (initmethod == 'perm') {
    if (sum(c('glm', 'survreg', 'coxph') %in% class(model)) > 0) {
      # initialize using quick randomization test of H0: theta = obs1,
      # as recommended in Garthwaite (1996)
      nperm_init <- ceiling((4 - alpha) / alpha)
      data$obs1 <- obs1
      perm.stat <- foreach::foreach(i = 1:nperm_init, .combine = c) %dorng% {
        data.tmp <- permute(data, trtname, runit, trtunit, strat,
                            restrictedDf, restrictedDfType) # permuted data
        model.tmp <- update(model, formula. = as.formula(paste0("~ . + offset(", trtname, ".obs * obs1)")), data = data.tmp) # fit
        as.numeric(coef(model.tmp)[trtname]) # return tx effect estimate
      }
      t1 <- sort(perm.stat)[2] # 2nd to smallest
      t2 <- sort(perm.stat)[nperm_init - 1] # 2nd to largest
      lower <- obs1 - ((t2 - t1) / 2)
      upper <- obs1 + ((t2 - t1) / 2)
      return(c(lower, upper))
    } else {
      stop(paste0("initmethod = '", initmethod,
                  "' not available for class(model)='", paste(class(model), collapse = ', '),
                  "'. Please provide starting values with 'init'"))
    }
  } else {
      print(initmethod)
      stop("'initmethod' not recognized")
  }
}

#' Get initial values for CI search (continuous outcome)
#'
#' \code{getinits} returns initial values for CI search, \code{c(lower, upper)}
#'
#' Initial bounds are based on the permutation approach used in Garthwaite (1996) with
#' \eqn{\hat{\theta}\pm \{(t_2 - t_1)/2\}}, where \eqn{t_1} and \eqn{t_2} denote
#' the second smallest and second largest estimates from the permutation test.
#'
#' This is a general utility function used within \code{permci_cont}.
#'
getInits_cont <- function(f, trtname, runit, trtunit, strat, data, alpha, obs1, restrictedDf, restrictedDfType) {
    # initialize using quick randomization test of H0: theta = obs1,
    # as recommended in Garthwaite (1996)
    nperm_init <- ceiling((4 - alpha) / alpha)
    data$obs1 <- obs1
    perm.stat <- foreach::foreach(i = 1:nperm_init, .combine = c) %dorng% {
      data.tmp <- permute(data, trtname, runit, trtunit, strat,
                          restrictedDf, restrictedDfType) # permuted data
      f(data.tmp)
    }
    t1 <- sort(perm.stat)[2] # 2nd to smallest
    t2 <- sort(perm.stat)[nperm_init - 1] # 2nd to largest
    lower <- obs1 - ((t2 - t1) / 2)
    upper <- obs1 + ((t2 - t1) / 2)
    return(c(lower, upper))
}
