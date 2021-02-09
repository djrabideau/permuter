#' Generate data from a parallel cluster randomized trial (CRT)
#'
#' \code{gendata_crt} generates a dataset from a parallel CRT based on a GLMM
#'
#' This function generates data from a generalized linear mixed model, a GLMM,
#' with a fixed intercept, fixed treatment effect, and random cluster intercept.
#' The data will be sorted increasing by cluster then subject. The GLMM is
#' \deqn{g(E[Y_{ki}]) = \mu + \alpha_k + \theta * X_k}
#' where \eqn{X_k} is the treatment indicator for the \eqn{k}th cluster and
#' \eqn{\alpha ~ N(0, Sigma)}. If only \code{sigma} is specified, then
#' \code{Sigma} is exchangeable with \code{sigma^2} on the diagonals and
#' \code{rho * sigma^2} on the off-diagonals. If redist = 'lognormal', then
#' \eqn{\alpha ~ logN(meanlog = 0, sdlog = sigma)}.
#' @param family a description of the error distribution and link function to
#'     be used in the data generation model. This can be a character string
#'     naming a family function, a family function or the result of a call to a
#'     family function. (See \code{\link[stats]{family}} for details.)
#'     Currently \code{gendata_crt} supports gaussian, binomial, and poisson.
#' @param nclus vector of length 2 specifying the number of clusters in
#' treatment group (trt == 1) and control group (trt == 0), respectively
#' @param size vector of length 2 specifying range of cluster sizes that will be
#'     drawn from a uniform distribution
#' @param theta treatment effect
#' @param sigma SD of random cluster effect; ignored if \code{Sigma} is specified
#' @param Sigma optional user-specified covariance matrix for multivariate
#'     normal random cluster effects
#' @param mu intercept in GLMM
#' @param rho between-cluster correlation (assuming exchangeable covariance); if
#'     non-zero, then there is corrrelation between clusters.
#' @param sd SD of residual error (only applies if family = gaussian)
#' @param redist assumed distribution for random effects (currently supports
#'     'normal' and 'lognormal')
#'
#' @examples
#' # generate data from GLMM with bernoulli outcome,
#' # treatment odds ratio of 1.5, baseline odds in
#' # control group of 0.3, 10 clusters (5 in each group)
#' # with sizes ranging from 100 to 200.
#'
#' set.seed(444)
#' ds <- gendata_crt(family = binomial, nclus = c(5, 5), size = c(100, 200),
#'     theta = log(1.5), sigma = 0.5, mu = log(0.3))
#' head(ds)
#' #   unique.id clusid id trt y
#' # 1       1.1      1  1   1 0
#' # 2       1.2      1  2   1 1
#' # 3       1.3      1  3   1 0
#' # 4       1.4      1  4   1 0
#' # 5       1.5      1  5   1 0
#' # 6       1.6      1  6   1 0
#' @export
gendata_crt <- function(family = gaussian, nclus, size, theta = 0,
                        sigma, Sigma, mu, rho = 0, sd = 1, redist = 'normal') {
  if (length(nclus) != 2)
    stop("length(nclus) should be 2")
  if (length(size) != 2)
    stop("length(size) should be 2")

  nclus_tot <- sum(nclus)
  nis <- round(runif(nclus_tot, size[1], size[2]))
  ntot <- sum(nis)
  mymu <- rep(mu, ntot)

  # random effects
  if (redist == 'normal') {
    # covariance matrix for random cluster effects
    if (missing(Sigma)) {
      Sigma <- diag(rep(sigma^2, nclus_tot))
      Sigma[lower.tri(Sigma)] <- Sigma[upper.tri(Sigma)] <- sigma^2 * rho
    } else {
      if (!identical(as.numeric(dim(Sigma)), c(nclus_tot, nclus_tot)))
        stop(paste0("dim(Sigma) should be c(", nclus_tot, ', ', nclus_tot, ')'))
    }
    # random cluster effects
    alpha <- MASS::mvrnorm(1, rep(0, nclus_tot), Sigma)
  } else if (redist == 'lognormal') {
    alpha <- rlnorm(nclus_tot, 0, sigma)
  } else {
    stop(paste0('redist = ', redist, ' not supported'))
  }
  myalpha <- rep(alpha, nis)

  # create an individual id formatted as cluster#.individual#
  id <- unlist(lapply(nis, function(x) 1:x))
  clusid <- rep(1:nclus_tot, nis)
  unique.id <- paste(clusid, id, sep = ".")

  # generate entire vector of treatment assignments where first nclus[1] are
  # trt = 1, next nclus[2] are trt = 0
  trt <- as.numeric(clusid <= nclus[1])

  # generate outcome
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
        stop(paste0("family '", family, "' not recognized"))
  }
  y.mean <- family$linkinv(mymu + myalpha + trt * theta)
  if (family$family == "gaussian") {
    y <- rnorm(ntot, y.mean, sd)
  } else if (family$family == "binomial") {
    y <- rbinom(ntot, 1, y.mean)
  } else if (family$family == "poisson") {
    y <- rpois(ntot, y.mean)
  } else {
    stop(paste0("family '", family$family, "' not supported by gendata_crt"))
  }

  out <- data.frame(unique.id, clusid, id, trt, y)
  return(out)
}

#' Generate data from a stepped wedge cluster randomized trial (SW-CRT)
#'
#' \code{gendata_swcrt} generates a dataset from a SW-CRT based on a GLMM
#'
#' This function generates data from a SW-CRT using a generalized linear mixed
#' model, a GLMM, with a fixed intercept, fixed treatment effect, fixed period
#' effect, random cluster effects, random cluster-period effects, and random
#' cluster-treatment effects.
#' The data will be sorted increasing by cluster, then period, then subject.The
#' GLMM is \deqn{g(E[Y_{kij}]) = \mu + \alpha_k + \beta_j + \eta_{jk} + \theta * X_{jk} + \omega_k * X_{jk}}
#' where \eqn{X_{kj}} is the treatment indicator for the \eqn{k}th cluster in
#' the \eqn{j}th time period (step).
#' @param family a description of the error distribution and link function to
#'     be used in the data generation model. This can be a character string
#'     naming a family function, a family function or the result of a call to a
#'     family function. (See \code{\link[stats]{family}} for details.)
#'     Currently \code{gendata_swcrt} supports gaussian, binomial, and poisson.
#' @param nclus total number of clusters
#' @param size if numeric, number of individuals per cluster-period; if vector
#'     (must be length 2), specifies the range of cluster-period sizes, each of
#'     which will be drawn independently from a uniform distribution; if matrix
#'     and ncol(size) == 2, the kth row specifies the range of
#'     cluster-period sizes for the kth cluster, similarly drawn from a uniform
#'     distribution; if matrix and ncol(size) == (nstep + 1) and sizeFixed == T,
#'     elements represent fixed cluster-period sizes.
#' @param sizeFixed if T, size matrix represents fixed cluster-period sizes
#'     (i.e. not drawn uniformly based on range of sizes, just fixed at provided
#'     sizes).
#' @param nstep number of randomization periods not including baseline, i.e.
#'     (nstep + 1) = number of periods.
#' @param theta treatment effect
#' @param sigma SD of random cluster effect; ignored if \code{Sigma} is specified
#' @param Sigma optional user-specified covariance matrix for multivariate
#'     normal random cluster effects
#' @param mu intercept in GLMM
#' @param beta vector of fixed period effects (length = nstep + 1)
#' @param nu SD of random cluster-period effects
#' @param phi SD of random cluster-treatment effects
#' @param rho between-cluster correlation (assuming exchangeable covariance); if
#'     non-zero, then there is correlation between clusters.
#' @param sd SD of residual error (only applies if family = gaussian)
#' @param redist assumed distribution for random cluster effects (currently supports
#'     'normal' and 'lognormal')
#'
#' @examples
#' # generate data from SW-CRT with continuous outcome,
#' # treatment effect of 1.5, baseline mean outcome in
#' # control group of 0, 10 clusters, 5 steps, period effects of 0:5,
#' # with sizes ranging from 20 to 40.
#'
#' set.seed(444)
#' ds <- gendata_swcrt(nclus = 10, size = c(20, 40), nstep = 5, theta = 1.5,
#'                     sigma = 0.5, mu = 0, beta = 0:5, nu = 0.2, sd = 1)
#' head(ds)
#' #   cluster period individual cluster.period.individual cluster.period trt          y
#' # 1       1      1          1                     1.1.1            1.1   0 -0.5186888
#' # 2       1      1          2                     1.1.2            1.1   0 -0.3058491
#' # 3       1      1          3                     1.1.3            1.1   0 -1.0040078
#' # 4       1      1          4                     1.1.4            1.1   0 -0.6668997
#' # 5       1      1          5                     1.1.5            1.1   0 -0.5208075
#' # 6       1      1          6                     1.1.6            1.1   0 -1.7543536
#' @export
gendata_swcrt <- function(family = gaussian, nclus, size, sizeFixed = F, nstep, theta = 0, sigma, Sigma, mu,
                     beta, nu = 0, phi = 0, rho = 0, sd = 1, redist = 'normal') {
  # check inputs
  if (nclus %% nstep != 0)
    stop('nclus must be a multiple of nstep')
  if (length(beta) != nstep + 1)
    stop('length(beta) must be equal to nstep + 1')

  # setup size matrix
  if (class(size) == 'numeric') {
    if (length(size) > 2) {
      stop('length(size) != 2')
    } else if (length(size) == 1) {
        size <- c(size, size)
    }
    size <- matrix(rep(size, nclus), nrow = nclus, byrow = T)
  } else if (class(size) == 'matrix') {
    if ((sum(dim(size) - c(nclus, 2) != c(0, 0)) > 1) & sizeFixed == F)
      stop('nrow(size) != nclus')
    if ((sum(dim(size) - c(nclus, nstep + 1) != c(0, 0)) > 1) & sizeFixed == T)
      stop('nrow(size) != nclus or ncol(size) != (nstep + 1)')
  }

  nperiod <- (nstep + 1)

  if (sizeFixed) {
    nis <- split(size, rep(1:nrow(size), each = ncol(size)))
  } else {
    nis <- list()
    for (k in 1:nclus) {
      nis[[k]] <- round(runif(nperiod, size[k, 1], size[k, 2]))
    }
  }

  ntotk <- unlist(lapply(nis, sum))
  ntot <- sum(ntotk)

  mymu <- rep(mu, ntot)

  # random cluster effects
  if (redist == 'normal') {
    # covariance matrix for random cluster effects
    if (missing(Sigma)) {
      Sigma <- diag(rep(sigma^2, nclus))
      Sigma[lower.tri(Sigma)] <- Sigma[upper.tri(Sigma)] <- sigma^2 * rho
    } else {
      if (!identical(as.numeric(dim(Sigma)), c(nclus, nclus)))
        stop(paste0("dim(Sigma) should be c(", nclus, ', ', nclus, ')'))
    }
    # random cluster effects
    alpha <- MASS::mvrnorm(1, rep(0, nclus), Sigma)
  } else if (redist == 'lognormal') {
    alpha <- rlnorm(nclus, 0, sigma)
  } else {
    stop(paste0('redist = ', redist, ' not supported'))
  }
  myalpha <- rep(alpha, ntotk)

  # fixed period effects
  mybeta <- unlist(lapply(nis, function(x) rep(beta, x)))

  # random cluster-period effect
  eta <- rnorm(nclus * nperiod, 0, nu)
  myeta <- rep(eta, unlist(nis))

  # random cluster-treatment effects
  omega <- rnorm(nclus, 0, phi)
  myomega <- rep(omega, ntotk)

  # create an individual id formatted as cluster#.period#.individual#
  cluster <- rep(1:nclus, ntotk)
  period <- unlist(lapply(nis, function(x) rep(1:nperiod, x)))
  individual <- unlist(lapply(unlist(nis), function(x) 1:x))
  cluster.period.individual <- paste(cluster, period, individual, sep='.')
  cluster.period <- paste(cluster, period, sep='.')

  # generate sw treatment matrix
  trt.mat <- matrix(1, nrow = nstep, ncol = nperiod)
  trt.mat[lower.tri(trt.mat, diag = T)] <- 0

  # generate entire vector of treatment assignments
  nperstep <- nclus / nstep
  if (nperstep - round(nperstep) > 0) stop('nclus not a multiple of nstep')
  trt <- c()
  for (row in 1:nrow(trt.mat)) {
    for (r in 1:nperstep) {
      tmp <- rep(trt.mat[row, ], nis[[(row - 1) * nperstep + r]])
      trt <- c(trt, tmp)
    }
  }

  # generate outcome
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
  }
  y.mean <- family$linkinv(mymu + myalpha + mybeta + myeta + trt * theta + trt * myomega)
  if (family$family=='gaussian') {
    y <- rnorm(ntot, y.mean, sd)
  } else if (family$family=='binomial') {
    y <- rbinom(ntot, 1, y.mean)
  } else if (family$family=='poisson') {
    y <- rpois(ntot, y.mean)
  } else {
    print(family)
    stop("'family' not yet supported by this function")
  }

  out <- data.frame(cluster, period, individual, cluster.period.individual, cluster.period, trt, y)
  return(out)
}

#' @export
gendata_swcrtStrat <- function(family = gaussian, nclus, size, nstep, theta = 0, sigma, Sigma, mu,
                               beta, nu = 0, phi = 0, rho = 0, sd = 1, redist = 'normal',
                               gamma = 0) {
  # check inputs
  if (nclus %% nstep != 0)
    stop('nclus must be a multiple of nstep')
  if (length(beta) != nstep + 1)
    stop('length(beta) must be equal to nstep + 1')
  if (nclus == nstep)
    stop('nclus must be a multiple (>=2) for stratified SW-CRT')

  # setup size matrix
  if (class(size) == 'numeric') {
    if (length(size) > 2) {
      stop('length(size) != 2')
    } else if (length(size) == 1) {
      size <- c(size, size)
    }
    size <- matrix(rep(size, nclus), nrow = nclus, byrow = T)
  } else if (class(size) == 'matrix') {
    if (sum(dim(size) - c(nclus, 2) != c(0, 0)) > 1)
      stop('nrow(size) != nclus')
  }

  nperiod <- (nstep + 1)

  nis <- list()
  for (k in 1:nclus) {
    nis[[k]] <- round(runif(nperiod, size[k, 1], size[k, 2]))
  }

  ntotk <- unlist(lapply(nis, sum))
  ntot <- sum(ntotk)

  mymu <- rep(mu, ntot)

  # random cluster effects
  if (redist == 'normal') {
    # covariance matrix for random cluster effects
    if (missing(Sigma)) {
      Sigma <- diag(rep(sigma^2, nclus))
      Sigma[lower.tri(Sigma)] <- Sigma[upper.tri(Sigma)] <- sigma^2 * rho
    } else {
      if (!identical(as.numeric(dim(Sigma)), c(nclus, nclus)))
        stop(paste0("dim(Sigma) should be c(", nclus, ', ', nclus, ')'))
    }
    # random cluster effects
    alpha <- MASS::mvrnorm(1, rep(0, nclus), Sigma)
  } else if (redist == 'lognormal') {
    alpha <- rlnorm(nclus, 0, sigma)
  } else {
    stop(paste0('redist = ', redist, ' not supported'))
  }
  myalpha <- rep(alpha, ntotk)

  # fixed period effects
  mybeta <- unlist(lapply(nis, function(x) rep(beta, x)))

  # random cluster-period effect
  eta <- rnorm(nclus * nperiod, 0, nu)
  myeta <- rep(eta, unlist(nis))

  # random cluster-treatment effects
  omega <- rnorm(nclus, 0, phi)
  myomega <- rep(omega, ntotk)

  # create an individual id formatted as cluster#.period#.individual#
  cluster <- rep(1:nclus, ntotk)
  period <- unlist(lapply(nis, function(x) rep(1:nperiod, x)))
  individual <- unlist(lapply(unlist(nis), function(x) 1:x))
  cluster.period.individual <- paste(cluster, period, individual, sep='.')
  cluster.period <- paste(cluster, period, sep='.')

  # generate sw treatment matrix
  trt.mat <- matrix(1, nrow = nstep, ncol = nperiod)
  trt.mat[lower.tri(trt.mat, diag = T)] <- 0

  # generate entire vector of treatment assignments
  nperstep <- nclus / nstep
  if (nperstep - round(nperstep) > 0) stop('nclus not a multiple of nstep')
  trt <- c()
  for (row in 1:nrow(trt.mat)) {
    for (r in 1:nperstep) {
      tmp <- rep(trt.mat[row, ], nis[[(row - 1) * nperstep + r]])
      trt <- c(trt, tmp)
    }
  }

  # generate stratification variable (binary for now)
  strat <- rep(rep(0:1, nclus / 2), ntotk)

  # generate outcome
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  y.mean <- family$linkinv(mymu + myalpha + mybeta + myeta + trt * theta + trt * myomega + strat * gamma)
  if (family$family=='gaussian') {
    y <- rnorm(ntot, y.mean, sd)
  } else if (family$family=='binomial') {
    y <- rbinom(ntot, 1, y.mean)
  } else if (family$family=='poisson') {
    y <- rpois(ntot, y.mean)
  } else {
    print(family)
    stop("'family' not yet supported by this function")
  }

  out <- data.frame(cluster, period, individual, cluster.period.individual, cluster.period, trt, strat, y)
  return(out)
}
