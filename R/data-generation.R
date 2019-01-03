#' Generate data from GLMM
#'
#' \code{gendata_crt} generates a dataset from a generalized linear mixed model
#'
#' This function generates data from a generalized linear mixed model, a GLMM,
#' with a fixed intercept, fixed treatment effect, and random cluster intercept.
#' The data will be sorted increasing by cluster then subject. The GLMM is
#' \deqn{g(E[Y_{ki}]) = \mu + \alpha_k + \theta * X_k}
#' where \eqn{X_k} is the treatment indicator for the \eqn{k}th cluster and
#' \eqn{\alpha ~ N(0, Sigma)}. If only \code{sigma} is specified, then
#' \code{Sigma} is exchangeable with \code{sigma^2} on the diagonals and
#' \code{rho * sigma^2} on the off-diagonals.
#'
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
#' @param Sigma optional user-specified covariance matrix for random cluster effects
#' @param mu intercept in GLMM
#' @param rho between-cluster correlation (assuming exchangeable covariance); if
#'     non-zero, then there is corrrelation between clusters.
#' @param sd SD of residual error (only applies if family = gaussian)
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
                        sigma, Sigma, mu, rho = 0, sd = 1) {
  if (length(nclus) != 2)
    stop("length(nclus) should be 2")
  if (length(size) != 2)
    stop("length(size) should be 2")

  nclus_tot <- sum(nclus)
  nis <- round(runif(nclus_tot, size[1], size[2]))
  ntot <- sum(nis)
  mymu <- rep(mu, ntot)
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
