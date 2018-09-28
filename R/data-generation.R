#' Generate data from GLMM
#'
#' \code{gendata_crt} generates a dataset from a generalized linear mixed model
#'
#' This function generates data from a generalized linear mixed model, a GLMM,
#' with a fixed intercept, fixed treatment effect, and random cluster intercept.
#' The data will be sorted increasing by cluster then subject. The GLMM is
#' \deqn{g(E[Y_{ki}]) = \mu + \alpha_k + \theta * X_k}
#' where \eqn{X_k} is the treatment indicator for the \eqn{k}th cluster and
#' \eqn{\alpha ~ N(0, \sigma^2)}.
#'
#' @param family a description of the error distribution and link function to
#'     be used in the data generation model. This can be a character string
#'     naming a family function, a family function or the result of a call to a
#'     family function. (See \code{\link[stats]{family}} for details.)
#'     Currently \code{gendata_crt} supports gaussian, binomial, and poisson.
#' @param nclus number of clusters
#' @param size vector of length 2 specifying range of cluster sizes that will be
#'     drawn from a uniform distribution
#' @param theta treatment effect
#' @param sigma SD of random cluster effect
#' @param mu intercept in GLMM
#' @param rho between-cluster correlation (assuming exchangeable covariance); if
#'     non-zero, then there is corrrelation between clusters.
#' @param sd SD of residual error (only applies if family = gaussian)
#'
#' @examples
#' # generate data from GLMM with bernoulli outcome,
#' # treatment odds ratio of 1.5, baseline odds in
#' # control group of 0.3, 10 clusters with sizes
#' # ranging from 100 to 200.
#'
#' set.seed(444)
#' ds <- gendata_crt(family = binomial, nclus = 10, size = c(100, 200),
#'     theta = log(1.5), sigma = 0.5, mu = log(0.3))
#' head(ds)
#' #   unique.id clusid id trt y
#' # 1       1.1      1  1   1 0
#' # 2       1.2      1  2   1 1
#' # 3       1.3      1  3   1 0
#' # 4       1.4      1  4   1 0
#' # 5       1.5      1  5   1 0
#' # 6       1.6      1  6   1 0
gendata_crt <- function(family = gaussian, nclus, size, theta = 0,
                        sigma, mu, rho = 0, sd = 1) {
  nis <- round(runif(nclus, size[1], size[2]))
  ntot <- sum(nis)
  mymu <- rep(mu, ntot)
  # covariance matrix for random cluster effects
  Sigma <- diag(rep(sigma^2, nclus))
  Sigma[lower.tri(Sigma)] <- Sigma[upper.tri(Sigma)] <- sigma^2 * rho
  # random cluster effects
  alpha <- MASS::mvrnorm(1, rep(0, nclus), Sigma)
  myalpha <- rep(alpha, nis)

  # create an individual id formatted as cluster#.individual#
  id <- unlist(lapply(nis, function(x) 1:x))
  clusid <- rep(1:nclus, nis)
  unique.id <- paste(clusid, id, sep = ".")

  # generate entire vector of treatment assignments
  if (nclus / 2 - round(nclus / 2) > 0) stop("nclus not even")
  trt <- clusid %% 2 # modulo division so that odd clusters trt=1, even trt=0

  # generate outcome
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
        print(family)
        stop("family not recognized")
  }
  y.mean <- family$linkinv(mymu + myalpha + trt * theta)
  if (family$family == "gaussian") {
    y <- rnorm(ntot, y.mean, sd)
  } else if (family$family == "binomial") {
    y <- rbinom(ntot, 1, y.mean)
  } else if (family$family == "poisson") {
    y <- rpois(ntot, y.mean)
  } else {
    print(family)
    stop("family not yet supported by this function")
  }

  out <- data.frame(unique.id, clusid, id, trt, y)
  return(out)
}
