# file:         gendata_crt
# programmer:   Dustin Rabideau
# description:  This function generates data from a parallel cluster
#               randomized trial (CRT) using a GLMM.
#               The data will be sorted increasing by cluster then subject.
#
# Simulates data from the following GLMM:
#   g(E[Y_ki]) = mu + alpha_k + theta * trt_k

library(MASS)

gendata.crt <- function(family = gaussian, nclus, size, theta = 0,
                        sigma, mu, rho = 0, sd = 1) {
  # family: outcome type (currently supports gaussian, binomial, poisson)
  # nclus:  number of clusters
  # size:   vector of length 2 specifying range of cluster sizes
  #           (to be drawn from uniform dist)
  # theta:  treatment effect
  # sigma:  SD of random cluster effect
  # mu:     intercept (i.e. mean outcome with covariates 0)
  # rho:    between-cluster correlation (assuming exchangeable covariance); if
  #         non-zero, then there is correlation between clusters.
  # sd:  SD of residual error (only applies if family=gaussian)

  nis <- round(runif(nclus, size[1], size[2]))
  ntot <- sum(nis)
  mymu <- rep(mu, ntot)
  # covariance matrix for random cluster effects
  Sigma <- diag(rep(sigma^2, nclus))
  Sigma[lower.tri(Sigma)] <- Sigma[upper.tri(Sigma)] <- sigma^2 * rho
  # random cluster effects
  alpha <- mvrnorm(1, rep(0, nclus), Sigma)
  myalpha <- rep(alpha, nis)

  # create an individual id formatted as cluster#.individual#
  id <- unlist(lapply(nis, function(x) 1:x))
  clusid <- rep(1:nclus, nis)
  unique.id <- paste(clusid, id, sep='.')

  # generate entire vector of treatment assignments
  if (nclus / 2 - round(nclus / 2) > 0) stop('nclus not even')
  trt <- clusid %% 2 # modulo division so that odd clusters trt=1, even trt=0

  # generate outcome
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
  }
  y.mean <- family$linkinv(mymu + myalpha + trt * theta)
  if (family$family == 'gaussian') {
    y <- rnorm(ntot, y.mean, sd)
  } else if (family$family == 'binomial') {
    y <- rbinom(ntot, 1, y.mean)
  } else if (family$family == 'poisson') {
    y <- rpois(ntot, y.mean)
  } else {
    print(family)
    stop("'family' not yet supported by this function")
  }

  out <- data.frame(unique.id, clusid, id, trt, y)
  return(out)
}
