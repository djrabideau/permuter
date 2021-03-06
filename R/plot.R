#' Plot randomization confidence interval
#'
#' This function takes a permci object and produces plots to monitor convergence
#' of each randomization-based confidence interval bound. The search values are
#' indicated by black lines, initial values an x, final update an o, estimated
#' CI values by red horizontal lines.
#'
#' @param x An object of class "permci" returned from one of the permci
#' functions (e.g. \code{\link[permuter]{permci_glm}}).
#' @export
plot.permci <- function(x) {
  xmax <- nrow(x$trace)
  utrace <- c(x$init[2], x$trace[, 2])
  ltrace <- c(x$init[1], x$trace[, 1])
  ulim <- range(utrace)
  llim <- range(ltrace)
  par(mfrow = c(2, 1), mar = c(2, 4, 1, 2) + 0.1, oma = c(2, 0, 0, 0))
  plot(0:xmax, utrace, type = 'n', las = 1, ylab = 'Upper', ylim = ulim)
    abline(h = x$ci[2], col = 'red')
    lines(0:xmax, utrace)
    points(0, x$init[2], pch = 4, lwd = 2)
    points(xmax, x$trace[xmax, 2], lwd = 2)
  plot(0:xmax, ltrace, type = 'n', las = 1, ylab = 'Lower', ylim = llim)
    abline(h = x$ci[1], col = 'red')
    lines(0:xmax, ltrace)
    points(0, x$init[1], pch = 4, lwd = 2)
    points(xmax, x$trace[xmax, 1], lwd = 2)
  mtext('Number of Permutations', 1, line = 2.5)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = rep(0, 4)) # reset par
}

#' Plot randomization distribution
#'
#' This function takes a permtest object and produces a histogram of the Monte
#' Carlo randomization distribution of the test statistic. The observed value is
#' indicated by a red vertical line.
#'
#' @param x An object of class "permtest" returned from one of the permtest
#' functions (e.g. \code{\link[permuter]{permtest_glm}}).
#' @param ... optional arguments to \code{\link[graphics]{hist}}
#' @export
plot.permtest <- function(x, ...) {
  hist(x$permCoefs, xlab = expression(hat(theta)^{(p)}),
       main = 'Randomization Distribution', ...)
  abline(v = x$coef, col = 'red', lwd = 2)
  axis(3, at = x$coef, labels = expression(hat(theta)), col = 'red',
       col.axis = 'red', lwd = 2, line = -1, tick = F)
}
