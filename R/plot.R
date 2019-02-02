#' Plot randomization confidence interval
#'
#' This function takes a permci object and produces plots to monitor convergence
#' of each randomization-based confidence interval bound.
#'
#' @param x An object of class "permci" returned from one of the permci functions
#' (e.g. \code{\link[permuter]{permci_glm}}).
#' @param ... optional arguments to \code{\link[graphics]{plot}}
#' @export
plot.permci <- function(x) {
  xmax <- nrow(x$trace)
  utrace <- c(x$init[2], x$trace[, 2])
  ltrace <- c(x$init[1], x$trace[, 1])
  ulim <- range(utrace)
  llim <- range(ltrace)
  par(mfrow = c(2, 1), mar = c(2, 4, 1, 2) + 0.1, oma = c(2, 0, 0, 0))
  plot(0:xmax, utrace, type = 'l', las = 1, ylab = 'Upper', ylim = ulim)
    points(0, x$init[2], pch = 4, lwd = 2)
    points(xmax, x$trace[xmax, 2], lwd = 2)
  plot(0:xmax, ltrace, type = 'l', las = 1, ylab = 'Lower', ylim = llim)
    points(0, x$init[1], pch = 4, lwd = 2)
    points(xmax, x$trace[xmax, 1], lwd = 2)
  mtext('Number of Permutations', 1, line = 2.5)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = rep(0, 4)) # reset par
}
