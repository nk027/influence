
#' Summary method for sensitivity objects
#'
#' @param x A \code{sensitivity} object, obtained from \code{\link{sens}}.
#' @param n Optional character vector with a maximum of removals.
#' @param threshold Optional numerical vector with influence thresholds.
#' @param ... Not used.
#'
#' @return Returns a list with summary results.
#'
#' @export
summary.sensitivity <- function(x, n = 0, threshold = qnorm(.975), ...) {

  n <- min(x$influence$N[1L], if(isTRUE(n <= 0)) Inf else n)

  x <- init.sensitivity(x)
  x$exact <- x$exact[!is.nan(x$exact)][seq(n)]
  x$initial <- x$initial[!is.nan(x$initial)][seq(n)]

  list("exact" = c(
    "zero" = which(diff(sign(x$exact)) != 0) + 0,
    "threshold" = which(diff(abs(x$exact) > threshold) != 0) + 0),
  "initial" = c(
    "zero" = which(diff(sign(x$initial)) != 0) + 0,
    "threshold" = which(diff(abs(x$initial) > threshold) != 0) + 0))
}
