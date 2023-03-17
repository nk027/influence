
#' Methods to compute sensitivity with `lm`, `ivreg`, and `influence` objects
#'
#' @param ... Options dispatched internally.
#'
#' @return Returns a \code{sensitivitiy} object.
#'
#' @export
sens <- function(x, ...) {{UseMethod("sens", x)}}

#' @noRd
#' @export
sens.lm <- function(x,
  lambda = set_lambda(), options = set_options(),
  cluster = NULL, verbose = TRUE) {

  sensitivity_lm(x, lambda = lambda, options = options,
    cluster = cluster, verbose = verbose)
}

#' @noRd
#' @export
sens.ivreg <- function(x,
  lambda = set_lambda(), options = set_options(),
  cluster = NULL, verbose = TRUE) {

  if(is.null(get_data(x)$Z)) {
    sensitivity_lm(x, lambda = lambda, options = options,
      cluster = cluster, verbose = verbose)
  } else {
    sensitivity_iv(x, lambda = lambda, options = options,
      cluster = cluster, verbose = verbose)
  }
}

#' @noRd
#' @export
sens.influence <- function(x,
  lambda = set_lambda(), options = set_options(),
  cluster = x$model$cluster, verbose = TRUE) {

  if(x$model$class == "lm") {
    sens.lm(x$model$x, lambda = lambda, options = options,
      cluster = cluster, verbose = verbose)
  } else {
    sens.ivreg(x$model$x, lambda = lambda, options = options,
      cluster = cluster, verbose = verbose)
  }
}
