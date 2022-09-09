
# Provide methods for `lm`, `ivreg`, and `influence` objects -----

sens <- function(x, ...) {{UseMethod("sens", x)}}


sens.lm <- function(x,
  lambda = set_lambda(), options = set_options(),
  cluster = NULL, verbose = TRUE) {

  sensitivity_lm(x, lambda = lambda, options = options,
    cluster = cluster, verbose = verbose)
}

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
