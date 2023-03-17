
#' Methods to compute influence for `lm`, `ivreg`, and `matrix` objects
#'
#' @param ... Options dispatched internally.
#'
#' @return Returns a \code{influence} object.
#'
#' @export
infl <- function(x, ...) {{UseMethod("infl", x)}}

#' @noRd
infl.matrix <- function(x, y, z,
  rm = NULL, options = set_compute(), cluster = NULL) {

  if(missing(z)) {
    influence_lm(list("x" = x, "y" = y),
      rm = rm, options = options, cluster = cluster)
  } else {
    influence_iv(list("x" = list("regressors" = x, "instruments" = z), "y" = y),
      rm = rm, options = options, cluster = cluster)
  }
}

#' @noRd
#' @export
infl.lm <- function(x,
  rm = NULL, options = set_compute(), cluster = NULL) {

  influence_lm(x, rm = rm, options = options, cluster = cluster)
}

#' @noRd
#' @export
infl.ivreg <- function(x,
  rm = NULL, options = set_compute(), cluster = NULL) {

  if(is.null(get_data(x)$Z)) {
    influence_lm(x, rm = rm, options = options, cluster = cluster)
  } else {
    influence_iv(x, rm = rm, options = options, cluster = cluster)
  }
}
