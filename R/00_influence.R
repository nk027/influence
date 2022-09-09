
# Provide methods for `lm`, `ivreg`, and `matrix` objects -----

infl <- function(x, ...) {{UseMethod("infl", x)}}


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

infl.lm <- function(x,
  rm = NULL, options = set_compute(), cluster = NULL) {

  influence_lm(x, rm = rm, options = options, cluster = cluster)
}

infl.ivreg <- function(x,
  rm = NULL, options = set_compute(), cluster = NULL) {

  if(is.null(get_data(x)$Z)) {
    influence_lm(x, rm = rm, options = options, cluster = cluster)
  } else {
    influence_iv(x, rm = rm, options = options, cluster = cluster)
  }
}
