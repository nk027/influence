
lambda <- function(type, variables, combination,
  f = function(x) {x$influence$beta_i[, 1L]}) {

  # Allow choosing predefined functions

  # Check custom functions

  return(f)
}