
get_data <- function(x, ...) {{UseMethod("get_data", x)}}

get_data.list <- function(x) {
  y <- x$y
  X <- if(is.null(x$x)) {
    if(is.null(x$X)) {stop("No explanatories found.")} else {as.matrix(x$X)}
  } else if(is.null(x$x$regressors)) {
    as.matrix(x$x)
  } else {as.matrix(x$x$regressors)}
  Z <- if(is.null(x$x$instruments)) NULL else {as.matrix(x$x$instruments)}
  return(list("y" = y, "X" = X, "Z" = Z))
}

get_data.lm <- function(x) {
  if(!is.null(mf <- x$model)) {
    mf <- x$model
    y <- model.response(mf, "numeric")
    X <- model.matrix(attr(mf, "terms"), mf, contrasts = NULL)
  } else if(is.null(y <- x$y) || is.null(X <- x$x)) {
    stop("Please run `lm()` with `model` or `y` and `x` set to `TRUE`.")
  }
  return(list("y" = y, "X" = as.matrix(X)))
}

get_data.ivreg <- function(x) {
  if(!is.null(mf <- x$model)) {
    mf <- x$model
    y <- model.response(mf, "numeric")
    X <- model.matrix(x$terms$regressors, mf, contrasts = NULL)
    Z <- model.matrix(x$terms$instruments, mf, contrasts = NULL)
  } else if(is.null(y <- x$y) ||
    is.null(X <- x$x$regressors || is.null(Z <- x$x$instruments))) {
    stop("Please run `ivreg()` with `model` or `y` and `x` set to `TRUE`.")
  }
  return(list("y" = y, "X" = as.matrix(X), "Z" = as.matrix(Z)))
}


#' Check numeric scalar
#'
#' Check whether an object is bounded and coercible to a numeric value.
#'
#' @param x Numeric scalar.
#' @param min Numeric scalar. Minimum value of \emph{x}.
#' @param max Numeric scalar. Maximum value of \emph{x}.
#' @param fun Function to apply to \emph{x} before returning.
#' @param msg String fed to \code{\link[base]{stop}} if an error occurs.
#'
#' @return Returns \code{fun(x)}.
#'
#' @noRd
num_check <- function(
  x, min = 0, max = Inf,
  msg = "Please check the numeric parameters.",
  fun = as.numeric) {

  if(!is.numeric(x) || length(x) != 1 || x < min || x > max) {stop(msg)}

  return(fun(x))
}

#' @noRd
int_check <- function(
  x, min = 0L, max = Inf,
  msg = "Please check the integer parameters.") {

  num_check(x, min, max, msg, fun = as.integer)
}


#' @noRd
check_cluster <- function(cluster, N) {
  if(!is.null(cluster)) {
    cluster <- as.data.frame(cluster)
    if(NROW(cluster) != N) {stop("Size of 'cluster' does not match the data.")}
    if(anyNA(cluster)) {stop("No missing 'cluster' values are allowed.")}
  }
  return(cluster)
}


#' @noRd
check_iterations <- function(N, n_max, p_max) {
  return(min(N - 1L, n_max, floor(N * p_max)))
}
