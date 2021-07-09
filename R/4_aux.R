
rank_influence <- function(x, lambda) {
  value <- lambda(x)
  order <- order(value, decreasing = FALSE, method = "radix")
  cbind("value" = value, "order" = order)
}

mdl_to_mat <- function(x) {
  if(inherits(x, "lm")) {
    return(lm_to_mat(x))
  } else if(inherits(x, "ivreg")) {
    if(is.null(x$terms$instruments)) {
      return(lm_to_mat(x))
    } else {
      return(iv_to_mat(x))
    }
  } else {
    stop("No method found to obtain data from the model object.")
  }
}

# Get matrices from `lm()` fit
lm_to_mat <- function(x) {
  if(!is.null(mf <- x$model)) {
    mf <- x$model
    y <- model.response(mf, "numeric")
    X <- model.matrix(attr(mf, "terms"), mf, contrasts = NULL)
  } else if(is.null(y <- x$y) || is.null(X <- x$x)) {
    stop("Please run `lm()` with `model` or `y` and `x` set to `TRUE`.")
  }
  return(list("y" = y, "X" = as.matrix(X)))
}

# Get matrices from `ivreg()` fit
iv_to_mat <- function(x) {
  if(!is.null(mf <- x$model)) {
    mf <- x$model
    y <- model.response(mf, "numeric")
    X <- model.matrix(x$terms$regressors, mf, contrasts = NULL)
    Z <- model.matrix(x$terms$instruments, mf, contrasts = NULL)
  } else if(is.null(y <- x$y) || is.null(X <- x$x$regressors ||
    is.null(Z <- x$x$instruments))) {
    stop("Please run `ivreg()` with `model` or `y` and `x` set to `TRUE`.")
  }
  return(list("y" = y, "X" = as.matrix(X), "Z" = as.matrix(Z)))
}

# Update inverse using the Sherman-Morrison formula
sherman_morrison <- function(XXi, X_rm) {
  XXi + (XXi %*% tcrossprod(X_rm) %*% XXi) /
    as.numeric(1 - crossprod(X_rm, XXi) %*% X_rm)
}

update_inv <- function(XX_inv, X_rm) {
  XX_inv + (XX_inv %*% crossprod(X_rm) %*% XX_inv) /
    as.numeric(1 - X_rm %*% tcrossprod(XX_inv, X_rm))
}

# Update crossproduct
update_cp <- function(XY, X_rm, Y_rm = X_rm) {
  XY - crossprod(X_rm, Y_rm)
}

# Solve using an upper triangular matrix
solve_cholesky <- function(R, b) {
  backsolve(R, forwardsolve(R, b, upper.tri = TRUE, transpose = TRUE))
}

# Update data using Frisch-Waugh-Lovell theorem
frisch_waugh_lovell <- function(X, y, variables) {
  if(!any(variables == 0)) {
    Q_fwl <- qr.Q(qr(X[, -variables, drop = FALSE]))
    y <- y - Q_fwl %*% crossprod(Q_fwl, y)
    X <- X[, variables, drop = FALSE] - Q_fwl %*%
      crossprod(Q_fwl, X[, variables, drop = FALSE])
  }
  return(list("y" = y, "X" = X))
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
