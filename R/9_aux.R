
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


#' Create object for `compute_initial()`
#'
#' @noRd
create_object <- function(x, rank, lambda) {

  is_lm <- isTRUE(x$meta$class == "lm")
  K <- length(x$model$beta)

  out <- list(
    "model" = as.data.frame(matrix(NA_real_,
      1L, 2 + length(x$model$beta) * 2 + if(is_lm) {3} else {4},
      dimnames = list(NULL, c("N", "sigma",
        paste0("beta_", seq.int(K)), paste0("se_", seq.int(K)),
        if(is_lm) {c("R2", "F", "LL")} else {c("R2", "F", "R2_1st", "F_1st")}))
    )),
    "initial" = data.frame(
      "id" = rank[, "order"], "lambda" = rank[rank[, "order"], "value"]
    ),
    "meta" = list("lambda" = lambda)
  )
  out$model[1, ] <- c(NROW(x$hat), x$model$sigma, x$model$beta, x$model$se,
    if(is_lm) {
      c(x$model$r2, x$model$fstat, x$model$ll)
    } else {
      c(x$model$r2, x$model$fstat, x$model$r2_first, x$model$fstat_first)
    })

  return(out)
}


#' Get id to retrieve lambda in `compute_initial()`
#'
#' @noRd
check_id <- function(id = NULL, lambda) {
  if(is.null(id)) {
    type <- gsub("^([a-z]+).*", "\\1", attr(lambda, "type"))
    if(grepl("custom", type)) {return("custom")}
    if(grepl("sigma", type)) {return("sigma")}
    position <- attr(lambda, "position")
    return(paste0(type, "_", position))
  }
  if(!is.character(id)) {stop("Please provide a character scalar.")}
  return(id)
}


#' Obtain exact values of lambda for `compute_initial()`
#'
#' @noRd
get_exact <- function(x, id) {
  if(grepl("tstat", id)) {
    x$model[[paste0("beta_", gsub(".*_([0-9]+)", "\\1", id))]] /
      x$model[[paste0("se_", gsub(".*_([0-9]+)", "\\1", id))]]
  } else {
    x$model[[id]] # Potentially NULL
  }
}


#' Recalculate model quantities quickly for `compute_goal()`
#'
#' @noRd
re_infl <- function(x, rm) {
  if(x$meta$class == "lm") {
    re <- influence_lm(x$meta$model, rm = rm,
      options = list("just_model" = TRUE), cluster = x$meta$cluster)
  } else {
    re <- influence_iv(x$meta$model, rm = rm,
      options = list("just_model" = TRUE), cluster = x$meta$cluster)
  }
  return(re)
}
