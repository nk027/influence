
#' Methods to compute an initial approximation of influence
#'
#' @param ... Options dispatched internally.
#'
#' @return Returns a \code{initial} object.
#'
#' @export
init <- function(x, ...) {{UseMethod("init", x)}}

#' @noRd
init.default <- function(x,
  lambda = set_lambda(), start = NULL,
  options = set_compute(), cluster = NULL) {

  # Cluster for clustered standard errors
  n <- if(!is.null(x$n)) x$n else x$rank + x$df.resid
  cluster <- check_cluster(cluster, n)

  x <- infl(x, options = options, cluster = cluster)
  init(x, lambda = lambda, start = start)
}

# init.ivreg <- function(x,
#   lambda = set_lambda(), start = NULL,
#   options = set_compute(), cluster = NULL) {

#   # Cluster for clustered standard errors
#   cluster <- check_cluster(cluster, x$n)

#   x <- infl.ivreg(x, options = options, cluster = cluster)
#   init.influence(x, lambda = lambda, start = start)
# }

#' @noRd
init.influence <- function(x, lambda = set_lambda(), start = NULL) {

  rank <- rank_influence(x, lambda = lambda)
  out <- create_object(x, rank = rank, lambda = lambda)

  compute_initial(out, start = start)
}

#' @noRd
init.sensitivity <- function(x, start = NULL) {

  compute_initial(x, start = start)
}


#' @noRd
compute_initial <- function(x, start = NULL) {

  lambda_id <- check_id(NULL, lambda = x$meta$lambda)
  exact <- get_exact(x, lambda_id)

  if(is.null(start)) {
    if(all(!grepl(lambda_id, names(x$model))) &&
      !grepl("tstat_[0-9]+", lambda_id) && !grepl("sigma", lambda_id)) {
      warning("Cannot determine starting value for the requested 'lambda'",
        "automatically. Set to zero, consider providing a value via 'start'.")
      start <- 0
    } else {
      start <- exact[1L]
    }
  }

  initial <- if(attr(x$meta$lambda, "sign") == -1L) { # Initial approximation
    cumsum(c(start, -(start + x$initial$lambda)))
  } else {
    cumsum(c(start, -(start - x$initial$lambda)))
  }

  structure(list(
    "initial" = initial, "id" = x$initial$id,
    "exact" = exact, "lambda_id" = lambda_id
  ), class = "initial")
}
