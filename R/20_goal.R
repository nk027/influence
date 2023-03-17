
#' Methods to achieve a target sensitivity
#'
#' @param ... Options dispatched internally.
#'
#' @return Returns a \code{goal} object.
#'
#' @export
goal <- function(x, ...) {{UseMethod("goal", x)}}

#' @noRd
#' @export
goal.default <- function(x,
  lambda = set_lambda(),
  target = set_target(), n_upper = NULL, n_lower = 1L,
  options = set_compute(), cluster = NULL) {

  # Still buggy
  warning("This method is still very buggy.")

  # Cluster for clustered standard errors
  N <- if(!is.null(x$n)) x$n else x$rank + x$df.resid
  cluster <- check_cluster(cluster, N)

  x <- infl(x, options = options, cluster = cluster)
  goal(x, lambda = lambda, target = target,
    n_upper = n_upper, n_lower = n_lower)
}

# goal.ivreg <- function(x,
#   lambda = set_lambda(), target = set_target(),
#   n_upper = NULL, n_lower = 0L, options = set_compute(), cluster = NULL) {

#   x <- infl.ivreg(x, options = options, cluster = cluster)
#   goal.influence(x, lambda = lambda, target = target,
#     n_upper = n_upper, n_lower = n_lower)
# }

#' @noRd
#' @export
goal.influence <- function(x,
  lambda = set_lambda(), target = set_target(),
  n_upper = NULL, n_lower = 1L) {

  compute_goal(x, lambda = lambda, target = target,
    n_upper = n_upper, n_lower = n_lower)
}


#' @noRd
#' @export
compute_goal <- function(x,
  lambda = set_lambda(), target = set_target(),
  n_upper = NULL, n_lower = 1L) {

  value <- attr(target, "target")
  N <- NROW(x$hat)

  if(is.null(n_upper)) {
    initial <- init(x, lambda = lambda)
    n_upper <- which(target(initial$initial, value))[1L] - 1L
    if(is.na(n_upper)) {
      stop("Target not within the initial approximation's reach. ",
        "Consider setting 'n_upper' manually.")
    }
  } else {
    n_upper <- num_check(n_upper, 1L, N,
      msg = "Choose a valid upper bound for observations to remove.")
  }

  rank <- rank_influence(x, lambda)
  rm <- rank[, "order"]

  # First step -- check the target is achievable
  re <- re_infl(x, rm[seq.int(n_upper)])
  achievable <- target(init(re, lambda)$exact[1L], value)
  if(!achievable) {stop("Target change not achieved at the upper bound.")}

  while(n_lower <= n_upper) {
    n_consider <- floor((n_lower + n_upper) / 2)
    re <- re_infl(x, rm[seq.int(n_consider)])
    if(!target(init(re, lambda)$exact[1L], value)) {
      n_lower <- n_consider + 1L
    } else {
      n_upper <- n_consider - 1L
    }
  }

  structure(list(
    "n_removed" = n_consider, "p_removed" = n_consider / N, "target" = target,
    "lambda" = rank[, "value"], "id" = rm,
    "model" = re
  ), class = "goal")
}
