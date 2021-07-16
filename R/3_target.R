
target <- function(x, ...) {{UseMethod("target", x)}}

target.influence <- function(x,
  lambda = set_lambda(),
  target = NULL, n_upper = NULL, n_lower = 0L, ...) {

  if(is.null(target)) {
    target <- attr(lambda, "target")
  } else {
    target <- num_check(target, -Inf, Inf, msg = "Please check the target.")
  }

  # How many to try removing?
  initial <- init(x, lambda = lambda)
  N <- length(initial$initial)
  if(is.null(n_upper)) {
    n_upper <- which(initial$initial <= target)[1L] - 1L
    if(is.na(n_upper)) {stop("Target not within the approximation's reach.")}
  } else {
    n_upper <- num_check(n_upper, 1L, N, msg = "Choose a valid maximum.")
  }
  n_lower <- 0L
  # Which ones to remove?
  rank <- rank_influence(x, lambda)
  rm <- rank[, "order"]

  # First step
  re <- re_infl(x, rm[seq.int(n_upper)])
  re_initial <- init(re, lambda)
  if(re_initial$exact[1L] > target) {
    stop("Upper bound not sufficient to achieve the target change.")
  }
  while(n_lower < n_upper) {
    n_consider <- floor((n_lower + n_upper) / 2)
    re <- re_infl(x, rm[seq.int(n_consider)])
    if(init(re, lambda)$exact[1L] > target) {
      n_lower <- n_consider
    } else {
      n_upper <- n_consider - 1L
    }
  }

  return(list(
    "n_removed" = n_consider, "p_removed" = n_consider / N, "target" = target,
    "id" = rm, "lambda" = rank[, "value"], "result" = re
  ))
}

re_infl <- function(x, rm) {
  if(x$meta$class == "lm") {
    re <- influence_lm(x$meta$X[-rm, ], x$meta$y[-rm],
      options = set_options("none"), cluster = x$meta$cluster[-rm, ])
  } else {
    re <- influence_iv(x$meta$X[-rm, ], x$meta$Z[-rm, ], x$meta$y[-rm],
      options = set_options("none"), cluster = x$meta$cluster[-rm, ])
  }
  return(re)
}
