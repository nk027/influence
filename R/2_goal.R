
goal <- function(x, ...) {{UseMethod("goal", x)}}

goal.lm <- function(x,
  lambda = set_lambda(), delta = set_delta(),
  n_upper = NULL, n_lower = 0L, target = NULL,
  options = set_compute(), cluster = NULL) {

  x <- infl.lm(x, options = options, cluster = cluster)
  goal(x, lambda = lambda, delta = delta,
    n_upper = n_upper, n_lower = n_lower, target = target)
}

goal.ivreg <- function(x,
  lambda = set_lambda(), delta = set_delta(),
  n_upper = NULL, n_lower = 0L, target = NULL,
  options = set_compute(), cluster = NULL) {

  x <- infl.ivreg(x, options = options, cluster = cluster)
  goal(x, lambda = lambda, delta = delta,
    n_upper = n_upper, n_lower = n_lower, target = target)
}

goal.influence <- function(x,
  lambda = set_lambda(), delta = set_delta(),
  n_upper = NULL, n_lower = 0L, target = NULL) {

  target(x, lambda = lambda, delta = delta,
    n_upper = n_upper, n_lower = n_lower, target = target)
}


target <- function(x,
  lambda = set_lambda(), delta = set_delta(),
  n_upper = NULL, n_lower = 0L, target = NULL) {

  if(is.null(target)) {
    target <- attr(lambda, "target")
  } else {
    target <- num_check(target, -Inf, Inf, msg = "Please check the target.")
  }

  if(is.null(n_upper)) {
    initial <- init(x, lambda = lambda)
    n_upper <- which(delta(initial$initial, target))[1L] - 1L
    if(is.na(n_upper)) {
      stop("Target not within the initial approximation's reach. ",
        "Consider setting 'n_upper' manually.")
    }
  } else {
    n_upper <- num_check(n_upper, 1L, NROW(x$hat),
      msg = "Choose a valid upper bound for observations to remove.")
  }

  rank <- rank_influence(x, lambda)
  rm <- rank[, "order"]

  # First step -- check the target is achievable
  re <- re_infl(x, rm[seq.int(n_upper)])
  if(!delta(init(re, lambda)$exact[1L], target)) {
    stop("Target change not achieved at the upper bound.")
  }

  while(n_lower <= n_upper) {
    n_consider <- floor((n_lower + n_upper) / 2)
    re <- re_infl(x, rm[seq.int(n_consider)])
    if(!delta(init(re, lambda)$exact[1L], target)) {
      n_lower <- n_consider + 1L
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
    re <- influence_lm(x$meta$model, rm = rm,
      options = list("just_model" = TRUE), cluster = x$meta$cluster)
  } else {
    re <- influence_iv(x$meta$model, rm = rm,
      options = list("just_model" = TRUE), cluster = x$meta$cluster)
  }
  return(re)
}
