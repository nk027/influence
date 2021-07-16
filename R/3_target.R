
init <- function(x, ...) {{UseMethod("init", x)}}

init.influence <- function(x, lambda = set_lambda(), id = NULL, ...) {

  rank <- rank_influence(x, lambda)

  # Rebuild a plain version of the influence object
  out <- sens_dummy(x, rank, lambda)
  out$model[1, ] <- c(NROW(x$hat), x$lm$sigma, x$lm$beta, x$lm$se)

  init.sensitivity(out, id = id, ...)
}

init.sensitivity <- function(x, id = NULL, ...) {

  if(is.null(id)) { # Try to guess
    type <- gsub("^([a-z]+).*", "\\1", attr(x$meta$lambda, "type"))
    position <- attr(x$meta$lambda, "position")
    id <- paste0(type, "_", position)
  }
  if(!any(grepl(type, names(x$model)))) {
    stop("Type not supported for automatic calculation.")
  }

  exact <- if(grepl("tstat", id)) { # Exact
    x$model[[paste0("beta_", gsub(".*_([0-9]+)", "\\1", id))]] /
      x$model[[paste0("se_", gsub(".*_([0-9]+)", "\\1", id))]]
  } else {
    x$model[[id]]
  }
  initial <- if(attr(x$meta$lambda, "sign") == -1L) { # Initial
    cumsum(c(exact[1L], -(exact[1L] + x$initial$lambda)))
  } else {
    cumsum(c(exact[1L], -(exact[1L] - x$initial$lambda)))
  }

  list("initial" = initial, "exact" = exact, "id" = id)
}


target <- function(x, ...) {{UseMethod("target", x)}}

target.influence <- function(x,
  lambda = set_lambda(), id = NULL,
  target = NULL, n_upper = NULL, ...) {

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

sens_dummy <- function(x, rank, lambda) {
  list(
    "initial" = data.frame(
      "id" = rank[, "order"], "lambda" = rank[rank[, "order"], "value"]
    ),
    "model" = as.data.frame(matrix(
      NA_real_, 1L, 2 + length(x$lm$beta) * 2,
      dimnames = list(NULL, c("N", "sigma",
        paste0("beta_", seq.int(length(x$lm$beta))),
        paste0("se_", seq.int(length(x$lm$beta))))
      )
    )),
    "meta" = list("lambda" = lambda)
  )
}
