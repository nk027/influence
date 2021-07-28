
init <- function(x, ...) {{UseMethod("init", x)}}


init.influence <- function(x, lambda = set_lambda(), start = NULL) {

  rank <- rank_influence(x, lambda = lambda)
  out <- create_object(x, rank = rank, lambda = lambda)

  initial(out, start = start)
}

init.sensitivity <- function(x, start = NULL) {

  initial(x, start = start)
}


initial <- function(x, start = NULL) {

  id <- check_id(NULL, lambda = x$meta$lambda)
  exact <- get_exact(x, id)

  if(is.null(start)) {
    if(all(!grepl(id, names(x$model))) &&
      !grepl("tstat_[0-9]+", id) && !grepl("sigma", id)) {
      warning("Cannot determine starting value for the requested 'id' ",
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

  return(list("initial" = initial, "exact" = exact, "id" = id))
}


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


get_exact <- function(x, id) {
  if(grepl("tstat", id)) {
    x$model[[paste0("beta_", gsub(".*_([0-9]+)", "\\1", id))]] /
      x$model[[paste0("se_", gsub(".*_([0-9]+)", "\\1", id))]]
  } else {
    x$model[[id]] # Potentially NULL
  }
}
