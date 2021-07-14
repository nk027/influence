
sens <- function(x, ...) {{UseMethod("sens", x)}}


sens.lm <- function(x,
  lambda = set_lambda(),
  options = set_options(),
  cluster = NULL,
  verbose = TRUE) {

  # Inputs ---

  verbose <- isTRUE(verbose)

  data <- mdl_to_mat(x)
  y <- data$y
  X <- data$X

  N <- NROW(y)
  K <- NCOL(X)

  # Iterations
  n_max <- check_iterations(N, options$n_max, options$p_max)

  # Cluster for clustered standard errors
  cluster <- check_cluster(cluster, N)

  # Reduce covariates using the Frisch-Waugh-Lovell theorem
  if(!any(options$fwl == 0)) {
    if(any(options$fwl > K)) {
      warning("No variables to marginalise using FWL found.")
    } else {
      fwl <- frisch_waugh_lovell(data$X, data$y, variables = options$fwl)
      y <- fwl$y
      X <- fwl$X
      K <- NCOL(X)
    }
  } # Reapplication later is determined by options$fwl_re


  # Start ---

  step <- influence_lm(X, y, options = options, cluster = cluster)
  rank <- rank_influence(step, lambda)

  # Prepare outputs
  idx <- seq.int(N)
  rm <- vector("numeric", n_max)

  out <- structure(list(
    "influence" = data.frame(
      "N" = seq.int(N, N - n_max),
      "id" = c(rank[1L, "order"], rep(NA_integer_, n_max)),
      "lambda" = c(rank[rank[1L, "order"], "value"], rep(NA_real_, n_max))
    ),
    "model" = as.data.frame(matrix(
      NA_real_, n_max + 1L, 2 + K * 2 + 3,
      dimnames = list(NULL, c("N", "sigma",
        paste0("beta_", seq.int(K)), paste0("se_", seq.int(K)),
        "R2", "F", "LL")
      )
    )),
    "initial" = data.frame(
      "id" = rank[, "order"], "lambda" = rank[rank[, "order"], "value"]
    ),
    "meta" = list("lambda" = lambda, "options" = options, "cluster" = cluster)
  ), class = "sensitivity")

  rm[1L] <- rank[1L, "order"]
  out$model[1, ] <- c(N,
    step$lm$sigma, step$lm$beta, step$lm$se,
    step$lm$r2, step$lm$fstat, step$lm$ll)

  if(n_max == 1L) {return(out)}


  # Iterate ---

  start <- Sys.time()
  if(verbose) {pb <- txtProgressBar(min = 2L, max = n_max, style = 3L)}

  # Loop
  for(i in seq.int(2L, n_max + 1L)) {

    # Reorthogonalise FWL
    if(!any(options$fwl == 0) && (i - 1L) %% options$fwl_re == 0) {
      fwl <- frisch_waugh_lovell(data$X[-rm, , drop = FALSE],
        data$y[-rm, drop = FALSE], variables = options$fwl)
      y[-rm] <- fwl$y
      X[-rm, ] <- fwl$X
    }

    # Update using QR or Sherman-Morrison
    if((i - 1L) %% options$sm_re == 0) {
      step <- tryCatch(influence_lm(
        X[-rm, , drop = FALSE], y[-rm, drop = FALSE],
        options = options, cluster = cluster[-rm, , drop = FALSE]),
        error = function(e) {
          warning("Computation failed at step ", i, " with: ", e); e
      })
    } else {
      step <- tryCatch(influence_lm(
        X[-rm, , drop = FALSE], y[-rm, drop = FALSE],
        options = options,
        XX_inv = update_inv(step$lm$XX_inv, X[rm[i - 1L], , drop = FALSE]),
        cluster = cluster[-rm, , drop = FALSE]),
        error = function(e) {
          warning("Computation failed at step ", i, " with: ", e); e
      })
    }
    if(inherits(step, "error")) {break}
    rank <- rank_influence(step, lambda)

    # Store results
    if(options$adaptive) {
      rm[i] <- idx[-rm][rank[1L, "order"]] # Index kept constant
      rm_val <- rank[rank[1L, "order"], "value"]
    } else {
      rm[i] <- out$initial$id[i]
      rm_val <- out$initial$lambda[i]
    }
    out$influence$id[i] <- rm[i]
    out$influence$lambda[i] <- rm_val
    out$model[i, ] <- c(N - i + 1,
      step$lm$sigma, step$lm$beta, step$lm$se,
      step$lm$r2, step$lm$fstat, step$lm$ll)

    if(verbose) {setTxtProgressBar(pb, i)}
  }

  timer <- format(Sys.time() - start)
  if(verbose) {close(pb); cat("Calculations took ", timer, ".\n", sep = "")}


  # Wrap up ---

  out$model <- out$model[!is.na(out$model$N), ]
  out$influence <- out$influence[!is.na(out$influence$id), ]

  return(out)
}


sens.ivreg <- function(x,
  lambda = set_lambda(),
  options = set_options(),
  cluster = NULL,
  verbose = TRUE) {

  # Inputs ---

  verbose <- isTRUE(verbose)

  data <- mdl_to_mat(x)
  if(is.null(data$Z)) {return(sens.lm(x, lambda, options, cluster))}
  y <- data$y
  X <- data$X
  Z <- data$Z

  N <- NROW(y)
  K <- NCOL(X)
  M <- NCOL(Z)

  # Iterations
  n_max <- check_iterations(N, options$n_max, options$p_max)

  # Cluster for clustered standard errors
  cluster <- check_cluster(cluster, N)

  if(!any(options$fwl == 0)) {
    warning("Frisch-Waugh-Lovell not implemented for IV.")
  }


  # Start ---

  step <- influence_iv(X, Z, y, options = options, cluster = cluster)
  rank <- rank_influence(step, lambda)

  # Prepare outputs
  idx <- seq.int(N)
  rm <- vector("numeric", n_max)

  obs <- rank[seq.int(0L, n_max + 1L), "order"]
  out <- structure(list(
    "influence" = data.frame(
      "N" = seq.int(N, N - n_max),
      "id" = c(rank[1L, "order"], rep(NA_integer_, n_max)),
      "lambda" = c(rank[rank[1L, "order"], "value"], rep(NA_real_, n_max))
    ),
    "model" = as.data.frame(matrix(
      NA_real_, n_max + 1L, 2 + K * 2 + 4,
      dimnames = list(NULL, c("N", "sigma",
        paste0("beta_", seq.int(K)), paste0("se_", seq.int(K)),
        "R2", "F", "R2_1st", "F_1st")
      )
    )),
    "initial" = data.frame(
      "id" = rank[, "order"], "lambda" = rank[rank[, "order"], "value"]
    ),
    "meta" = list("lambda" = lambda, "options" = options, "cluster" = cluster)
  ), class = "sensitivity")

  rm[1L] <- rank[1L, "order"]
  out$model[1, ] <- c(N,
    step$lm$sigma, step$lm$beta, step$lm$se,
    step$lm$r2, step$lm$fstat, step$lm$r2_first, step$lm$fstat_first)


  # Iterate ---

  start <- Sys.time()
  if(verbose) {pb <- txtProgressBar(min = 2L, max = n_max, style = 3L)}

  # Loop
  for(i in seq.int(2L, n_max + 1L)) {

    # No FWL or SM for IV models
    step <- tryCatch(influence_iv(X[-rm, , drop = FALSE],
      Z[-rm, , drop = FALSE], y[-rm, drop = FALSE],
      options = options, cluster = cluster[-rm, , drop = FALSE]),
      error = function(e) {
        warning("Computation failed at step ", i, " with: ", e); e
    })
    if(inherits(step, "error")) {break}
    rank <- rank_influence(step, lambda)

    # Store results
    if(options$adaptive) {
      rm[i] <- idx[-rm][rank[1L, "order"]] # Index kept constant
      rm_val <- rank[rank[1L, "order"], "value"]
    } else {
      rm[i] <- out$initial$id[i]
      rm_val <- out$initial$lambda[i]
    }
    out$influence$id[i] <- rm[i]
    out$influence$lambda[i] <- rm_val
    out$model[i, ] <- c(N - i + 1,
      step$lm$sigma, step$lm$beta, step$lm$se,
      step$lm$r2, step$lm$fstat, step$lm$r2_first, step$lm$fstat_first)

    if(verbose) {setTxtProgressBar(pb, i)}
  }

  timer <- format(Sys.time() - start)
  if(verbose) {close(pb); cat("Calculations took ", timer, ".\n", sep = "")}


  # Wrap up ---

  out$model <- out$model[!is.na(out$model$N), ]
  out$influence <- out$influence[!is.na(out$influence$id), ]

  return(out)
}
