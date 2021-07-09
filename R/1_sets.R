
sensitivity <- function(x, ...) {{UseMethod("sensitivity", x)}}

sensitivity.lm <- function(x,
  lambda = lambda(),
  options = set_options(),
  cluster = NULL) {

  # Inputs ---

  # Get matrices from lm, make sure X is a matrix
  data <- mdl_to_mat(x)
  y <- data$y
  X <- data$X

  N <- NROW(y)
  K <- NCOL(X)

  # Cluster for clustered standard errors
  if(!is.null(cluster)) {
    cluster <- as.data.frame(cluster)
    if(NROW(cluster) != N) {stop("Size of 'cluster' does not match the data.")}
    if(anyNA(cluster)) {stop("No missing 'cluster' values are allowed.")}
  }

  # Reduce covariates using the Frisch-Waugh-Lovell theorem
  if(!any(options$fwl == 0)) {
    if(any(options$fwl > K)) {
      warning("No variables to marginalise found.")
    } else {
      fwl <- frisch_waugh_lovell(data$X, data$y, variables = options$fwl)
      y <- fwl$y
      X <- fwl$X
      K <- NCOL(X)
    }
  } # Reapplication later is determined by options$fwl_re


  # Start ---

  step <- influence_lm(X, y, opt = options$sensitivity, cluster = cluster)
  rank <- rank_influence(step, lambda)

  # Prepare outputs
  idx <- seq.int(N)
  rm <- vector("numeric", n_max)
  obs <- rank[seq.int(0L, n_max + 1L), "order"]
  out <- structure(list(
    "influence" = data.frame(
      "N" = seq.int(N, N - n_max),
      "id" = c(obs[1L], rep(NA_integer_, n_max)),
      "lambda" = c(rank[obs[1L], "value"], rep(NA_real_, n_max)),
      "init_id" = obs, "init_lambda" = rank[obs, "value"]
    ),
    "model" = matrix(NA_real_, n_max + 1L, 2 + K * 2 + 3,
      dimnames = list(NULL, c("N", "sigma",
        paste0("beta-", seq.int(K)), paste0("se-", seq.int(K)),
        "R2", "F", "LL")
      )
    )
  ), class = "influence_sensitivity")
  rm[1L] <- rank[1L, "order"]
  out$model[1, ] <- c(N, step$lm$sigma, step$lm$beta, step$lm$se,
    step$lm$r2, step$lm$fstat, step$lm$ll)


  # Iterate ---

  start <- Sys.time()
  if(verbose) {pb <- txtProgressBar(min = 2L, max = n_max, style = 3L)}

  for(i in seq.int(2L, n_max + 1L)) {

    # Reorthogonalise using FWL
    if((i - 1L) %% options$fwl_re == 0 && !any(options$fwl == 0)) {
      fwl <- frisch_waugh_lovell(data$X[-rm, ], data$y[-rm],
        variables = options$fwl)
      y[-rm] <- fwl$y
      X[-rm, ] <- fwl$X
    }

    # Update using Sherman-Morrison or a full QR recalculation
    if((i - 1L) %% options$sm_re == 0) {
      step <- influence_lm(X[-rm, ], y[-rm],
        opt = options$sensitivity, cluster = cluster[-rm, ])
    } else {
      step <- influence_lm(X[-rm, ], y[-rm], opt = options$sensitivity,
        XX_inv = update_inv(step$lm$XX_inv, X[rm[i - 1L], , drop = FALSE]),
        cluster = cluster[-rm, ])
    }
    rank <- rank_influence(step, lambda)

    rm[i] <- idx[-rm][rank[1L, "order"]] # Index kept constant
    out$influence$id[i] <- rm[i]
    out$influence$lambda[i] <- rank[rm[i], "value"]
    out$model[i, ] <- c(N - i + 1, step$lm$sigma, step$lm$beta, step$lm$se,
      step$lm$r2, step$lm$fstat, step$lm$ll)

    if(verbose) {setTxtProgressBar(pb, i)}
  }

  timer <- format(Sys.time() - start)
  if(verbose) {close(pb); cat("Calculations took ", timer, ".\n", sep = "")}


  # Wrap up ---

  return(out)
}
