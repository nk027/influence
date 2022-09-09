
#' @noRd
sensitivity_lm <- function(x,
  lambda = set_lambda(), options = set_options(),
  cluster = NULL, verbose = TRUE) {

  # Inputs ---

  verbose <- isTRUE(verbose)
  meta <- list("lambda" = lambda, "options" = options, "cluster" = cluster,
    "model" = x, "class" = "lm")

  # Dimensions
  K <- x$rank
  N <- K + x$df.residual
  # Iterations
  n_max <- check_iterations(N, options$n_max, options$p_max)
  # Cluster for clustered standard errors
  cluster <- check_cluster(cluster, N)
  # Get data
  data <- get_data(x)

  # Reduce covariates using the Frisch-Waugh-Lovell theorem
  if(!any(options$fwl == 0)) {
    if(any(options$fwl > K)) {
      warning("No variables to marginalise using FWL found.")
    } else {
      x <- update_fwl(data$X, data$y, variables = options$fwl)
      K <- NCOL(x$X)
    }
  } # Reapplication later is determined by options$fwl_re


  # Start ---

  # First calculation
  step <- influence_lm(x, options = options, cluster = cluster)
  rank <- rank_influence(step, lambda = lambda)

  # Prepare outputs
  idx <- seq.int(N)
  rm <- vector("numeric", n_max)
  # List with values wrt lambda, the model, and the initial approximation
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
    "meta" = meta
  ), class = "sensitivity")

  # Fill for step one
  rm[1L] <- rank[1L, "order"]
  out$model[1, ] <- c(N,
    step$model$sigma, step$model$beta, step$model$se,
    step$model$r2, step$model$fstat, step$model$ll)

  # We're done if there's only one removal
  if(n_max == 1L) {return(out)}


  # Iterate ---

  start <- Sys.time()
  if(verbose) {pb <- txtProgressBar(min = 2L, max = n_max, style = 3L)}

  # Loop start >
  for(i in seq.int(2L, n_max + 1L)) {

    # Reorthogonalise FWL
    if(!any(options$fwl == 0) && (i - 1L) %% options$fwl_re == 0) {
      x <- update_fwl(data$X, data$y, variables = options$fwl, rm = rm)
    }
    # Update XX_inv if we're using Sherman-Morrison
    XX_inv <- if((i - 1L) %% options$sm_re != 0) {
      tryCatch(update_inv(step$model$XX_inv,
        X_rm = get_data(x)$X[rm[i - 1L], , drop = FALSE]),
        error = function(e) {
          chol2inv(qr.R(qr(X[-rm, , drop = FALSE]))) # More robust QR
        })
    } else {NULL}
    # Calculate new values
    step <- tryCatch(influence_lm(x,
      rm = rm, options = options, cluster = cluster, XX_inv = XX_inv),
      error = function(e) {
        message("\nComputation failed at step ", i, " with:\n", e); e
      })
    if(inherits(step, "error")) {break} # Exit loop
    rank <- rank_influence(step, lambda)
    # Find observation to remove next
    if(isTRUE(options$adaptive)) {
      rm[i] <- idx[-rm][rank[1L, "order"]] # Index is kept constant
      rm_val <- rank[rank[1L, "order"], "value"]
    } else {
      rm_val <- rank[which(idx[-rm] == out$initial$id[i]), "value"]
      rm[i] <- out$initial$id[i]
    }

    # Store results
    out$influence$id[i] <- rm[i]
    out$influence$lambda[i] <- rm_val
    out$model[i, ] <- c(N - i + 1,
      step$model$sigma, step$model$beta, step$model$se,
      step$model$r2, step$model$fstat, step$model$ll)

    if(verbose) {setTxtProgressBar(pb, i)}
  }
  # < Loop end

  timer <- format(Sys.time() - start)
  if(verbose) {close(pb); cat("Calculations took ", timer, ".\n", sep = "")}


  # Wrap up ---

  out$model <- out$model[!is.na(out$model$N), ]
  out$influence <- out$influence[!is.na(out$influence$id), ]

  return(out)
}
