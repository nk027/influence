
#' @noRd
sensitivity_iv <- function(x,
  lambda = set_lambda(), options = set_options(),
  cluster = NULL, verbose = TRUE) {

  # Inputs ---

  verbose <- isTRUE(verbose)
  meta <- list("lambda" = lambda, "options" = options, "cluster" = cluster,
    "model" = x, "class" = "ivreg")

  # Dimensions
  N <- x$n
  K <- if(is.null(x$p)) {x$rank} else {x$p}

  # Iterations
  n_max <- check_iterations(N, options$n_max, options$p_max)
  # Cluster for clustered standard errors
  cluster <- check_cluster(cluster, N)
  if(!any(options$fwl == 0) && options$sm_re != 1L) {
    warning("Frisch-Waugh-Lovell and Sherman-Morrison not implemented for IV.")
  }


  # Start ---

  # First calculation
  step <- influence_iv(x, options = options, cluster = cluster)
  rank <- rank_influence(step, lambda = lambda)

  # Prepare outputs
  idx <- seq.int(N)
  rm <- vector("numeric", n_max)
  obs <- rank[seq.int(0L, n_max + 1L), "order"]
  # List with values wrt lambda, the model, and the initial approximation
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
    "meta" = meta
  ), class = "sensitivity")

  # Fill for step one
  rm[1L] <- rank[1L, "order"]
  out$model[1, ] <- c(N,
    step$model$sigma, step$model$beta, step$model$se,
    step$model$r2, step$model$fstat,
    step$model$r2_first, step$model$fstat_first)


  # Iterate ---

  start <- Sys.time()
  if(verbose) {pb <- txtProgressBar(min = 2L, max = n_max, style = 3L)}

  # Loop start >
  for(i in seq.int(2L, n_max + 1L)) {

    # No FWL or SM for IV models
    step <- tryCatch(influence_iv(x,
      rm = rm, options = options, cluster = cluster),
      error = function(e) {
        message("\nComputation failed at step ", i, " with:\n", e); e
      })
    if(inherits(step, "error")) {break} # Exit loop
    rank <- rank_influence(step, lambda)
    # Find observation to remove next
    rm[i] <- idx[-rm][rank[1L, "order"]] # Index kept constant
    rm_val <- rank[rank[1L, "order"], "value"]

    # Store results
    out$influence$id[i] <- rm[i]
    out$influence$lambda[i] <- rm_val
    out$model[i, ] <- c(N - i + 1,
      step$model$sigma, step$model$beta, step$model$se,
      step$model$r2, step$model$fstat,
      step$model$r2_first, step$model$fstat_first)

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
