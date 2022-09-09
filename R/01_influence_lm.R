
influence_lm <- function(x,
  rm = NULL, options, cluster = NULL, XX_inv = NULL) {

  # Inputs ---

  meta <- list("model" = x, "cluster" = cluster, "class" = "lm")

  data <- get_data(x)
  if(is.null(rm)) { # Reuse existent quantities
    y <- data$y
    X <- data$X

    if(is.null(qr_x <- x$qr)) {qr_x <- qr(X)}
    R <- qr.R(qr_x)
    r_cond <- rcond(R, norm = "I")
  } else { # Compute without rows
    y <- data$y[-rm, drop = FALSE]
    X <- data$X[-rm, , drop = FALSE]
    if(!is.null(cluster)) {cluster <- cluster[-rm, , drop = FALSE]}
    if(is.null(XX_inv)) { # Avoid recomputation if XX_inv is provided
      qr_x <- qr(X)
      if(qr_x$rank != NCOL(X)) {stop("Removal resulted in rank deficiency.")}
      R <- qr.R(qr_x)
      r_cond <- rcond(R, norm = "I")
    } else {
      qr_x <- R <- r_cond <- NULL
    }
  }


  # Regression quantities ---

  N <- NROW(y)
  K <- NCOL(X)
  idx <- seq.int(N)

  # Coefficients
  if(!is.null(qr_x)) {
    beta <- as.numeric(solve_cholesky(R, crossprod(X, y)))
    XX_inv <- chol2inv(R)
  } else {
    beta <- as.numeric(XX_inv %*% crossprod(X, y))
  }
  res <- as.numeric(y - X %*% beta)
  rss <- sum(res^2)
  sigma <- sqrt(rss / (N - K))

  # Standard errors
  if(is.null(cluster)) {
    vcov <- XX_inv * sigma^2 # Plain
  } else {
    veggies <- veggiesCL(res, X, # Clustered
      cluster = cluster, type = "HC1")
    bread <- XX_inv * N
    vcov <- 1 / N * (bread %*% veggies %*% bread)
  }
  se <- sqrt(diag(vcov))
  tstat <- beta / se

  # Others -- all assuming there's an intercept
  r2 <- 1 - rss / sum((y - mean(y))^2)
  fstat <- r2 / (1 - r2) * (N - K) / min(K - 1, 1)
  ll <- 0.5 * (-N * (log(2 * pi) + 1 - log(N) + log(rss)))
  # aic <- -2 * ll + 2 * (K + 1)
  # bic <- -2 * ll + (K + 1) * log(N)

  model <- list("beta" = beta, "sigma" = sigma, "se" = se, "tstat" = tstat,
    "r2" = r2, "fstat" = fstat, "ll" = ll,
    "qr" = qr_x, "XX_inv" = XX_inv, "r_cond" = r_cond)

  if(isTRUE(options$just_model)) {
    return(structure(list("model" = model, "meta" = meta), class = "influence"))
  }


  # Influence quantities ---

  # Diagonal of the hat matrix
  hat <- if(isTRUE(options$hat)) {
    if(!is.null(qr_x)) {
      matrix(rowSums(qr.Q(qr_x)^2)) # Make use of QQ'
    } else {
      matrix(vapply(idx, function(i) { # Diagonal of X(X'X)⁻¹X'
        X[i, ] %*% XX_inv %*% X[i, ]}, numeric(1L)))
    }
  } else {
    matrix(0, N) # Skip calculation
  }

  # DFBETA
  beta_i <- if(isTRUE(options$beta)) {
    if(!is.null(R)) {
      t(t(t(-solve_cholesky(R, t(X))) *
        ifelse(hat == 1, 0, res / (1 - hat))[, 1L]) + beta)
    } else {
      t(t(t(-tcrossprod(XX_inv, X)) *
        ifelse(hat == 1, 0, res / (1 - hat))[, 1L]) + beta)
    }
  } else { # BGM's derivative
    if(!is.null(R)) {
      t(t(t(-solve_cholesky(R, t(X))) * res) + beta)
    } else {
      t(t(t(-tcrossprod(XX_inv, X)) * res) + beta)
    }
  }

  # Sigma when dropping observation i
  sigma_i <- if(isTRUE(options$sigma)) {
    matrix(sqrt((rss - res^2 / ifelse(hat == 1, 1, (1 - hat))) / (N - K - 1)))
  } else { # BGM's derivative
    matrix(res^2 / (N - K) + sigma)
  }

  # Standard errors
  if(isTRUE(options$se) && (is.null(cluster) || isFALSE(options$cluster))) {
    se_i <- vapply(idx, function(i) {
      inv_i <- tryCatch(update_inv(XX_inv, X[i, , drop = FALSE]),
        error = function(e) {matrix(NaN, nrow(X), ncol(X))})
      sqrt(diag(inv_i) * (sigma_i[i])^2)
    }, numeric(K))
    if(K != 1) {se_i <- t(se_i)} else {se_i <- matrix(se_i)}
  } else if(isTRUE(options$se) && !is.null(cluster)) {
    se_i <- vapply(idx, function(i) {
      res_i <- res[-i] + as.numeric(X[-i, ] %*% (beta - beta_i[i, ]))
      veggies_i <- veggiesCL(res_i, X[-i, , drop = FALSE],
        cluster = cluster[-i, , drop = FALSE], type = "HC1")
      inv_i <- tryCatch(update_inv(XX_inv, X[i, , drop = FALSE]),
        error = function(e) {matrix(NaN, nrow(X), ncol(X))})
      bread_i <- inv_i * (N - 1)
      sqrt(diag(1 / (N - 1) * (bread_i %*% veggies_i %*% bread_i)))
    }, numeric(K))
    if(K != 1) {se_i <- t(se_i)} else {se_i <- matrix(se_i)}
  } else { # BGM's derivative
    sigma_dbeta <- -2 * colSums(res * X) / (N - K)
    R_x <- if(!is.null(R)) {solve_cholesky(R, t(X))} else {XX_inv %*% t(X)}
    veggies_i <- (-R_x^2) * sigma^2 +
      outer(diag(XX_inv), (sigma_i[, 1] - sigma))
    sigma_d <- colSums(sigma_dbeta * (t(beta_i) - beta)) +
      (sigma_i[, 1] - sigma)
    bread_i <- veggies_i + outer(diag(XX_inv),
      colSums(sigma_dbeta * -(t(beta_i) - beta)))
    se_i <- t(-0.5 * bread_i / sqrt(diag(vcov)) + se)
  }

  # t value
  tstat_i <- if(isTRUE(options$tstat)) {
    matrix(beta_i / se_i, N, K)
  } else { # Doesn't enforce exact beta and se
    matrix(beta_i / se_i, N, K)
  }

  # DFFITS
  dffits <- if(isTRUE(options$dffits)) {
    matrix(res * sqrt(hat) / (sigma_i * (1 - pmin(1, hat))))
  } else {
    NULL
  }

  # Cook's distance
  cooksd <- if(isTRUE(options$cooksd)) {
    matrix(((res / ((1 - pmin(1, hat)) * sigma))^2 * hat) / K)
  } else {
    NULL
  }

  # Studentised residual
  rstudent <- if(isTRUE(options$rstudent)) {
    matrix(res / (sigma_i * sqrt(1 - pmin(1, hat))))
  } else {
    NULL
  }

  # Covratio
  covratio <- if(isTRUE(options$covratio)) {
    matrix(1 / ((1 - pmin(1, hat)) *
      ((N - K - 1 + (res / (sigma_i * sqrt(1 - pmin(1, hat))))^2) / (N - K))^K))
  } else {
    NULL
  }


  # Return ---

  structure(list(
    "model" = model,
    "hat" = hat,
    "beta_i" = beta_i, "sigma_i" = sigma_i, "se_i" = se_i, "tstat_i" = tstat_i,
    "cooksd" = cooksd, "dffits" = dffits,
    "rstudent" = rstudent, "covratio" = covratio,
    "meta" = meta
  ), class = "influence")
}
