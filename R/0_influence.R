
infl <- function(x, ...) {{UseMethod("infl", x)}}


infl.matrix <- function(x, y, z,
  options = set_compute(), cluster = NULL) {
  if(missing(z)) {
    return(influence_lm(x, y, options = options, cluster = cluster))
  } else {
    return(influence_iv(x, y, z, options = options, cluster = cluster))
  }
}

infl.lm <- function(x, options = set_compute(), cluster = NULL) {
  data <- mdl_to_mat(x)
  y <- data$y
  X <- data$X
  influence_lm(X, y, options = options, cluster = cluster)
}

infl.ivreg <- function(x, options = set_compute(), cluster = NULL) {
  data <- mdl_to_mat(x)
  if(is.null(data$Z)) {infl.lm(x, options = options, cluster = cluster)}
  y <- data$y
  X <- data$X
  Z <- data$Z
  influence_iv(X, Z, y, options = options, cluster = cluster)
}


influence_lm <- function(X, y,
  options, XX_inv = NULL, cluster = NULL) {

  # Regression quantities ---

  N <- NROW(y)
  K <- NCOL(X)
  idx <- seq.int(N)

  # Coefficients
  if(is.null(XX_inv)) {
    qr_x <- qr(X)
    R <- qr.R(qr_x)
    beta <- as.numeric(solve_cholesky(R, crossprod(X, y)))
    XX_inv <- chol2inv(R)
  } else {
    qr_x <- R <- NULL
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
    matrix(1, N) # Skip calculation
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
    # sweep(-t(solve_cholesky(R, t(X))) * # Cleaner but slower
    #   ifelse(hat == 1, 0, res / (1 - hat)), 2, beta, "+")
  } else {
    NULL
  }

  # Sigma when dropping observation i
  sigma_i <- if(isTRUE(options$sigma)) {
    matrix(sqrt((rss - res^2 / ifelse(hat == 1, 1, (1 - hat))) / (N - K - 1)))
  } else {
    matrix(rep(sigma, N))
  }

  # Standard errors
  if(isTRUE(options$se) && is.null(cluster)) {
    se_i <- vapply(idx, function(i) {
      sqrt(diag(update_inv(XX_inv, X[i, , drop = FALSE])) * (sigma_i[i])^2)
    }, numeric(K))
    if(K != 1) {se_i <- t(se_i)} else {se_i <- matrix(se_i)}
  } else if(isTRUE(options$se) && !is.null(cluster)) {
    se_i <- vapply(idx, function(i) {
      res_i <- res[-i] + as.numeric(X[-i, ] %*% (beta - beta_i[i, ]))
      veggies_i <- veggiesCL(res_i, X[-i, , drop = FALSE],
        cluster = cluster[-i, , drop = FALSE], type = "HC1")
      bread_i <- update_inv(XX_inv, X[i, , drop = FALSE]) * (N - 1)
      sqrt(diag(1 / (N - 1) * (bread_i %*% veggies_i %*% bread_i)))
    }, numeric(K))
    if(K != 1) {se_i <- t(se_i)} else {se_i <- matrix(se_i)}
  } else {
    NULL
  }

  # t value
  tstat_i <- if(isTRUE(options$tstat)) {
    matrix(beta_i / se_i, N, K)
  } else {
    NULL
  }

  # DFFITS
  dffits <- if(isTRUE(options$dffits)) {
    matrix(res * sqrt(hat) / (sigma_i * (1 - hat)))
  } else {
    NULL
  }

  # Cook's distance
  cooksd <- if(isTRUE(options$cooksd)) {
    matrix(((res / ((1 - hat) * sigma))^2 * hat) / K)
  } else {
    NULL
  }

  # Studentised residual
  rstudent <- if(isTRUE(options$rstudent)) {
    matrix(res / (sigma_i * sqrt(1 - hat)))
  } else {
    NULL
  }

  # Covratio
  covratio <- if(isTRUE(options$covratio)) {
    matrix(1 / ((1 - hat) *
      ((N - K - 1 + (res / (sigma_i * sqrt(1 - hat)))^2) / (N - K))^K))
  } else {
    NULL
  }


  # Return ---

  out <- structure(list(
    "lm" = list("beta" = beta, "sigma" = sigma,
      "se" = se, "tstat" = tstat,
      "r2" = r2, "fstat" = fstat, "ll" = ll,
      "qr" = qr_x, "XX_inv" = XX_inv),
    "hat" = hat,
    "beta_i" = beta_i, "sigma_i" = sigma_i,
    "se_i" = se_i, "tstat_i" = tstat_i,
    "cooksd" = cooksd, "dffits" = dffits,
    "rstudent" = rstudent, "covratio" = covratio
  ), class = "influence")
}


influence_iv <- function(X, Z, y,
  options, cluster = NULL) {

  # Regression quantities ---

  N <- NROW(y)
  K <- NCOL(X)
  M <- NCOL(Z)
  idx <- seq.int(N)

  # Coefficients
  qr_z <- qr(Z) # Needed for projection, residuals, and stage 1 hat
  X_proj <- qr.fitted(qr_z, X)
  X_resid <- X - X_proj
  qr_x <- qr(X_proj) # Needed for beta, vcov, and stage 2 hat
  beta <- qr.coef(qr_x, y)
  qr_a <- qr(crossprod(X, X_proj)) # Needed for influence
  XX_inv <- chol2inv(qr.R(qr_x))

  res <- as.numeric(y - X %*% beta)
  res_z <- qr.resid(qr_z, y)
  res_p <- as.numeric(y - X_proj %*% beta)
  rss <- sum(res^2)
  sigma <- sqrt(rss / (N - K))

  # Standard errors
  if(is.null(cluster)) {
    vcov <- XX_inv * sigma^2 # Plain
  } else {
    veggies <- veggiesCL(res, X_proj, # Clustered
      cluster = cluster, type = "HC0")
    bread <- XX_inv * N
    vcov <- 1 / N * (bread %*% veggies %*% bread)
  }
  se <- sqrt(diag(vcov))
  tstat <- beta / se

  # Others -- all assuming there's an intercept
  r2 <- 1 - rss / sum((y - mean(y))^2)
  fstat <- r2 / (1 - r2) * (N - K) / (K - 1)
  pos_endo <- apply(abs(X_resid), 2, max) > 1e-8
  r2_first <- 1 - sum(X_resid[, pos_endo]^2) /
    sum((X[, pos_endo] - colMeans(X[, pos_endo, drop = FALSE]))^2)
  fstat_first <- r2_first / (1 - r2_first) * (N - M) / (M - 1)


  # Influence quantities ---

  # Diagonal of hat matrices ( X(X'PX)⁻¹X', R(X'PX)⁻¹X', (X'PX)⁻¹R' )
  Ai_X <- qr.solve(qr_a, t(X))
  Ai_Xr <- qr.solve(qr_a, t(X_resid))
  h_X <- vapply(idx, function(i) {X[i, ] %*% Ai_X[, i]}, numeric(1L))
  h_XrX <- vapply(idx, function(i) {X_resid[i, ] %*% Ai_X[, i]}, numeric(1L))
  h_XrXr <- vapply(idx, function(i) {X_resid[i, ] %*% Ai_Xr[, i]}, numeric(1L))
  h_Pz <- rowSums(qr.Q(qr_z)^2) # Diagonal of the projection matrix P

  hat <- if(options$hat) {
    cbind("stage1" = h_Pz, "stage2" = rowSums(qr.Q(qr_x)^2), "projection" = h_X)
  } else {
    matrix(1, N, 3, dimnames = list(NULL, c("stage1", "stage2", "projection")))
  }

  # DFBETA, loosely following Phillips (1977)
  denom <- (1 - h_Pz + h_XrXr)
  delta <- 1 - h_X + (h_XrX^2) / denom
  h <- ((1 - h_X) * X_resid + (h_XrX) * X) / (denom * delta)
  j <- ((h_XrX * X_resid) / denom - X) / delta
  g <- h * as.numeric(res_z - res_p) + (h + j) * res
  dfb <- t(qr.solve(qr_a, t(g)))
  beta_i <- t(t(dfb) + beta)

  # Sigma when dropping observation i
  if(options$sigma) {
    XX <- crossprod(X)
    XR <- crossprod(X, res)
    rss_i <- rss + vapply(idx, function(i) {
        dfb[i, ] %*% update_cp(XX, X[i, , drop = FALSE]) %*% dfb[i, ]
      }, numeric(1L)) - 2 * vapply(idx, function(i) {
        dfb[i, ] %*% update_cp(XR, X[i, , drop = FALSE], res[i])
      }, numeric(1L)) - res^2
    sigma_i <- matrix(sqrt(rss_i / (N - K - 1L)))
  } else {
    sigma_i <- matrix(sigma, N)
  }

  # Standard errors
  if(isTRUE(options$se) && is.null(cluster)) {
    ZX <- crossprod(Z, X)
    ZZ_inv <- chol2inv(qr.R(qr_z))
    se_i <- vapply(idx, function(i) {
      proj_i <- Z[-i, ] %*% update_inv(ZZ_inv, Z[i, , drop = FALSE]) %*%
        update_cp(ZX, Z[i, , drop = FALSE], X[i, , drop = FALSE])
      sqrt(diag(chol2inv(chol(crossprod(proj_i)))) * sigma_i[i]^2)
    }, numeric(K))
    if(K != 1) {se_i <- t(se_i)} else {se_i <- matrix(se_i)}
  } else if(isTRUE(options$se) && !is.null(cluster)) {
    ZX <- crossprod(Z, X)
    ZZ_inv <- chol2inv(qr.R(qr_z))
    se_i <- vapply(idx, function(i) {
      res_i <- res[-i] + as.numeric(X[-i, ] %*% t(beta_i[i, , drop = FALSE]))
      proj_i <- Z[-i, ] %*% update_inv(ZZ_inv, Z[i, , drop = FALSE]) %*%
        update_cp(ZX, Z[i, , drop = FALSE], X[i, , drop = FALSE])
      veggies_i <- veggiesCL(res_i, proj_i,
        cluster = cluster[-i, , drop = FALSE], type = "HC0")
      bread_i <- chol2inv(chol(crossprod(proj_i))) * (N - 1)
      sqrt(diag(1 / (N - 1) * (bread_i %*% veggies_i %*% bread_i)))
    }, numeric(K))
    if(K != 1) {se_i <- t(se_i)} else {se_i <- matrix(se_i)}
  } else {
    se_i <- NULL
  }

  # t value
  tstat_i <- if(isTRUE(options$tstat)) {
    matrix(beta_i / se_i)
  } else {
    NULL
  }

  # DFFITS
  dffits <- if(isTRUE(options$dffits) || isTRUE(options$cooksd)) {
    matrix(vapply(idx, function(i) {
      X[i, ] %*% beta_i[i, ] / (sigma_i[i] * sqrt(h_X[i]))
    }, numeric(1L)))
  } else {
    NULL
  }

  # Cook's distance
  cooksd <- if(isTRUE(options$cooksd)) {
    matrix((sigma_i^2 / sigma^2) * dffits^2 / K)
  } else {
    NULL
  }

  # Studentised residual
  rstudent <- if(isTRUE(options$rstudent)) {
    matrix(res / (sigma_i * sqrt(1 - hat[, "stage2"])))
  } else {
    NULL
  }

  # Covratio
  covratio <- if(isTRUE(options$covratio)) {
    matrix(1 / ((1 - hat[, "stage2"]) *
      ((N - K - 1 + (res / (sigma_i * sqrt(1 - hat[, "stage2"])))^2) /
      (N - K))^K))
  } else {
    NULL
  }


  # Return ---

  out <- structure(list(
    "lm" = list("beta" = beta, "sigma" = sigma,
      "se" = se, "tstat" = tstat,
      "r2" = r2, "fstat" = fstat,
      "r2_first" = r2_first, "fstat_first" = fstat_first,
      "qr_z" = qr_z, "qr_x" = qr_x, "qr_a" = qr_a),
    "hat" = hat,
    "beta_i" = beta_i, "sigma_i" = sigma_i,
    "se_i" = se_i, "tstat_i" = tstat_i,
    "cooksd" = cooksd, "dffits" = dffits,
    "rstudent" = rstudent, "covratio" = covratio
  ), class = "influence")
}
