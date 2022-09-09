
#' @noRd
influence_iv <- function(x,
  rm = NULL, options, cluster = NULL) {

  # Inputs ---

  meta <- list("model" = x, "cluster" = cluster, "class" = "ivreg")

  data <- get_data(x)
  if(is.null(rm)) {
    y <- data$y
    X <- data$X
    Z <- data$Z

    if(is.null(qr_z <- x$qr1)) {qr_z <- qr(Z)}
    X_proj <- qr.fitted(qr_z, X)
    X_resid <- X - X_proj
    if(is.null(qr_x <- x$qr)) {qr_x <- qr(X_proj)}
    R <- qr.R(qr_x)
  } else {
    y <- data$y[-rm, drop = FALSE]
    X <- data$X[-rm, , drop = FALSE]
    Z <- data$Z[-rm, , drop = FALSE]

    qr_z <- qr(Z)
    if(qr_z$rank != NCOL(Z)) {stop("Removal resulted in loss of full rank.")}
    X_proj <- qr.fitted(qr_z, X)
    X_resid <- X - X_proj
    qr_x <- qr(X_proj)
    if(qr_x$rank != NCOL(X)) {stop("Removal resulted in loss of full rank.")}
    R <- qr.R(qr_x)
    if(!is.null(cluster)) {cluster <- cluster[-rm, , drop = FALSE]}
  }


  # Regression quantities ---

  N <- NROW(y)
  K <- NCOL(X)
  M <- NCOL(Z)
  idx <- seq.int(N)

  # Coefficients
  beta <- qr.coef(qr_x, y)
  qr_a <- qr(crossprod(X, X_proj), tol = 1e-8)
  XX_inv <- chol2inv(R)

  res <- as.numeric(y - X %*% beta)
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

  model <- list("beta" = beta, "sigma" = sigma, "se" = se, "tstat" = tstat,
    "r2" = r2, "fstat" = fstat,
    "r2_first" = r2_first,  "fstat_first" = fstat_first,
    "qr_z" = qr_z, "qr_x" = qr_x, "qr_a" = qr_a)

  if(isTRUE(options$just_model)) {
    return(structure(list("model" = model, "meta" = meta), class = "influence"))
  }


  # Influence quantities ---

  if(options$beta) {
    res_z <- qr.resid(qr_z, y)
    res_p <- as.numeric(y - X_proj %*% beta)
    # Diagonal of hat matrices ( X(X'PX)⁻¹X', R(X'PX)⁻¹X', (X'PX)⁻¹R' )
    Ai_X <- qr.solve(qr_a, t(X))
    Ai_Xr <- qr.solve(qr_a, t(X_resid))
    h_X <- vapply(idx, function(i) {X[i, ] %*% Ai_X[, i]}, numeric(1L))
    h_XrX <- vapply(idx, function(i) {X_resid[i, ] %*% Ai_X[, i]}, numeric(1L))
    h_XrXr <- vapply(idx, function(i) {X_resid[i, ] %*% Ai_Xr[, i]}, numeric(1L))
    h_Pz <- rowSums(qr.Q(qr_z)^2) # Diagonal of the projection matrix P
    # DFBETA, loosely following Phillips (1977)
    denom <- (1 - h_Pz + h_XrXr)
    delta <- 1 - h_X + (h_XrX^2) / denom
    h <- ((1 - h_X) * X_resid + (h_XrX) * X) / (denom * delta)
    j <- ((h_XrX * X_resid) / denom - X) / delta
    g <- h * as.numeric(res_z - res_p) + (h + j) * res
    dfb <- t(qr.solve(qr_a, t(g)))
    beta_i <- t(t(dfb) + beta)
  } else { # BGM's derivative
    dfb <- t(-solve(qr(crossprod(Z, X)), t(Z * res)))
    beta_i <- t(t(dfb) + beta)
  }

  hat <- if(options$hat) {
    if(!exists("h_Pz") || !exists("h_X")) {
      h_Pz <- rowSums(qr.Q(qr_z)^2)
      Ai_X <- qr.solve(qr_a, t(X))
      h_X <- vapply(idx, function(i) {X[i, ] %*% Ai_X[, i]}, numeric(1L))
    }
    cbind("stage1" = h_Pz, "stage2" = rowSums(qr.Q(qr_x)^2), "projection" = h_X)
  } else {
    matrix(0, N, 3, dimnames = list(NULL, c("stage1", "stage2", "projection")))
  }

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
  } else { # BGM's derivative
    sigma_i <- matrix(res^2 / (N - K) + sigma, N)
  }

  # Standard errors
  if(isTRUE(options$se) && (is.null(cluster) || isFALSE(options$cluster))) {
    ZX <- crossprod(Z, X)
    ZZ_inv <- chol2inv(qr.R(qr_z))
    se_i <- vapply(idx, function(i) {
      inv_i <- tryCatch(update_inv(ZZ_inv, Z[i, , drop = FALSE]),
        error = function(e) {matrix(NaN, nrow(X), ncol(X))})
      proj_i <- Z[-i, ] %*% inv_i %*%
        update_cp(ZX, Z[i, , drop = FALSE], X[i, , drop = FALSE])
      sqrt(diag(chol2inv(chol(crossprod(proj_i)))) * sigma_i[i]^2)
    }, numeric(K))
    if(K != 1) {se_i <- t(se_i)} else {se_i <- matrix(se_i)}
  } else if(isTRUE(options$se) && !is.null(cluster)) {
    ZX <- crossprod(Z, X)
    ZZ_inv <- chol2inv(qr.R(qr_z))
    se_i <- vapply(idx, function(i) {
      res_i <- res[-i] + as.numeric(X[-i, ] %*% t(beta - beta_i[i, , drop = FALSE]))
      inv_i <- tryCatch(update_inv(ZZ_inv, Z[i, , drop = FALSE]),
        error = function(e) {matrix(NaN, nrow(X), ncol(X))})
      proj_i <- Z[-i, ] %*% inv_i %*%
        update_cp(ZX, Z[i, , drop = FALSE], X[i, , drop = FALSE])
      veggies_i <- veggiesCL(res_i, proj_i,
        cluster = cluster[-i, , drop = FALSE], type = "HC0")
      bread_i <- chol2inv(chol(crossprod(proj_i))) * (N - 1)
      sqrt(diag(1 / (N - 1) * (bread_i %*% veggies_i %*% bread_i)))
    }, numeric(K))
    if(K != 1) {se_i <- t(se_i)} else {se_i <- matrix(se_i)}
  } else { # BGM's derivative
    sigma_dbeta <- -2 * colSums(res * X) / (N - K)
    R_x <- if(!is.null(R)) {solve_cholesky(R, t(X))} else {XX_inv %*% t(X)}
    R_z <- solve(qr(crossprod(Z, X)), t(Z))
    veggies_i <- (R_z^2 - 2 * R_x * R_z) * sigma^2 +
      outer(diag(XX_inv), (sigma_i[, 1L] - sigma))
    sigma_d <- colSums(sigma_dbeta * -(t(beta_i) - beta)) +
      (sigma_i[, 1L] - sigma)
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

  structure(list(
    "model" = model,
    "hat" = hat,
    "beta_i" = beta_i, "sigma_i" = sigma_i, "se_i" = se_i, "tstat_i" = tstat_i,
    "cooksd" = cooksd, "dffits" = dffits,
    "rstudent" = rstudent, "covratio" = covratio,
    "meta" = meta
  ), class = "influence")
}
