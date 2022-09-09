
# Update inverse
update_inv <- function(XX_inv, X_rm) {

  if(NROW(X_rm) == 1) {
    out <- XX_inv + (XX_inv %*% crossprod(X_rm) %*% XX_inv) /
      as.numeric(1 - X_rm %*% tcrossprod(XX_inv, X_rm))
  } else {
    out <- XX_inv + tcrossprod(XX_inv, X_rm) %*%
      solve(diag(NROW(X_rm)) - X_rm %*% tcrossprod(XX_inv, X_rm),
      X_rm %*% XX_inv)
  }
  if(abs(norm(XX_inv, "I") - norm(out, "I")) > 1e12) {
    stop("Inverse update likely to suffer from numerical inaccuracy.")
  }
  return(out)
}

# Update crossproduct
update_cp <- function(XY, X_rm, Y_rm = X_rm) {
  XY - crossprod(X_rm, Y_rm)
}


# Solve using an upper triangular matrix
solve_cholesky <- function(R, b) {
  backsolve(R, forwardsolve(R, b, upper.tri = TRUE, transpose = TRUE))
}


# Marginalise using Frisch-Waugh-Lovell theorem
update_fwl <- function(X, y, variables, rm = NULL) {
  if(!any(variables == 0)) {
    if(is.null(rm)) {
      Q_fwl <- qr.Q(qr(X[, -variables, drop = FALSE]))
      y <- y - Q_fwl %*% crossprod(Q_fwl, y)
      X <- X[, variables, drop = FALSE] - Q_fwl %*%
        crossprod(Q_fwl, X[, variables, drop = FALSE])
    } else {
      Q_fwl <- qr.Q(qr(X[-rm, -variables, drop = FALSE]))
      y[-rm] <- y[-rm] - Q_fwl %*% crossprod(Q_fwl, y[-rm])
      X[-rm, variables] <- X[-rm, variables, drop = FALSE] - Q_fwl %*%
        crossprod(Q_fwl, X[-rm, variables, drop = FALSE])
      X <- X[, variables]
    }
  }
  return(list("y" = y, "X" = X))
}


# Adaptation of sandwich::meatCL
veggiesCL <- function(residual, X,
  cluster = NULL,
  type = c("HC0", "HC1"), # ll and LM
  ...) {

  type <- match.arg(type)

  ef <- residual * X
  K <- NCOL(X)
  N <- NROW(X)

  # Allow multi-way clustering
  P <- NCOL(cluster)
  if(P > 1L) {
    cl <- unlist(lapply(seq_len(P), function(i) {
      combn(seq_len(P), i, simplify = FALSE)
    }), recursive = FALSE)
    sign <- vapply(cl, function(i) {(-1L)^(length(i) + 1L)}, numeric(1L))
    paste_ <- function(...) {paste(..., sep = "_")}
    for(i in seq.int(P + 1L, length(cl))) {
      cluster <- cbind(cluster, Reduce(paste_, unclass(cluster[, cl[[i]]])))
    }
  } else {
    cl <- list(1)
    sign <- 1
  }

  # Number of clusters and cluster interactions
  G <- sapply(seq_along(cl), function(i) {
    if(is.factor(cluster[[i]])) {
      length(levels(cluster[[i]]))
    } else {
      length(unique(cluster[[i]]))
    }
  })

  out <- matrix(0, nrow = K, ncol = K)

  for (i in seq_along(cl)) {
    adj <- G[i] / (G[i] - 1L) # Cluster adjustment
    # Aggregate within cluster levels
    efi <- if(G[i] < N) {apply(ef, 2L, rowsum, cluster[[i]])} else {ef}
    # Aggregate across cluster variables
    out <- out + sign[i] * adj * crossprod(efi) / N
  }
  # HC1 adjustment with residual degrees of freedom
  if(type == "HC1") {out <- (N - 1L) / (N - K) * out}

  return(out)
}
