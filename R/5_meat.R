
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
