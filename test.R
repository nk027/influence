
sapply(list.files("R/", ".R$"), function(x) source(paste0("R/", x)))

N <- 1000
K <- 30
M <- 30
Z <- matrix(rnorm(N * M), N)
X <- Z %*% (b0 <- matrix(rnorm(M * K), M)) + (e0 <- rnorm(N * K))
y <- X %*% (b1 <- rnorm(K)) + (e1 <- X[, 1] * 2 * (u1 <- rnorm(N)))
