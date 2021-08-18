
set.seed(42)

# Simulate some data ---
N <- 1000
K <- 10

beta <- rnorm(K)
X <- cbind(1, matrix(rnorm(N * (K - 1)), N))
Z <- cbind(1, matrix(rnorm(N * (K - 1)), N))
y <- X %*% beta + rnorm(N, 0, 1)

# Regress ---
mdl1 <- lm(y ~ X - 1)
mdl2 <- ivreg::ivreg(y ~ X - 1 | Z - 1)
mdl3 <- AER::ivreg(y ~ X - 1 | Z - 1)

# Check influential observations ---
lapply(list.files("R", "R$"), function(x) {source(paste0("R/", x))})

# Work-horse -- yields model quantities and some after removal
options <- set_compute()
infl(mdl1, options = options)
infl(mdl2, options = options)
infl(mdl3, options = options)

# Initial approximation -- requires a target quantity
lambda <- set_lambda("beta", pos = 2, sign = -1)
init(mdl1, lambda = lambda)
init(mdl2, lambda = lambda)
init(mdl3, lambda = lambda)

# Divide and conquer -- requires a target value
goal(mdl1, lambda = lambda, target = set_target(-.25, "geq"))
goal(mdl2, lambda = lambda, target = set_target(0, "geq"))
goal(mdl3, lambda = lambda, target = set_target(0, "geq"))

# Adaptive identification -- requires additional options wrt computation
options <- set_options()
sens(mdl1, lambda = lambda, options = options)
sens(mdl2, lambda = lambda, options = options)
sens(mdl3, lambda = lambda, options = options)

# Compare adaptive and initial
plot(sens(mdl1, lambda = lambda, options = options)$influence$lambda[1:100],
  type = "l", col = 5)
lines(abs(init(mdl1, lambda = lambda)$initial[2:101]))
