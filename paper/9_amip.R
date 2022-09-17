
# AMIP performance -----

set.seed(27)

# devtools::install_github("rgiordan/zaminfluence", subdir = "zaminfluence")
library("zaminfluence") # v0.0.0.9000 / commit 43b53f7 and recently 2d29fbb
devtools::load_all()

# Leverage & error make up influence on coefficients
h <- \(x) diag(x %*% solve(crossprod(x), t(x)))


# We start by demonstrating the lack of leverage ---
N <- 1000 # Simulate data for a univariate regression
x <- rt(N, df = 8) # t for some leverage
y <- x + rt(N, df = 8) # Normal errors

# Fit a linear model
m <- lm(y ~ x - 1, x = TRUE, y = TRUE)
summary(m)
# plot(m)

# Show the data and regression
cairo_pdf("paper/output/amip_sim.pdf", height = 4.2, width = 4)
op <- par(mar = c(2, 3, 0, .5), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(x) + c(-.1, .1), ylim = range(y) + c(-.1, .1))
points(x, y, cex = 1.25)
axis(1, at = round(c(min(x), max(x)), 2))
axis(2, at = round(c(min(y), 0, max(y)), 2), las = 1)
# Add regression line
abline(m, col = "darkgray", lty = 1, lwd = 2)
dev.off()

# Compute the AMIP and its gradients
amip <- ComputeModelInfluence(m) |> AppendTargetRegressorInfluence("x")
# Compute our equivalent
infl_app <- infl.lm(m, options = set_compute("none"))

# Coefficient influence ---
# Check equivalence to (X'X)⁻¹ x e
dfbad <- t(solve(crossprod(x), t(x * resid(m))))
norm(t(amip$param_grad) - dfbad) < 1e-12
norm(t(-amip$param_grad) - (infl_app$beta_i - infl_app$model$beta)) < 1e-12

# Compare against the exact values
dfbex <- dfbeta(m)

# Plot errors vs. leverage and influence
y_values <- scale(dfbex - dfbad)
ylim <- range(y_values) + c(-.5, .5)
# First part
cairo_pdf("paper/output/amip_b-i.pdf", height = 4, width = 4)
op <- par(mar = c(4, 4, 1, 1), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(dfbex) + c(-.0005, .0005), ylim = ylim)
points(dfbex, y_values, cex = 1.25, pch = 4)
axis(1, at = round(c(min(dfbex), 0, max(dfbex)), 3))
axis(2, at = round(c(min(y_values), 0,
  max(y_values)), 1), las = 1)
abline(v = 0, col = "darkgray", lty = 2)
abline(h = 0, col = "darkgray", lty = 2)
title(ylab = "Standardized error", xlab = "Influence", line = 2.5, font.lab = 2)
dev.off()
# Second part
cairo_pdf("paper/output/amip_b-h.pdf", height = 4, width = 4)
op <- par(mar = c(4, 4, 1, 1), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(h(x)) + c(-1e-3, 5e-3), ylim = ylim)
points(h(dfbex), y_values, cex = 1.25, pch = 4)
axis(1, at = round(c(0, max(h(x))), 3))
axis(2, at = round(c(min(y_values), 0,
  max(y_values)), 1), las = 1)
abline(h = 0, col = "darkgray", lty = 2)
title(ylab = "Standardized error", xlab = "Leverage", line = 2.5, font.lab = 2)
dev.off()


# Error influence ---
# Compute exact values
infl <- infl.lm(m, options = set_compute(se = TRUE))
se_i <- infl$se_i[, 1] - infl$model$se

# Pull BGM values
se_bgm <- -amip$se_grad[1, ]
# Check equivalence to our implementation
all(abs(se_bgm - (infl_app$se_i[, 1] - infl$model$se)) < 1e-12)

# We get a good fit, on average we're on point
m_se <- lm(se_i ~ se_bgm - 1)
summary(m_se)
# plot(m_se)

# Plot the fit
s_res <- scale(resid(m_se))
# First part
cairo_pdf("paper/output/amip_s-f.pdf", height = 4, width = 4)
op <- par(mar = c(4, 4, 1, 1), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(fitted(m_se)) + c(-1e-4, 1e-4), ylim = range(s_res) + c(-.1, .1))
points(x = fitted(m_se), y = s_res, pch = 3, cex = 1.25)
axis(1, at = round(c(min(se_i), 0, max(se_i)), 4))
axis(2, at = round(c(min(s_res), 0, max(s_res)), 1), las = 1)
abline(v = 0, col = "darkgray", lty = 2)
abline(h = 0, col = "darkgray", lty = 2)
title(ylab = "Standardized residuals", xlab = "Fitted values", line = 2.5, font.lab = 2)
dev.off()
# Second part
cairo_pdf("paper/output/amip_s-h.pdf", height = 4, width = 4)
op <- par(mar = c(4, 4, 1, 1), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(h(se_i)) + c(-1e-3, 1e-3), ylim = range(s_res) + c(-.1, .1))
points(x = h(se_i), y = s_res, pch = 3, cex = 1.25)
axis(1, at = round(c(0, max(h(se_i))), 4))
axis(2, at = round(c(min(s_res), 0, max(s_res)), 1), las = 1)
abline(v = 0, col = "darkgray", lty = 2)
abline(h = 0, col = "darkgray", lty = 2)
title(ylab = "Standardized residuals", xlab = "Leverage", line = 2.5, font.lab = 2)
dev.off()


# Consider t values ---
# This is their approximation, see microcredit application
amip_s <- amip |> GetInferenceSignals()
# Influence on the significance and a significant signflip
all(abs(amip_s$x$sig$qoi$infl -
  (amip$param_grad[1, ] - qnorm(.975) * amip$se_grad[1, ])) < 1e-12)
all(abs(amip_s$x$both$qoi$infl -
  (amip$param_grad[1, ] + qnorm(.975) * amip$se_grad[1, ])) < 1e-12)
# Significance is flagged beyond the base value
abs(amip_s$x$sig$qoi$base_value -
  (infl$model$beta - qnorm(.975) * infl$model$se)) < 1e-12
abs(amip_s$x$both$qoi$base_value -
  (infl$model$beta + qnorm(.975) * infl$model$se)) < 1e-12


# Monte Carlo simulation to check order and underestimation ---
N <- 100
n_draw <- 100000

s_order1 <- s_order5 <- s_abs <- s_rel <- vector("numeric", n_draw)
s_bias <- matrix(NA_real_, n_draw, 3)
for(i in seq(n_draw)) {
  # Simulate data for a univariate regression
  x <- rt(N, df = 8)
  y <- x + rnorm(N)
  m <- lm(y ~ x - 1, x = TRUE, y = TRUE)
  # AMIP estimate and exact one
  dfbad <- t(solve(crossprod(x), t(x * resid(m))))
  dfbex <- dfbeta(m)
  o_dfbad <- order(dfbad)
  o_dfbex <- order(dfbex)
  # Compare the top and top 5 observations identified
  s_order1[i] <- sum(o_dfbad[seq(N * .01)] != o_dfbex[seq(N * .01)]) / (N * .01)
  s_order5[i] <- sum(o_dfbad[seq(N * .05)] != o_dfbex[seq(N * .05)]) / (N * .05)
  # Bias (induced by the lack of leverage)
  bias <- (abs(dfbad) / abs(dfbex))
  s_bias[i, ] <- c(min(bias), mean(bias), max(bias))
  # Difference between AMIP and the influence
  diff <- abs(sum(dfbex[o_dfbex[seq(N * .01)]]) - sum(dfbad[o_dfbad[seq(N * .01)]]))
  s_abs[i] <- diff
  s_rel[i] <- diff / abs(sum(dfbex[o_dfbex[seq(N * .01)]]))
}

summary(s_abs) # These need some context
summary(s_rel) # Usually, we are 5.4% off, at most 68.5%
summary(s_order1) # We miss the most influential one 4% of the time
summary(s_order5) # We miss the five most influential ones 10.3% of the time
summary(s_bias) # The worst estimate is on average at 88.8% of the true value

# Summarise in a table
f <- \(x) c(mean(x), min(x), max(x), quantile(x, .1), median(x), quantile(x, .9))
tbl <- as.data.frame(t(data.frame(
  "most influential, error / influence" = round(f(s_rel), 4),
  "worst AMIP / influence" = round(f(s_bias[, 1]), 4)
  "average AMIP / influence" = round(f(s_bias[, 2]), 4),
  "most influential missed" = round(f(s_order1), 4),
  "top 5 missed" = round(f(s_order5), 4),
)))
colnames(tbl) <- c("mean", "minimum", "maximum", "1st decile", "median", "9th decile")

memisc:::toLatex.data.frame(tbl, digits = 4)
