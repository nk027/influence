

# Simulations -----

devtools::load_all()


# Simulate some data ---
set.seed(42)
N <- 100L
N_a <- N * .02
N_b <- N * .03
M <- 1L
n_save <- 10L
n_draw <- 1000L

# Track sensitvity and coefficients
a0 <- a1 <- a2 <- b0 <- matrix(NA_real_, n_draw, n_save + 1L)
coef_f <- coef_s1 <- coef_s2 <- coef_m <- vector("numeric", 1000L)
for(i in seq(n_draw)) {
x <- c(
  rt(N_a, df = 8) + 7,
  rt(N_b, df = 8) + 5,
  rt((N - N_a - N_b), df = 8))
y <- c(
  x[seq(N_a)] * 3 + rnorm(N_a),
  x[seq(N_a + 1, N_a + N_b)] * 2 + rnorm(N_b),
  x[seq(N_a + N_b + 1, N)] * 1 + rnorm(N - N_a - N_b))
(m <- lm(y ~ x - 1))
coef_f[i] <- coef(m)
coef_s1[i] <- coef(lm(y[-seq(N_a)] ~ x[-seq(N_a)] - 1))
coef_s2[i] <- coef(lm(y[-seq(N_a + N_b)] ~ x[-seq(N_a + N_b)] - 1))
coef_m[i] <- coef(MASS::rlm(y ~ x - 1))
out1 <- init.default(m)
out2 <- sens.lm(m, options = set_options(n_max = n_save, adaptive = FALSE), verbose = FALSE)
out3 <- sens.lm(m, options = set_options(n_max = n_save), verbose = FALSE)
out4 <- init.default(m, options = set_compute("none"))
a0[i, ] <- out1$initial[seq(n_save + 1L)]
a1[i, ] <- c(coef(m)[1L], out2$influence$lambda[seq(n_save)])
a2[i, ] <- c(coef(m)[1L], out3$influence$lambda[seq(n_save)])
b0[i, ] <- out4$initial[seq(n_save + 1L)]
}

# Interestingly, robust estimation kicks the first set
cbind(summary(coef_f), summary(coef_s1), summary(coef_s2), summary(coef_m))

# Plot one instance of the data ---
cairo_pdf("paper/output/sim_influential-set.pdf", height = 4.2, width = 4, pointsize = 16)
op <- par(mar = c(2, 2.75, 0, .5), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(x) + c(-.1, .1), ylim = range(y) + c(-.1, .1))
points(x, y, cex = 1.25)
axis(1, at = round(c(min(x), 0, max(x)), 1))
axis(2, at = round(c(min(y), 0, max(y)), 1), las = 1)
# Add regression line
abline(lm(y ~ x - 1), col = "darkgray", lty = 1, lwd = 2)
abline(lm(y[-seq(N_a + N_b)] ~ x[-seq(N_a + N_b)] - 1), col = "black", lty = 2, lwd = 2)
# Highlight sets
points(x[seq(N_a + N_b)], y[seq(N_a + N_b)], cex = 1.5, col = "darkgray",
pch = 21, bg = (viridisLite::inferno(n = 5, begin = 0.5)))
points(x[seq(N_a)], y[seq(N_a)], cex = 1.75, pch = 4)
points(x[seq(N_a + 1, N_a + N_b)], y[seq(N_a + 1, N_a + N_b)], cex = 1.75, pch = 3)
dev.off()

# Plot simulation results ---
cairo_pdf("paper/output/sim_results.pdf", height = 2.4, width = 8, pointsize = 12)
# png("paper/output/sim_results.png", height = 240, width = 800, pointsize = 24)
op <- par(mfrow = c(1, 4), mar = c(2, .5, .5, 0), family = "Noto Sans")
plot.new()
plot.window(xlim = c(-.1, 10.1), ylim = c(0.5, 2.2))
axis(1, at = c(0, 2, 5, 10))
axis(2, at = c(2.2, 2, 1.8, 1.5, 1, .5, 0), las = 1)
for(i in seq_len(nrow(a0))) {
lines(0:10, b0[i, ], col = "#66666622")
}
abline(h = c(2, 1.5, 1), lty = 3, col = "darkgray")
abline(v = c(2, 5), lty = 3, col = "darkgray")
lines(c(0.1, 1:9, 9.9), colMeans(a2), col = "#00C0C0", lwd = 3, lty = 1)
lines(c(0.1, 1:9, 9.9), colMeans(a1), col = "#C000C0", lwd = 3, lty = 2)
lines(c(0.1, 1:9, 9.9), colMeans(a0), col = "#60C060", lwd = 3, lty = 1)
lines(c(0.1, 1:9, 9.9), colMeans(b0), col = "#C0C0C0", lwd = 3, lty = 2)
abline(v = c(0, 10), lty = 1, lwd = 1, col = "black")
title(main = "B0", line = -1.5, adj = .8,
  font.main = 2, cex.main = 1.75, family = "Merriweather")
plot.new()
plot.window(xlim = c(-.1, 10.1), ylim = c(0.5, 2.2))
axis(1, at = c(0, 2, 5, 10))
# axis(2, at = c(2.2, 2, 1.8, 1.5, 1, .5, 0), las = 1)
for(i in seq_len(nrow(a0))) {
lines(0:10, a0[i, ], col = "#66666622")
}
abline(h = c(2, 1.5, 1), lty = 3, col = "darkgray")
abline(v = c(2, 5), lty = 3, col = "darkgray")
lines(c(0.1, 1:9, 9.9), colMeans(a2), col = "#00C0C0", lwd = 3, lty = 1)
lines(c(0.1, 1:9, 9.9), colMeans(a1), col = "#C000C0", lwd = 3, lty = 2)
lines(c(0.1, 1:9, 9.9), colMeans(a0), col = "#60C060", lwd = 3, lty = 1)
lines(c(0.1, 1:9, 9.9), colMeans(b0), col = "#C0C0C0", lwd = 3, lty = 2)
abline(v = c(0, 10), lty = 1, lwd = 1, col = "black")
title(main = "A0", line = -1.5, adj = .8,
  font.main = 2, cex.main = 1.75, family = "Merriweather")
plot.new()
plot.window(xlim = c(-.1, 10.1), ylim = c(0.5, 2.2))
axis(1, at = c(0, 2, 5, 10))
# axis(2, at = c(2.2, 2, 1.8, 1.5, 1, .5, 0), las = 1)
for(i in seq_len(nrow(a0))) {
lines(0:10, a1[i, ], col = "#66666622")
}
abline(h = c(2, 1.5, 1), lty = 3, col = "darkgray")
abline(v = c(2, 5), lty = 3, col = "darkgray")
lines(c(0.1, 1:9, 9.9), colMeans(a2), col = "#00C0C0", lwd = 3, lty = 1)
lines(c(0.1, 1:9, 9.9), colMeans(a1), col = "#C000C0", lwd = 3, lty = 2)
lines(c(0.1, 1:9, 9.9), colMeans(a0), col = "#60C060", lwd = 3, lty = 1)
lines(c(0.1, 1:9, 9.9), colMeans(b0), col = "#C0C0C0", lwd = 3, lty = 2)
abline(v = c(0, 10), lty = 1, lwd = 1, col = "black")
title(main = "A1", line = -1.5, adj = .8,
  font.main = 2, cex.main = 1.75, family = "Merriweather")
plot.new()
plot.window(xlim = c(-.1, 10.1), ylim = c(0.5, 2.2))
axis(1, at = c(0, 2, 5, 10))
# axis(2, at = c(2.2, 2, 1.8, 1.5, 1, .5, 0), las = 1)
for(i in seq_len(nrow(a0))) {
lines(0:10, a2[i, ], col = "#66666622")
}
abline(h = c(2, 1.5, 1), lty = 3, col = "darkgray")
abline(v = c(2, 5), lty = 3, col = "darkgray")
lines(c(0.1, 1:9, 9.9), colMeans(a2), col = "#00C0C0", lwd = 3, lty = 1)
lines(c(0.1, 1:9, 9.9), colMeans(a1), col = "#C000C0", lwd = 3, lty = 2)
lines(c(0.1, 1:9, 9.9), colMeans(a0), col = "#60C060", lwd = 3, lty = 1)
lines(c(0.1, 1:9, 9.9), colMeans(b0), col = "#C0C0C0", lwd = 3, lty = 2)
abline(v = c(0, 10), lty = 1, lwd = 1, col = "black")
title(main = "A2", line = -1.5, adj = .8,
  font.main = 2, cex.main = 1.75, family = "Merriweather")
dev.off()
# Inkscape magic to fix the axis


# 2SLS versus OLS ---
set.seed(42)
n_draw <- 1000L
N <- 100L
s_lm <- s_iv <- vector("list", n_draw)

for(i in seq(n_draw)) {
  w <- rt(N, 8)
  x <- w + rt(N, 8)
  y <- x + rt(N, 8)
  s_lm[[i]] <- sens.lm(lm(y ~ x), lambda = set_lambda("beta", pos = 2),
    options = set_options(hat = TRUE, beta = TRUE, n_max = 95))
  s_iv[[i]] <- sens.ivreg(ivreg::ivreg(y ~ x | w),
    lambda = set_lambda("beta", pos = 2),
    options = set_options(hat = TRUE, beta = TRUE, n_max = 95))
}

z_lm <- sapply(s_lm, \(s) {x <- summary.sensitivity(s)$exact; x[grepl("zero[1]*$", names(x))]})
z_iv <- sapply(s_iv, \(s) {x <- summary.sensitivity(s)$exact; x[grepl("zero[1]*$", names(x))]})

len_lm <- sapply(s_lm, \(s) length(s$influence$lambda))
len_iv <- sapply(s_iv, \(s) length(s$influence$lambda))

# Highlight a few
mean_lm <- sapply(s_lm, \(s) s$influence$lambda) |> apply(1, mean, na.rm = TRUE)
mean_iv <- sapply(s_iv, \(s) {
  v <- rep(NA_real_, 96)
  v[seq(length(s$influence$lambda))] <- s$influence$lambda
  v
}) |> apply(1, mean, na.rm = TRUE)
med_lm <- sapply(s_lm, \(s) s$influence$lambda) |>
apply(1, quantile, c(.05, .5, .95), na.rm = TRUE)
med_iv <- sapply(s_iv, \(s) {
  v <- rep(NA_real_, 96)
  v[seq(length(s$influence$lambda))] <- s$influence$lambda
  v
}) |> apply(1, quantile, c(.05, .5, .95), na.rm = TRUE)


cairo_pdf("paper/output/sim_2sls.pdf", height = 5, width = 8, pointsize = 8)
op <- par(mfrow = c(2, 1), mar = c(3, 3, .5, .5), bg = "white", family = "Noto Sans")
# OLS on top
plot.new()
plot.window(xlim = c(100, 5), ylim = c(-20, 2.5))
axis(1, at = c(100, 80, 50, 20, 5))
axis(2, at = c(1, 0, -1, -10, -20))
abline(h = c(1, 0, -1, -10, -20), col = "darkgray", lty = 3)
abline(v = c(100, 80, 50, 20), col = "darkgray", lty = 3)
for(i in sample(seq_along(s_lm), 250))
  lines(x = s_lm[[i]]$influence$N, y = s_lm[[i]]$influence$lambda, col = "#66666622")
lines(x = 100:5, y = med_lm[1, ], col = "#00C0E0", lwd = 1.5, lty = 2)
lines(x = 100:5, y = med_lm[2, ], col = "#00C0E0", lwd = 2)
lines(x = 100:5, y = med_lm[3, ], col = "#00C0E0", lwd = 1.5, lty = 2)
lines(x = 100:5, y = mean_lm, col = "#E000C0", lwd = 2, lty = 3)
abline(v = 5, col = "black", lty = 1)
# title(xlab = "Observations", line = 1, adj = .35, font.lab = 2, cex.lab = 1.25)
title(main = "OLS", line = -1, adj = .1, font.main = 2, cex.main = 1.75, family = "Merriweather")
# 2SLS next
plot.new()
plot.window(xlim = c(100, 5), ylim = c(-20, 2.5))
axis(1, at = c(100, 80, 50, 20, 5))
axis(2, at = c(1, 0, -1, -10, -20))
abline(h = c(1, 0, -1, -10, -20), col = "darkgray", lty = 3)
abline(v = c(100, 80, 50, 20), col = "darkgray", lty = 3)
for(i in sample(seq_along(s_iv), 500))
  lines(x = s_iv[[i]]$influence$N, y = s_iv[[i]]$influence$lambda, col = "#66666622")
lines(x = 100:5, y = med_iv[1, ], col = "#00C0E0", lwd = 1.5, lty = 2)
lines(x = 100:5, y = med_iv[2, ], col = "#00C0E0", lwd = 2)
lines(x = 100:5, y = med_iv[3, ], col = "#00C0E0", lwd = 1.5, lty = 2)
lines(x = 100:5, y = mean_iv, col = "#E000C0", lwd = 2, lty = 3)
for(i in which(len_iv != 96))
points(x = tail(s_iv[[i]]$influence$N, 1), y = 2.1 + rnorm(1L, 0, .2),
  pch = 4, col = "#60A000", lwd = 1.5, cex = 1.5)
abline(v = 5, col = "black", lty = 1)
title(xlab = "Observations", line = 1, adj = .35, font.lab = 2, cex.lab = 1.25)
title(main = "2SLS", line = -1, adj = .1, font.main = 2, cex.main = 1.75, family = "Merriweather")
dev.off()