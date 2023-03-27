
# Masking example -----

set.seed(753)
devtools::load_all()

# Simulate the data
N <- 54
x <- c(rnorm(N), rnorm(3, 6, 0.25), rnorm(3, 8, 0.25))
y <- c(
  x[seq(N)] * -0.5 + rnorm(N, 0, 1),
  x[seq(N + 1, N + 3)] * 0.1 + rnorm(3, 0, 0.1),
  x[seq(N + 4, N + 6)] * 0.4 + rnorm(3, 0, 0.1)
)

# Fit a model and briefly check it out
m <- lm(y ~ x - 1)
summary(m)
plot(x, y)
abline(m, col = "darkgray", lwd = 2)

# Re-run with removals
m_ab <- lm(y[1:54] ~ x[1:54] - 1)
m_a <- lm(y[1:57] ~ x[1:57] - 1)
abline(m_a, col = "darkgray", lty = 2)
abline(m_ab, col = "darkgray", lty = 3)

summary(m_a)
summary(m_ab)


# Plot the data ---

plot_base <- \(xax = FALSE, yax = FALSE) {
  plot.new()
  plot.window(xlim = range(x) + c(-.12, .12), ylim = range(y) + c(-.1, .1))
  if(xax) axis(1, at = c(min(x), 0, max(x)), labels = FALSE, tck = -.02)
  if(yax) axis(2, at = c(min(y), 0, max(y)), labels = FALSE, tck = -.02)
}
plot_rms <- \(rm, offset = .1, xpos = c(8.2, 8.1, 8.0),
  col = viridisLite::inferno(n = 3, begin = 0.2, end = 0.8)) {
  cf <- vector("numeric", length(rm) + 1)
  cf[1] <- coef(lm(y ~ x - 1))
  for(i in seq_along(rm)) {
    cf[i + 1] <- coef(lm(y[-rm[[i]]] ~ x[-rm[[i]]] - 1))
    arrows(x0 = xpos[i], x1 = xpos[i], col = col[i],
      y0 = xpos[i] * cf[i] - offset, y1 = xpos[i] * cf[i + 1] + offset, length = 0.05)
    segments(2, 2 * cf[i + 1], xpos[i], xpos[i] * cf[i + 1],
      col = "black", lwd = .5, lty = 3)
    text(xpos[i] + .4, y = xpos[i] * cf[i + 1] + offset, label = paste0(i, "."))
  }
}

# cairo_pdf("paper/output/masking/masked_data.pdf", height = 4.2, width = 4)
cairo_pdf("paper/output/masking/masked_data.pdf", height = 5, width = 7, pointsize = 12)
op <- par(mar = c(.8, .5, 0, .5), fig = c(0, .5, .49, 1),
  bg = "transparent", family = "Noto Sans")
plot_base(xax = TRUE, yax = TRUE)
points(x, y, cex = 1.25)
# Add regression line
abline(m, col = "#333333", lty = 1, lwd = 2)
text(x = -1.8, y = 3.3, adj = 0, labels = "Influential sets",
  col = "black", cex = 1.5, font = 2, family = "Merriweather")

par(fig = c(.5, 1, .49, 1), new = TRUE)
plot_base(xax = TRUE, yax = TRUE)
points(x, y, cex = 1.25)
points(x[58:60], y[58:60], cex = 1.25, pch = 21, col = "darkgray",
  bg = (viridisLite::inferno(n = 3, begin = 0.2, end = 0.8)))
# Add regression line
abline(m, col = "#333333", lty = 1, lwd = 2)
plot_rms(list(60, 60:59, 60:58))
text(x = -1.8, y = 3.3, adj = 0, labels = "Joint influence",
  col = "black", cex = 1.2, font = 2, family = "Merriweather")

par(fig = c(0, .5, 0, .52), new = TRUE, mar = c(.5, .5, .5, .5))
plot_base(xax = TRUE, yax = TRUE)
# Add regression line
abline(m, col = "darkgray", lty = 3, lwd = 2)
points(x[seq(N + 3)], y[seq(N + 3)], cex = 1.25)
points(x[seq(N + 4, N + 6)], y[seq(N + 4, N + 6)], cex = 1.25,
  pch = 4, col = "gray")
symbols(add = TRUE, mean(x[seq(N + 4, N + 6)]), mean(y[seq(N + 4, N + 6)]),
  circle = .8, inches = FALSE, fg = "#800080", lwd = 5)
text(x = 6.6, y = mean(y[seq(N + 4, N + 6)]) - .5, labels = "a", cex = 1.5, font = 2)
abline(m_a, col = "#333333", lty = 1, lwd = 2)

par(fig = c(.5, 1, 0, .52), new = TRUE)
plot_base(xax = TRUE, yax = TRUE)
# Add regression line
abline(m, col = "darkgray", lty = 3, lwd = 2)
abline(m_a, col = "darkgray", lty = 3, lwd = 2)
points(x[seq(N)], y[seq(N)], cex = 1.25)
points(x[seq(N + 1, N + 3)], y[seq(N + 1, N + 3)], cex = 1.25,
  pch = 4, col = "gray")
abline(m_ab, col = "#333333", lty = 1, lwd = 2)
symbols(add = TRUE, mean(x[seq(N + 1, N + 3)]), mean(y[seq(N + 1, N + 3)]),
  circle = .8, inches = FALSE, fg = "#800080", lwd = 5)
text(x = -1.8, y = 3.3, adj = 0, labels = "Masking",
  col = "black", cex = 1.2, font = 2, family = "Merriweather")
text(x = 7.4, y = mean(y[seq(N + 1, N + 3)]) - .5, labels = "b", cex = 1.5, font = 2)
abline(m_ab, col = "#333333", lty = 1, lwd = 2)
dev.off()


# Applied masking example ---

# Influence calculation
a0 <- init.default(m)
id0 <- a0$id[1:7]
id2 <- sens.lm(m)$influence$id[1:7]

# Panel 1
cairo_pdf("paper/output/masking/masked_alg0.pdf", height = 4.2, width = 4)
op <- par(mar = c(.8, .8, 0, .5), bg = "transparent", family = "Noto Sans")
plot_base(xax = TRUE, yax = TRUE)
text(2.6, 3.4, # Title
  labels = expression("Algorithms " * phantom("0") * " & " * phantom("1")),
  col = "black", cex = 1.5, font = 2, family = "Merriweather")
text(2.6, 3.4, labels = expression(phantom("Algorithms ") * "0" * phantom(" & 1")),
  col = "#800080", cex = 1.5, font = 2, family = "Merriweather")
text(2.6, 3.4, labels = expression(phantom("Algorithms 0 & ") * "1"),
  col = "#408040", cex = 1.5, font = 2, family = "Merriweather")
# Add regression line
abline(m, col = "darkgray", lty = 3, lwd = 2)
points(x, y, cex = 1.25)
# Highlight the set
points(x[id0], y[id0], cex = 1.5, col = "darkgray",
  pch = 21, bg = (viridisLite::inferno(n = length(id0), begin = 0.2)))
points(x[id0[seq(3)]], y[id0[seq(3)]], cex = 1.75, pch = 4)
points(x[id0[-seq(3)]], y[id0[-seq(3)]], cex = 1.75, pch = 3)
# Add algorithms
for(i in seq(2, 3))
  segments(2, 2 * a0$initial[i], 8, 8 * a0$initial[i],
    col = "#800080", lwd = .1, lty = 3)
abline(lm(y[-id0[1:3]] ~ x[-id0[1:3]] - 1), col = "#408040", lty = 1, lwd = 2)
abline(lm(y[-id0] ~ x[-id0] - 1), col = "#408040", lty = 2, lwd = 2)
abline(0, a0$initial[4], col = "#800080", lty = 1, lwd = 2)
abline(0, a0$initial[8], col = "#800080", lty = 2, lwd = 2)
# Add a customised legend (cut for the PDF)
rect(xleft = 2.4, ybottom = -3.5, xright = 8.05, ytop = -2.6, col = "white")
legend(1.55, -2.82 + .05, horiz = TRUE, bty = "n",
  fill = (viridisLite::inferno(n = length(id0), begin = .2)), cex = 2,
  border = NA, x.intersp = -1, y.intersp = 0.4, legend = rep(NA, length(id0)))
legend(1.35 + .15, -2.87, horiz = TRUE, bty = "n", pch = c(4, 4, 4, 3, rep(3, 3)),
  col = c("lightgray", (viridisLite::inferno(n = length(id0), begin = .2))[2:3],
    "black", (viridisLite::inferno(n = length(id0), begin = .2))[5:7]),
  fill = rep("transparent", 7), cex = 1.6,
  border = NA, x.intersp = -0.8, y.intersp = 0.4, legend = rep(NA, length(id0)))
text(5.05, -2.80, labels = "⟶ set order ⟶", cex = 1.2)
axis(1, at = c(min(x), 0, max(x)), labels = FALSE, tck = -.03)
dev.off()

# Panel 2
cairo_pdf("paper/output/masking/masked_alg2.pdf", height = 4.2, width = 4)
op <- par(mar = c(.8, .8, 0, .5), bg = "transparent", family = "Noto Sans")
plot_base(xax = TRUE, yax = TRUE)
text(2.6, 3.4, labels = expression("Algorithm " * phantom("2")), # Title
  col = "black", cex = 1.5, font = 1, family = "Merriweather")
text(2.6, 3.4, labels = expression(phantom("Algorithm ") * "2"),
  col = "#008080", cex = 1.5, font = 1, family = "Merriweather")
# Add regression line
abline(m, col = "darkgray", lty = 3, lwd = 2)
points(x, y, cex = 1.25)
# Highlight the set
points(x[id2], y[id2], cex = 1.5, col = "darkgray",
  pch = 21, bg = (viridisLite::inferno(n = length(id2), begin = 0.2)))
points(x[id2[seq(3)]], y[id2[seq(3)]], cex = 1.75, pch = 4)
points(x[id2[-seq(3)]], y[id2[-seq(3)]], cex = 1.75, pch = 3)
# Add algorithms
for(i in c(1, 2, 4, 5, 6)) {
  cf <- coef(lm(y[-id2[seq(i)]] ~ x[-id2[seq(i)]] - 1))
  segments(2, 2 * cf, 8, 8 * cf, col = "#008080", lwd = .1, lty = 3)
}
abline(lm(y[-id2[1:3]] ~ x[-id2[1:3]] - 1), col = "#008080", lty = 1, lwd = 2)
abline(lm(y[-id2] ~ x[-id2] - 1), col = "#008080", lty = 2, lwd = 2)
# Add labels for the number of removals
text(7.2, -.6, labels = "3 removed", col = "black", cex = 1.2, font = 1)
text(3.2, -3.2, labels = "7", col = "black", cex = 1.2, font = 1)
dev.off()
