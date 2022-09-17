
# Masking example -----

set.seed(42)
devtools::load_all()

# Simulate the data
N <- 54
x <- c(rnorm(N), rnorm(3, 6, 0.25), rnorm(3, 8, 0.1))
y <- c(
  x[seq(N)] * -0.5 + rnorm(N, 0, 1),
  x[seq(N + 1, N + 3)] * 0.1 + rnorm(3, 0, 0.1),
  x[seq(N + 4, N + 6)] * 0.4 + rnorm(3, 0, 0.1)
)

# Fit a model and briefly check it out
m <- lm(y ~ x - 1)
plot(x, y)
abline(m, col = "darkgray", lwd = 2)
summary(m)


# Plot the data ---

cairo_pdf("paper/output/masking/masked_data.pdf", height = 4.2, width = 4)
op <- par(mar = c(2, 2.5, 0, .5), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(x) + c(-.1, .2), ylim = range(y) + c(-.1, .1))
# grid()
points(x, y, cex = 1.25)
axis(1, at = round(c(min(x), 0, max(x)), 1))
axis(2, at = round(c(min(y), 0, max(y)), 1), las = 1)
text(x = 2, y = 3.4, labels = "Masking",
  col = "black", cex = 1.5, font = 1, family = "Merriweather")
# Add regression line
abline(m, col = "darkgray", lty = 1, lwd = 2)
abline(lm(y[seq(N)] ~ x[seq(N)] - 1), col = "black", lty = 2, lwd = 2)
# Highlight sets
symbols(add = TRUE,
  c(mean(x[seq(N + 1, N + 3)]), mean(x[seq(N + 4, N + 6)]) - .1),
  c(mean(y[seq(N + 1, N + 3)]), mean(y[seq(N + 4, N + 6)])),
  circle = c(1, 1), inches = FALSE,
  fg = c("#400080", "#400080"), lwd = 2.5)
text(x = 6.4, y = 2.8, labels = "a", cex = 1.5)
text(x = 4.9, y = 0, labels = "b", cex = 1.5)
dev.off()
# There's some Inkscape voodoo afterwards


# Applied masking example ---

a0 <- init.default(m)
id0 <- a0$id[1:7]
id2 <- sens.lm(m)$influence$id[1:7]

# Panel 1
cairo_pdf("paper/output/masking/masked_alg0.pdf", height = 4.2, width = 4)
op <- par(mar = c(1, 1, 0, .5), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(x) + c(-.1, .2), ylim = range(y) + c(-.1, .1))
# grid()
points(x, y, cex = 1.25)
axis(1, at = round(c(min(x), 0, max(x)), 1), labels = FALSE)
axis(2, at = round(c(min(y), 0, max(y)), 1), las = 1, labels = FALSE)
# Add regression line
text(2.6, 3.4,
  labels = expression("Algorithms " * phantom("0") * " & " * phantom("1")),
  col = "black", cex = 1.5, font = 1, family = "Merriweather")
text(2.6, 3.4, labels = expression(phantom("Algorithms ") * "0" * phantom(" & 1")),
  col = "#800080", cex = 1.5, font = 1, family = "Merriweather")
text(2.6, 3.4, labels = expression(phantom("Algorithms 0 & ") * "1"),
  col = "#408040", cex = 1.5, font = 1, family = "Merriweather")
# Highlight the set
points(x[id0], y[id0], cex = 1.5, col = "darkgray",
  pch = 21, bg = (viridisLite::inferno(n = length(id0), begin = 0.5)))
points(x[id0[seq(3)]], y[id0[seq(3)]], cex = 1.75, pch = 4)
points(x[id0[-seq(3)]], y[id0[-seq(3)]], cex = 1.75, pch = 3)
# Add regression lines
abline(m, col = "lightgray", lty = 1, lwd = 2)
abline(lm(y[-id0[1:3]] ~ x[-id0[1:3]] - 1), col = "#408040", lty = 3, lwd = 2)
abline(lm(y[-id0] ~ x[-id0] - 1), col = "#408040", lty = 5, lwd = 2)
abline(0, a0$initial[4], col = "#800080", lty = 3, lwd = 2)
abline(0, a0$initial[8], col = "#800080", lty = 5, lwd = 2)
text(7.8, -0.05, labels = "0", col = "#800080", cex = 1.2, font = 1)
text(7, -0.8, labels = "1", col = "#408040", cex = 1.2, font = 1)
# box()
dev.off()

# Panel 2
cairo_pdf("paper/output/masking/masked_alg2.pdf", height = 4.2, width = 4)
op <- par(mar = c(1, 1, 0, .5), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(x) + c(-.1, .2), ylim = range(y) + c(-.1, .1))
# grid()
points(x, y, cex = 1.25)
axis(1, at = round(c(min(x), 0, max(x)), 1), labels = FALSE)
axis(2, at = round(c(min(y), 0, max(y)), 1), las = 1, labels = FALSE)
# plot(x, y, cex = 1.25)
text(2.6, 3.4, labels = expression("Algorithm " * phantom("2")),
  col = "black", cex = 1.5, font = 1, family = "Merriweather")
text(2.6, 3.4, labels = expression(phantom("Algorithm ") * "2"),
  col = "#008080", cex = 1.5, font = 1, family = "Merriweather")
# Highlight the set
points(x[id2], y[id2], cex = 1.5, col = "darkgray",
  pch = 21, bg = (viridisLite::inferno(n = length(id2), begin = 0.5)))
points(x[id2[seq(3)]], y[id2[seq(3)]], cex = 1.75, pch = 4)
points(x[id2[-seq(3)]], y[id2[-seq(3)]], cex = 1.75, pch = 3)
# Add regression lines
abline(m, col = "lightgray", lty = 1, lwd = 2)
abline(lm(y[-id2[1:3]] ~ x[-id2[1:3]] - 1), col = "#008080", lty = 3, lwd = 2)
abline(lm(y[-id2] ~ x[-id2] - 1), col = "#008080", lty = 5, lwd = 2)
# box()
# Add labels for the number of removals
text(6.5, -0.5, labels = "3 removed", col = "black", cex = 1.2, font = 1)
text(4.4, -1.3, labels = "7 removed", col = "black", cex = 1.2, font = 1)
# Add a customised legend (cut for the PDF)
rect(xleft = 3.05, ybottom = -2.475, xright = 8.65, ytop = -1.675, col = "white")
legend(2.2, -1.825, horiz = TRUE, bty = "n",
  fill = (viridisLite::inferno(n = length(id0), begin = .5)), cex = 2,
  border = NA, x.intersp = -1, y.intersp = 0.4, legend = rep(NA, length(id0)))
legend(2, -1.88, horiz = TRUE, bty = "n", pch = c(4, 4, 4, 3, rep(3, 3)),
  col = c("black", (viridisLite::inferno(n = length(id0), begin = .5))[2:3],
    "black", (viridisLite::inferno(n = length(id0), begin = .5))[5:7]),
  fill = rep("transparent", 7), cex = 1.75,
  border = NA, x.intersp = -0.875, y.intersp = 0.4, legend = rep(NA, length(id0)))
text(5.8, -1.9, labels = "⟶ set order ⟶", cex = 1.2)
dev.off()
# Again, there's some Inkscape magic afterwards
