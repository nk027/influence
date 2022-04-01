
plot.influence <- function(x,
  type = c("beta_i", "se_i"), position,
  ...) {

  type <- match.arg(type)

  if(missing(position)) {position <- seq_len(NCOL(x[[type]]))}
  n_plots <- length(position)
  boxplot(x[[type]][, position], main = "Variation after one removal")
  xl <- seq(n_plots) - 0.45
  xr <- seq(n_plots) + 0.45
  segments(x0 = xl, x1 = xr,
    y0 = x$lm[[gsub("([a-z]+)_i", "\\1", type)]][position],
    col = "#800000", lty = 3, lwd = par("lwd") * 2)

  invisible(x)
}


plot.sensitivity <- function(x,
  type = c("path", "masking"), ...) {

  type <- match.arg(type)
  if(type == "masking") {
    return(.plot_masking(x, ...))
  } else if(type == "path") {
    return(.plot_path(x, ...))
  }
}

.plot_masking <- function(x) {

  masking <- vector("numeric", nrow(x$influence))
  for(i in seq_len(nrow(x$influence))) {
    masking[i] <- i -
      sum(x$influence$id[seq(i)] %in% x$initial$id[seq(i)])
  }
  masked <- 1L + which(diff(masking) > 0)
  plot(x$influence$lambda, type = "l",
    xlab = "Index / Number masked", ylab = "Influence")
  axis(3L, at = masked, labels = masking[masked])
  abline(h = 0)
  grid()

  invisible(x)
}

.plot_path <- function(x, n = 0L, threshold = 2) {

  z <- init(x)
  z$exact <- z$exact[!is.nan(z$exact)]
  z$initial <- z$initial[!is.nan(z$initial)]

  N <- x$influence$N[1L]

  poi <- list("exact" = c(which(diff(sign(z$exact)) != 0) + 0,
      which(diff(abs(z$exact) > threshold) != 0) + 0),
    "initial" = c(which(diff(sign(z$initial)) != 0) + 0,
      which(diff(abs(z$initial) > threshold) != 0) + 0))

  if(n > 0L) {
    ylim <- c(min(z$exact[seq.int(n)], z$initial[seq.int(n)]),
      max(z$exact[seq.int(n)], z$initial[seq.int(n)]))
    plot(z$initial[seq.int(n)], x = seq_along(z$exact) - 1,
      type = "l", col = "gray", lty = 2, ylim = ylim,
      ylab = "Value", xlab = "Index / Percent")
    lines(z$exact, x = seq_along(z$exact) - 1)
    axis_at <- axTicks(3L)
    axis_lab <- round((seq(0, n) / N)[axTicks(3L)], 2)
    if(any(axis_at == 0)) {axis_lab <- c(0, axis_lab)}
    axis(3L, at = axis_at, labels = axis_lab)
  } else {
    ylim <- c(min(z$exact, z$initial), max(z$exact, z$initial))
    plot(z$exact, x = seq_along(z$exact) - 1, type = "l", ylab = "Value",
      ylim = ylim, xlab = "Index / Percent")
    lines(z$initial, x = seq_along(z$initial) - 1, col = "darkgray", lty = 2)
    axis_at <- axTicks(3L)
    axis_lab <- round((seq(0, N) / N)[axTicks(3L)], 2)
    if(any(axis_at == 0)) {axis_lab <- c(0, axis_lab)}
    axis(3L, at = axis_at, labels = axis_lab)
  }
  grid()
  abline(v = poi$exact, lty = 1)
  abline(v = poi$initial + .05, col = "darkgray", lty = 2)
  axis(1L, at = poi$initial, labels = TRUE, tick = FALSE, padj = 1)
  axis(1L, at = poi$exact, labels = TRUE, font = 2, padj = -1)
  abline(h = 0)

  invisible(x)
}
