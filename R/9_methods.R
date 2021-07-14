
init <- function(x, ...) {{UseMethod("init", x)}}

init.sensitivity <- function(x, id = NULL, ...) {

  if(is.null(id)) { # Try to guess
    type <- gsub("^([a-z]+).*", "\\1", attr(x$meta$lambda, "type"))
    position <- attr(x$meta$lambda, "position")
    id <- paste0(type, "_", position)
  }

  exact <- if(grepl("tstat", id)) { # Exact
    x$model[[paste0("beta_", gsub(".*_([0-9]+)", "\\1", id))]] /
      x$model[[paste0("se_", gsub(".*_([0-9]+)", "\\1", id))]]
  } else {
    x$model[[id]]
  }
  initial <- if(attr(x$meta$lambda, "sign") == -1L) { # Initial
    cumsum(c(exact[1L], -(exact[1L] + x$initial$lambda)))
  } else {
    cumsum(c(exact[1L], -(exact[1L] - x$initial$lambda)))
  }

  list("initial" = initial, "exact" = exact, "id" = id)
}


print.influence <- function(x, ...) {
  cat("Influence object\n")
  print(str(x))
  invisible(x)
}

print.sensitivity <- function(x, ...) {
  cat("Sensitivity object\n")
  print(str(x[c("influence", "model")]))
  invisible(x)
}


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
  type = c("path", "influence"), ...) {

  type <- match.arg(type)
  if(type == "influence") {
    return(.plot_influence(x, ...))
  } else if(type == "path") {
    return(.plot_path(x, ...))
  }
}

.plot_influence <- function(x) {

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

.plot_path <- function(x, n = 0L) {

  z <- init(x)
  z$exact <- z$exact[!is.nan(z$exact)]
  z$initial <- z$initial[!is.nan(z$initial)]

  N <- x$influence$N[1L]

  poi <- list("exact" = c(which(diff(sign(z$exact)) != 0) + 1,
      which(diff(abs(z$exact) > 2) != 0) + 1),
    "initial" = c(which(diff(sign(z$initial)) != 0) + 1,
      which(diff(abs(z$initial) > 2) != 0) + 1))

  if(n > 0L) {
    ylim <- c(min(z$exact[seq.int(n)], z$initial[seq.int(n)]),
      max(z$exact[seq.int(n)], z$initial[seq.int(n)]))
    plot(z$initial[seq.int(n)], type = "l", col = "gray", lty = 2,
      ylim = ylim, ylab = "Value", xlab = "Index / Percent")
    lines(z$exact)
    axis_at <- axTicks(3L)
    axis_lab <- round((seq(0, n) / N)[axTicks(3L)], 2)
    if(any(axis_at == 0)) {axis_lab <- c(0, axis_lab)}
    axis(3L, at = axis_at, labels = axis_lab)
  } else {
    ylim <- c(min(z$exact, z$initial), max(z$exact, z$initial))
    plot(z$exact, type = "l", ylab = "Value",
      ylim = ylim, xlab = "Index / Percent")
    lines(z$initial, col = "gray", lty = 2)
    axis_at <- axTicks(3L)
    axis_lab <- round((seq(0, N) / N)[axTicks(3L)], 2)
    if(any(axis_at == 0)) {axis_lab <- c(0, axis_lab)}
    axis(3L, at = axis_at, labels = axis_lab)
  }
  grid()
  abline(v = poi$initial, col = "gray", lty = 2)
  abline(v = poi$exact)
  axis(1L, at = poi$exact, labels = TRUE)
  abline(h = 0)

  invisible(x)
}
