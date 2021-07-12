
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

plot.sensitivity <- function(x) {

  masking <- vector("numeric", nrow(x$influence))
  for(i in seq_len(nrow(x$influence))) {
    masking[i] <- i -
      sum(x$influence$id[seq(i)] %in% x$influence$init_id[seq(i)])
  }
  masked <- 1L + which(diff(masking) > 0)
  plot(x$influence$lambda, type = "l",
    xlab = "Index / Number masked", ylab = "Influence")
  axis(3L, at = masked, labels = masking[masked])

  invisible(x)
}
