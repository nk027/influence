
divergence <- function(a, b) {

  x <- vector("numeric", length(a))
  for(i in seq_along(a)) {
    x[i] <- i - sum(a[seq(i)] %in% b[seq(i)])
  }
  x
}
