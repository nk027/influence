
rank_influence <- function(x, lambda) {
  value <- lambda(x)
  order <- order(value, decreasing = FALSE, method = "radix")
  cbind("value" = value, "order" = order)
}


set_lambda <- function(
  type = c("beta_i", "sigma_i", "se_i", "tstat_i",
    "cooksd", "dffits", "rstudent", "covratio", "BKW"),
  sign = 1L,
  position = 1L,
  f = function(x, ...) {NULL}) {

  # Check custom functions
  if(!is.null(f())) {return(f)}

  # Allow choosing some predefined functions
  type <- match.arg(type)
  sign <- sign(sign)
  position <- int_check(position, 1L, 1e6L, "Choose a valid position.")

  scalars <- c("sigma_i", "cooksd", "dffits", "rstudent", "covratio", "BKW")
  if(type %in% scalars && position != 1L) {
    warning("Scalar type chosen, setting position to one.")
    position <- 1L
  }

  if(type == "BKW") {
    f <- function(x, ...) {
      vapply(seq_len(NROW(x$influence$beta_i)), function(i) {
        crossprod(-x[["beta_i"]][i, ] + x[["lm"]][["beta"]])
      }, numeric(1L)) * sign
    }
  } else {
    f <- function(x, ...) {x[[type]][, position] * sign}
  }

  return(f)
}
