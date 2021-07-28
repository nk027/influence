
set_lambda <- function(
  type = c("beta_i", "sigma_i", "se_i", "tstat_i",
    "cooksd", "dffits", "rstudent", "covratio", "BKW"),
  position = 1L,
  sign = 1L,
  target = 0L,
  f = function(x, ...) {NULL}) {

  target <- num_check(target, -Inf, Inf, msg = "Please check the target.")

  # Check custom functions
  if(!is.null(f())) {
    attr(f, "type") <- "custom"
    attr(f, "target") <- target
    return(f)
  }

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
      vapply(seq_len(NROW(x$beta_i)), function(i) {
        crossprod(-x[["beta_i"]][i, ] + x[["model"]][["beta"]])
      }, numeric(1L)) * sign
    }
    attr(f, "type") <- "BKW"
    attr(f, "sign") <- sign
    attr(f, "target") <- target
  } else {
    f <- function(x, ...) {x[[type]][, position] * sign}
    attr(f, "type") <- type
    attr(f, "position") <- position
    attr(f, "sign") <- sign
  }
  attr(f, "target") <- target

  return(f)
}


set_delta <- function(type = c("less", "leq", "geq", "greater"),
  f = function(x, y, ...) {NULL}) {

  # Check custom functions
  if(!is.null(f())) {
    return(f)
  }

  # Basic stuff otherwise
  type <- match.arg(type)

  f <- list("less" = `<`, "leq" = `<=`, "geq" = `>=`, "greater" = `>` )[[type]]

  return(f)
}
