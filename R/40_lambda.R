
#' Options for the influence type
#'
#' @param type Character scalar with short-hand influence types.
#' @param position Integer scalar with the position of the variable of interest.
#' @param sign Integer scalar with the sign to apply to the influence.
#' @param f Optional custom function to calculate influence manually.
#'
#' @return Returns a function to compute lambda, with information in attributes.
#'
#' @export
set_lambda <- function(
  type = c("beta_i", "sigma_i", "se_i", "tstat_i",
    "cooksd", "dffits", "rstudent", "covratio", "BKW"),
  position = 1L,
  sign = 1L,
  f = NULL) {

  # Check custom functions
  if(!is.null(f)) {
    attr(f, "type") <- "custom"
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
  } else {
    f <- function(x, ...) {x[[type]][, position] * sign}
    attr(f, "type") <- type
    attr(f, "position") <- position
    attr(f, "sign") <- sign
  }

  return(f)
}


#' Options to set a target influence
#'
#' @param target Numerical vector with the target value(s).
#' @param type Character scalar with the tpye of inequality.
#' @param f Optional custom function to compare values.
#'
#' @return Returns a function to check whether a target value is reached.
#'
#' @export
set_target <- function(target = 0,
  type = c("less", "leq", "geq", "greater"),
  f = function(x, y, ...) {NULL}) {

  target <- num_check(target, -Inf, Inf, msg = "Please check the target.")

  # Check custom functions
  if(is.null(f())) {
    type <- match.arg(type)
    f <- list("less" = `<`, "leq" = `<=`, "geq" = `>=`, "greater" = `>` )[[type]]
  }

  attr(f, "target") <- target

  return(f)
}
