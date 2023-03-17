
#' @noRd
#' @export
print.influence <- function(x, ...) {
  cat("Influence object\n")
  print(str(x))
  invisible(x)
}

#' @noRd
#' @export
print.init <- function(x, ...) {
  cat("Sensitivity object\n")
  print(str(x[c("influence", "model")]))
  invisible(x)
}

#' @noRd
#' @export
print.goal <- function(x, ...) {
  cat("Sensitivity object\n")
  print(str(x[c("influence", "model")]))
  invisible(x)
}

#' @noRd
#' @export
print.sens <- function(x, ...) {
  cat("Sensitivity object\n")
  print(str(x[c("influence", "model")]))
  invisible(x)
}
