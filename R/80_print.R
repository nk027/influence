
print.influence <- function(x, ...) {
  cat("Influence object\n")
  print(str(x))
  invisible(x)
}

print.init <- function(x, ...) {
  cat("Sensitivity object\n")
  print(str(x[c("influence", "model")]))
  invisible(x)
}

print.goal <- function(x, ...) {
  cat("Sensitivity object\n")
  print(str(x[c("influence", "model")]))
  invisible(x)
}

print.sens <- function(x, ...) {
  cat("Sensitivity object\n")
  print(str(x[c("influence", "model")]))
  invisible(x)
}
