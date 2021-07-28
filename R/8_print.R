
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
