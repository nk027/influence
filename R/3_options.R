
set_options <- function(
  p_max = NULL, n_max = NULL,
  sm_re = 5L, fwl = 0L, fwl_re = 1e6L,
  ...) {

  out <- set_compute(...)

  if(is.null(p_max)) {
    if(is.null(n_max)) {
      p_max <- .5
      n_max <- 100L
    } else {
      p_max = 1
    }
  } else if(is.null(n_max)) {
    n_max <- Inf
  }

  out$p_max = num_check(p_max, 0, 1, msg = "Choose a valid maximum percentage.")
  out$n_max = int_check(n_max, 1L, Inf, msg = "Choose a valid maximum number.")
  out$sm_re = int_check(sm_re, 0L, 1e6L,
    msg = "Choose a valid step size for Sherman-Morrison updates.")
  out$fwl = vapply(fwl, int_check, integer(1L), 0L, Inf,
    msg = "Choose valid indices for variables to retain after FWL.")
  out$fwl_re = int_check(fwl_re, 1L, 1e6L,
    msg = "Choose a valid step size for recalculating FWL.")

  return(out)
}


set_compute <- function(x = c("all", "some", "none"),
  hat, beta, sigma, se, tstat, cooksd, dffits, rstudent, covratio) {

  x <- match.arg(x)

  out <- if(x == "all") {
    list("hat" = TRUE, "beta" = TRUE, "sigma" = TRUE,
      "se" = TRUE, "tstat" = TRUE,
      "cooksd" = TRUE, "dffits" = TRUE, "rstudent" = TRUE, "covratio" = TRUE)
  } else if(x == "some") {
    list("hat" = TRUE, "beta" = TRUE, "sigma" = TRUE,
      "se" = FALSE, "tstat" = FALSE, "cooksd" = FALSE, "dffits" = FALSE,
      "rstudent" = FALSE, "covratio" = FALSE)
  } else if(x == "none") {
    list("hat" = FALSE, "beta" = FALSE, "sigma" = FALSE,
      "se" = FALSE, "tstat" = FALSE, "cooksd" = FALSE, "dffits" = FALSE,
      "rstudent" = FALSE, "covratio" = FALSE)
  }

  if(!missing(hat)) {out$hat <- isTRUE(hat)}
  if(!missing(beta)) {out$beta <- isTRUE(beta)}
  if(!missing(sigma)) {out$hat <- isTRUE(sigma)}
  if(!missing(se)) {out$hat <- isTRUE(se)}
  if(!missing(tstat)) {out$hat <- isTRUE(tstat)}
  if(!missing(cooksd)) {out$hat <- isTRUE(cooksd)}
  if(!missing(dffits)) {out$hat <- isTRUE(dffits)}
  if(!missing(rstudent)) {out$hat <- isTRUE(rstudent)}
  if(!missing(covratio)) {out$hat <- isTRUE(covratio)}

  if(out$se) {out$beta <- out$sigma <- TRUE}
  if(out$tstat) {out$se <- out$beta <- out$sigma <- TRUE}
  if(out$cooksd || out$dffits || out$rstudent || out$covratio) {out$hat <- TRUE}
  if((out$dffits || out$rstudent || out$covratio )&& !out$sigma) {
    warning("Consider computing observation-specific sigma.")
  }

  return(out)
}
