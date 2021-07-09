
set_options <- function(
  p_max = NULL, n_max = NULL,
  sm_re = 5L, fwl = 0L, fwl_re = 1e6L,
  hat = TRUE, beta = TRUE, sigma = TRUE, se = TRUE, tstat = TRUE,
  cooksd = TRUE, dffits = TRUE, rstudent = TRUE, covratio = TRUE) {

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

  out <- list(
    sensitivity = list("hat" = isTRUE(hat),
      "beta" = isTRUE(beta), "sigma" = isTRUE(sigma),
      "se" = isTRUE(se), "tstat" = isTRUE(tstat),
      "cooksd" = isTRUE(cooksd), "dffits" = isTRUE(dffits),
      "rstudent" = isTRUE(rstudent), "covratio" = isTRUE(covratio)
    ),
    p_max = num_check(p_max, 0, 1, msg = "Choose a valid maximum percentage."),
    n_max = int_check(n_max, 1L, Inf, msg = "Choose a valid maximum number."),
    sm_re = int_check(sm_re, 0L, 1e6L,
      msg = "Choose a valid step size for Sherman-Morrison updates."),
    fwl = vapply(fwl, int_check, integer(1L), 0L, Inf,
      msg = "Choose valid indices for variables to retain after FWL."),
    fwl_re = int_check(fwl_re, 1L, 1e6L,
      msg = "Choose a valid step size for recalculating FWL.")
  )
}
