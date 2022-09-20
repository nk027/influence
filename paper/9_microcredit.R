
# Microcredit ---------

devtools::load_all()


# Prepare ---

# Load data
x <- readRDS("paper/data/microcredit.rds")

# Short-hand to prepare sensitivity summaries
get_sm <- \(s, threshold = qnorm(.975)) {
  out <- unlist(summary.sensitivity(s, threshold = threshold))
  if(is.na(threshold)) {
    out[grepl("zero[1]*$", names(out))]
  } else {
    out[grepl("threshold[1]*$", names(out))]
  }
}
# Short-hand for the BGM approach to significance
approx_signif <- \(m, approx = TRUE, threshold = qnorm(.975)) {
  infl <- infl.lm(m, options = set_options(if(approx) "none" else "all"))
  sign <- if(coef(m)[2] >= 0) 1 else 2
  comparison <- c(`>`, `<`)[[sign]]
  base <- infl$model$beta[2] + sign(coef(m)[2]) * threshold * infl$model$se[2]
  signif_i <- (infl$model$beta[2] - infl$beta_i[, 2]) +
    sign(coef(m)[2]) * threshold * (infl$model$se[2] - infl$se_i[, 2])
  # The number to remove is the first one past the base
  which(comparison(cumsum(sort(signif_i, decreasing = (sign == 1))), base))[1]
}


# Reproduce -----

# Estimate the models ---
mdls <- lapply(x, \(x)
  lm(y ~ D, data = x[complete.cases(x), ], x = TRUE, y = TRUE))
lapply(mdls, summary)

# First, we focus on sign-switches (with and without approximation) ---
exact_sign <- lapply(mdls, \(m) sens.lm(m,
  lambda = set_lambda("beta", position = 2L, sign = sign(coef(m)[2])),
  options = set_options(hat = TRUE, beta = TRUE)))
approx_sign <- lapply(mdls, \(m) sens.lm(m,
  lambda = set_lambda("beta", position = 2L, sign = sign(coef(m)[2])),
  options = set_options("none")))
# Check results (we only want the zeros)
ns_ex <- sapply(exact_sign, \(s) get_sm(s, threshold = NA))
ns_ap <- sapply(approx_sign, \(s) get_sm(s, threshold = NA))
# Alg2 improves in two cases, there's no issue with the approximation here

# Second, we focus on significance (with and without approximations) ---
exact_t <- lapply(mdls, \(m) sens.lm(m,
  lambda = set_lambda("tstat", position = 2L, sign = sign(coef(m)[2])),
  options = set_options(hat = TRUE, beta = TRUE, se = TRUE, tstat = TRUE)))
approx_t <- lapply(mdls, \(m) sens.lm(m,
  lambda = set_lambda("tstat", position = 2L, sign = sign(coef(m)[2])),
  options = set_options("none")))
# Check results (a bit messier to get the right thresholds)
nt_ex <- sapply(exact_t, \(s) get_sm(s, threshold = qnorm(.975)))
nt_ap <- sapply(approx_t, \(s) get_sm(s, threshold = qnorm(.975)))
# Alg2 improves everywhere, the approximation is not too bad

# To reproduce BGM, we also need to use their significance check ---
nt2_ex <- sapply(mdls, approx_signif, approx = FALSE, threshold = qnorm(.975))
nt2_ap <- sapply(mdls, approx_signif, approx = TRUE, threshold = qnorm(.975))
# Interestingly, this appr7oach overshoots, which balances out the underestimates

# Outputs -----

# Table with main results ---
tbl <- data.frame(i = c("Estimate$^\\dag$", "", "Sign-switch", "", "",
  "Significant sign-switch", "", "", "Observations"),
  j = c(rep("", 9)))
# Add information
for(i in seq_along(x)) {
  tbl[[names(x)[i]]] <- c(
    round(coef(mdls[[i]])[2], 3),
    paste0("(", round(summary(mdls[[i]])$coefficients[2, 3], 2), ")"),
    ns_ex[1, i], ns_ex[2, i], ns_ap[2, i],
    nt_ex[1, i], nt_ex[2, i], nt2_ap[i],
    length(mdls[[i]]$y)
  )
}

kableExtra::kbl(tbl, format = "latex", booktabs = TRUE, align = "llccccccc",
  escape = FALSE, label = "microcredit", position = "th", row.names = FALSE,
  linesep = c(rep("", 1), "\\addlinespace", rep("", 2), "\\addlinespace", rep("", 1)),
  caption = "Sensitivity of the average treatment effect of microcredits")


# Plot the data and regression lines ---
mode <- function(v) {uv <- unique(v); uv[which.max(tabulate(match(v, uv)))]}
set.seed(42)
cairo_pdf("paper/output/microcredit_data.pdf", height = 2, width = 4, pointsize = 8)
op <- par(mfrow = c(1, 2), mar = c(2, 4, 2, .5), bg = "white", family = "Noto Sans")
for(j in c("ETH", "MEX")) {
  y_vals <- mdls[[j]]$y
  x_vals <- mdls[[j]]$x[, 2]
  plot.new()
  plot.window(xlim = c(-.2, 1.2), ylim = range(mdls[[j]]$y) + c(-10, 10))
  abline(v = c(0, 1), col = "gray", lty = 3)
  abline(h = range(y_vals), col = "gray", lty = 3)
  id <- (y_vals - mode(y_vals))^2 > 10000 # We can sample close to the median
  points(x_vals[id] + rnorm(sum(id), 0, .03), y_vals[id], pch = 1, lwd = 1)
  # We don't to plot all the values close to the mode
  points(x_vals[sample(which(!id), min(sum(!id), 100))] + rnorm(min(100, sum(!id)), 0, .03),
    y_vals[sample(which(!id), min(sum(!id), 100))], pch = 1, lwd = 1)
  # Indicate the mean and median -- they're basically invisible
  # x_0 <- x_vals == 0
  # segments(-.1, quantile(y_vals[x_0], .1), .1, col = "#00A0A0", lty = 2, lwd = 2)
  # segments(-.1, quantile(y_vals[x_0], .5), .1, col = "#00A0A0", lty = 1, lwd = 3)
  # segments(-.1, quantile(y_vals[x_0], .9), .1, col = "#00A0A0", lty = 2, lwd = 2)
  # segments(-.1, mean(y_vals[x_0]), .1, col = "#A000A0", lty = 1, lwd = 3)
  # x_1 <- x_vals == 1
  # segments(.9, quantile(y_vals[x_1], .1), 1.1, col = "#00A0A0", lty = 2, lwd = 2)
  # segments(.9, quantile(y_vals[x_1], .5), 1.1, col = "#00A0A0", lty = 1, lwd = 3)
  # segments(.9, quantile(y_vals[x_1], .9), 1.1, col = "#00A0A0", lty = 2, lwd = 2)
  # segments(.9, mean(y_vals[x_1]), 1.1, col = "#A000A0", lty = 1, lwd = 3)
  abline(mdls[[j]], col = "darkgray", lwd = 2)
  axis(1, at = c(0, 1))
  axis(2, at = round(c(min(mdls[[j]]$y), 0, max(mdls[[j]]$y)), 0),
    labels = round(c(min(mdls[[j]]$y), 0, max(mdls[[j]]$y)) / 1000, 1), las = 1)
  title(main = j, line = 0, font.main = 2, cex.main = 1.5, family = "Merriweather")
}
dev.off()
