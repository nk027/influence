
# Tsetse fly ---------

library("sandwich")
library("lmtest")
library("MASS")
library("robustbase")
library("robust")
devtools::load_all()


# Prepare ---

# Load data
x <- foreign::read.dta("paper/data/tsetse_data.dta")
cluster <- x$province

# Short-hand to target t-values
get_sens <- \(m, pos) {
  sens.lm(m,
    lambda = set_lambda("tstat", pos = pos, sign = sign(coef(m)[pos])),
    options = set_options(n_max = 100, tstat = TRUE, hat = TRUE, cluster = TRUE),
    cluster = cluster[-m$na.action])
}

# Short-hand function to link observations and other info
get_info <- \(m, id, info = "isocode") { # province, name, c("lon", "lat")
  x[rownames(m$model)[id], info]
}

# Short-hand to get robust t-values
get_t <- \(m, pos = 2) {
  coeftest(m, vcov = vcovCL(m, "HC1", cluster = cluster))[pos, 3]
}

# Short-hand to summarise thresholds
get_thresholds <- \(s) {
  pos <- summary.sensitivity(s)$exact
  if(length(pos) == 3) {
    pos <- pos[c("threshold1", "zero", "threshold2")]
  } else {pos <- c("--", pos[c("zero", "threshold")])}
  paste0(pos[1], "[", pos[2], "]", "{", pos[3], "}")
}


# Reproduce results of Alsan (2015) -----

# Prepare the formula to loop over
lhs <- paste(c("animals", "intensive", "plow", "female_ag",
  "ln_popd_murdock", "slavery_indigenous", "central"))
rhs <- "~ TSI + prop_tropics + meantemp + meanrh + itx"
controls <- "+ malaria_index + coast + river + lon + abslat + meanalt + SI"

# Reproduce results ---
mdl_base <- vector("list", length(lhs))
# mdl_plain <- vector("list", length(lhs))
for(i in seq_along(lhs)) {
  mdl_base[[i]] <- lm(as.formula(paste(lhs[i], rhs, controls)), data = x)
  # mdl_plain[[i]] <- lm(as.formula(paste(lhs[i], rhs)), data = x)
}
lapply(mdl_base, \(m) coeftest(m, vcov = vcovCL(m, "HC1", cluster = cluster)))

# Check robust estimation ---
r_base <- vector("list", length(lhs))
for(i in seq_along(lhs)) {
  r_base[[i]] <- rlm(as.formula(paste(lhs[i], rhs, controls)), data = x,
    method = "M", psi = psi.huber, maxit = 1000)
}
lapply(r_base, \(m) coeftest(m, vcov = vcovCL(m, "HC1", cluster = cluster)))
# Pathological behaviour for the plow and slavery as outcomes

# Check S-estimation ---
rs_base <- vector("list", length(lhs))
for(i in seq_along(lhs)) {
  rs_base[[i]] <- lmRob(as.formula(paste(lhs[i], rhs, controls)), data = x)
}
lapply(rs_base, summary) # Pathological for all but population
# A lot of observations are removed as potential outliers

# Check sensitivity of the main results ---
sens_base <- vector("list", length(mdl_base))
# sens_plain <- vector("list", length(mdl_plain))
for(i in seq_along(lhs)) {
  sens_base[[i]] <- get_sens(mdl_base[[i]], 2)
  # sens_plain[[i]] <- get_sens(mdl_plain[[i]], 2)
}
lapply(sens_base, summary.sensitivity)


# Outputs -----

# Table with main results ---
tbl <- data.frame(i = c("Baseline", rep("", 2), "Observations", "Clusters", "$R^2$"))

# Add columns based on the separate models
names <- c("Animals", "Intensive", "Plow", "Female", "Urbanisation", "Slavery", "Centralization")
for(i in seq_along(names)) {
  tbl[[names[i]]] <- c(
    round(coef(mdl_base[[i]])[2], 3), paste0("(", round(get_t(mdl_base[[i]]), 2), ")"),
    get_thresholds(sens_base[[i]]),
    # round(coef(mdl_plain[[i]])[2], 3), paste0("(", round(get_t(mdl_plain[[i]]), 2), ")"),
    # get_thresholds(sens_plain[[i]]),
    length(mdl_base[[i]]$residuals),
    cluster[as.integer(rownames(mdl_base[[i]]$model))] |> unique() |> length(),
    round(summary(mdl_base[[i]])$r.squared, 3)
  )
}

kableExtra::kbl(tbl, format = "latex", booktabs = TRUE, align = "lccccccc",
  escape = FALSE, label = "tsetse", position = "ht", row.names = FALSE,
  linesep = c(rep("", 1), "\\addlinespace", rep("", 3), "\\addlinespace", rep("", 3)),
  caption = "The effects of the Tsetse fly.")

# Table with robust results ---
tbl <- data.frame(i = c("M-estimation", "", "S-estimation", "", "Observations"))

# Add columns based on the separate models
names <- c("Animals", "Intensive", "Plow", "Female", "Urbanisation", "Slavery", "Centralization")
for(i in seq_along(names)) {
  tbl[[names[i]]] <- c(
    round(coef(r_base[[i]])[2], 3),
    paste0("(", round(get_t(r_base[[i]]), 2), ")"),
    round(coef(rs_base[[i]])[2], 3),
    paste0("(", round(summary(rs_base[[i]])$coefficients[2, 3], 2), ")"),
    length(mdl_base[[i]]$residuals)
  )
}

# M-estimates are screwed for the plow and slavery
tbl[1:2, c(4, 7)] <- matrix("--", 2, 2)
# S-estimates are screwed for all, but urbanisation
tbl[3:4, c(-1, -6)] <- matrix("--", 2, 6)

kableExtra::kbl(tbl, format = "latex", booktabs = TRUE, align = "lccccccc",
  escape = FALSE, label = "tsetse_rob", position = "ht", row.names = FALSE,
  linesep = c(rep("", 1), "\\addlinespace", rep("", 1), "\\addlinespace", rep("", 3)),
  caption = "Robust estimation for the effects of the Tsetse fly")
