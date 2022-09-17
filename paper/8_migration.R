
# Migration ---------

library("dplyr")
library("ivreg")
library("sandwich")
library("lmtest")
devtools::load_all()


# Prepare -----

# Shorthand to target t-values
get_sens <- \(m, pos = 2L, target = c("tstat", "beta", "se")) {
  target <- match.arg(target)
  sign <- if(target == "se") -1 else sign(coef(m)[pos])
  if(class(m) == "ivreg") {
    sens.ivreg(m,
      lambda = set_lambda(target, pos = pos, sign = sign),
      options = set_options(n_max = 80L, tstat = TRUE, cluster = TRUE),
      cluster = seq_along(resid(m)))
  } else {
    sens.lm(m,
      lambda = set_lambda(target, pos = pos, sign = sign),
      options = set_options(n_max = 80L, tstat = TRUE, cluster = TRUE),
      cluster = seq_along(resid(m)))
  }
}

# Reproduce -----

# Reproduce baseline with HC0 standard errors ---
mdl_base <- readRDS("paper/data/migration_reg1.rds")
coeftest(mdl_base, vcov = vcovHC(mdl_base, "HC0"))

# Check sensitivity
sens_base <- get_sens(mdl_base)
plot.sensitivity(sens_base)
summary.sensitivity(sens_base) # 5/11/27 (Alg0 needs 7/30/-)

# Others are more robust to influence (Cribari-Neto)
coeftest(mdl_base, vcov = vcovHC(mdl_base, "HC4"))

# Compare adaptive and initial sets
overlap <- vector("numeric", 20L)
for(i in seq_along(overlap)) {
  overlap[i] <- sum(sens_base$influence$id[seq(i)] %in%
    sens_base$initial$id[seq(i)] / i)
}
overlap # 60% for 5, 60% for 10

# Reproduce plain model with HC0 standard errors ---
mdl_plain <- readRDS("paper/data/migration_reg2.rds")
coeftest(mdl_plain, vcov = vcovHC(mdl_plain, "HC0"))

# Check sensitivity
sens_plain <- get_sens(mdl_plain)
plot.sensitivity(sens_plain)
summary.sensitivity(sens_plain) # 6/11/32 (Alg0 needs 11/37/-)

# Others are more robust to influence (Cribari-Neto)
coeftest(mdl_plain, vcov = vcovHC(mdl_plain, "HC4"))

# Compare adaptive and initial sets
overlap <- vector("numeric", 20L)
for(i in seq_along(overlap)) {
  overlap[i] <- sum(sens_plain$influence$id[seq(i)] %in%
    sens_plain$initial$id[seq(i)] / i)
}
overlap # 40% for 5, 50% for 10

# First stages ---

# Base
data <- get_data.ivreg(mdl_base)
mdl_base_f <- lm(X[, 2] ~ Z[, -1], data = data)
coeftest(mdl_base_f, vcov = vcovHC(mdl_base_f, "HC0"))

# Plain
mdl_plain_f <- lm(X[, 2] ~ Z[, -1], data = get_data.ivreg(mdl_plain))

# Sensitivities
sens_base_f <- get_sens(mdl_base_f)
sens_plain_f <- get_sens(mdl_plain_f)

op <- par(mfrow = c(2, 1))
plot.sensitivity(sens_base_f)
plot.sensitivity(sens_plain_f)
par(op)
summary.sensitivity(sens_base_f) # 10/18/26 (Alg0 needs 15/34/-)
summary.sensitivity(sens_plain_f) # 10/18/31 (Alg0 needs 14/47/-)

# Numerics ---
sens_se1 <- get_sens(mdl_base, target = "se") # We can break it with 18
sens_se2 <- get_sens(mdl_plain, target = "se") # We can break it with 17
# Essentially, we're getting rid of the IV
sens_num <- sens.lm(mdl_base_f,
  lambda = set_lambda("custom", f = function(x, ...) {abs(x[["tstat_i"]][, 2])}),
  options = set_options(n_max = 80L, tstat = TRUE, cluster = TRUE),
  cluster = seq_along(resid(mdl_base_f)))
v <- sens_num$model[, c("beta_2", "se_2")]
v[, 1] / v[, 2]


# Investigate influential observations ---

df <- cbind("id" = seq_along(data$y), "gdp" = data$y,
  "pct_eu" = data$X[, 2], data$Z[, -1]) |> as_tibble()

# Remove observations
rm <- sens_base$influence$id[seq(4)]
m_rm <- ivreg(gdp ~ pct_eu + dist_ba + land_qu + rail_dens +
  pct_agr + pop_dens + rate_urb +
  rain_m + temp_m + elev_m + rel_m + ba + sf + er |
  iv1 + dist_ba + land_qu + rail_dens +
  pct_agr + pop_dens + rate_urb +
  rain_m + temp_m + elev_m + rel_m + ba + sf + er, data = df[-rm, ])
coeftest(m_rm, vcov = vcovHC(m_rm, "HC4")) # Needs 5 with HC0, 4 with HC4

# Full models
m_iv <- ivreg(gdp ~ pct_eu + dist_ba + land_qu + rail_dens +
  pct_agr + pop_dens + rate_urb + rain_m + temp_m + elev_m + rel_m +
  ba + sf + er | iv1 + dist_ba + land_qu + rail_dens + pct_agr + pop_dens +
  rate_urb + rain_m + temp_m + elev_m + rel_m + ba + sf + er, data = df)
m_lm <- lm(pct_eu ~ iv1 + dist_ba + land_qu + rail_dens + pct_agr + pop_dens +
  rate_urb + rain_m + temp_m + elev_m + rel_m + ba + sf + er, data = df)

# Condition index
svd(scale(df[, -1:-3]))$d^2 # max(d) / d following BKW
svd(scale(cbind(predict(m_lm), df[, c(-1:-4)])))$d^2
# Condition númber
rcond(as.matrix(df[, -1:-3]))
rcond(as.matrix(cbind(predict(m_lm), df[, c(-1:-4)])))
# Variance inflation
sort(car::vif(m_iv), decreasing = TRUE)
sort(car::vif(m_lm), decreasing = TRUE)

rm <- get_sens(mdl_base, target = "se")$influence$id
m_lm <- lm(pct_eu ~ iv1 + dist_ba + land_qu + rail_dens + pct_agr + pop_dens +
  rate_urb + rain_m + temp_m + elev_m + rel_m + ba + sf + er, data = df[-rm, ])
# Condition index
svd(scale(df[-rm, -1:-3]))$d^2
svd(scale(cbind(predict(m_lm), df[-rm, c(-1:-4)])))$d^2
# Condition númber
rcond(as.matrix(df[-rm, -1:-3]))
rcond(as.matrix(cbind(predict(m_lm), df[-rm, c(-1:-4)])))


# Outputs -----

tbl <- data.frame(
  v = c(
    "European share$^\\dag$", "", "Instrument$^\\dag$", "",
    "Geographic controls", "Socioeconomic controls", "Province fixed effects",
    "Thresholds$^\\dag$", "Observations", "$R^2$"
  ),
  m1_first = c(
    "", "",
    coeftest(mdl_base_f, vcov = vcovHC(mdl_base_f, "HC1"))[2, c(1, 3)] |>
    {\(x) c(round(x[1], 3), paste0("(", round(x[2], 2), ")"))}(),
    "Yes", "Yes", "Yes",
    "10[18]\\{26\\}", len(mdl_base$y), round(summary(mdl_base_f)$r.squared, 3)
  ),
  m1 = c(
    coeftest(mdl_base, vcov = vcovHC(mdl_base, "HC0"))[2, c(1, 3)] |>
    {\(x) c(round(x[1], 3), paste0("(", round(x[2], 2), ")"))}(),
    "", "",
    "Yes", "Yes", "Yes",
    "5[11]\\{27\\}", len(mdl_base$y), round(summary(mdl_base)$r.squared, 3)
  ),
  m2_first = c(
    "", "",
    coeftest(mdl_plain_f, vcov = vcovHC(mdl_plain_f, "HC1"))[2, c(1, 3)] |>
    {\(x) c(round(x[1], 3), paste0("(", round(x[2], 2), ")"))}(),
    "Yes", "No", "Yes",
    "10[18]\\{31\\}", len(mdl_plain$y), round(summary(mdl_plain_f)$r.squared, 3)
  ),
  m2 = c(
    coeftest(mdl_plain, vcov = vcovHC(mdl_plain, "HC0"))[2, c(1, 3)] |>
    {\(x) c(round(x[1], 3), paste0("(", round(x[2], 2), ")"))}(),
    "", "",
    "Yes", "No", "Yes",
    "6[11]\\{32\\}", len(mdl_plain$y), round(summary(mdl_plain)$r.squared, 3)
  )
)
tbl <- rbind(c("", rep(c("1st stage", "2nd stage"), 2)), tbl)
kableExtra::kbl(tbl, format = "latex", booktabs = TRUE, align = "lrrrr",
  escape = FALSE, label = "migration", position = "ht", row.names = FALSE,
  col.names = c("$\\log \\text{GDP}/\\text{capita} \\sim$",
    "Baseline", "", "Plain", ""),
  linesep = c(rep("", 4), "\\addlinespace", rep("", 2), "\\addlinespace",
    rep("", 3)),
  caption = "Sensitivity of long-term migration impacts.")


# Plot marginalised instrument ----
d_fwl <- update_fwl(data$Z, data$X[, 2], 2)
r_fwl <- update_fwl(data$Z[-rm, ], data$X[-rm, 2], 2)

cairo_pdf("paper/output/migration-iv.pdf", height = 2.4, width = 7)
op <- par(mfrow = c(1, 2), mar = c(3, 4, 0, 1), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(d_fwl$X) + c(-.01, .01), ylim = range(d_fwl$y) + c(-.01, .01))
points(x = d_fwl$X, y = d_fwl$y, pch = 1, cex = 1.25)
points(x = d_fwl$X[rm], y = d_fwl$y[rm], pch = 4, cex = 1.75, lwd = 2, col = "#008080")
axis(1, at = round(c(min(d_fwl$X), 0, max(d_fwl$X)), 2))
axis(2, at = round(c(min(d_fwl$y), 0, max(d_fwl$y)), 2), las = 1)
abline(v = c(min(d_fwl$X), 0, max(d_fwl$X)), col = "darkgray", lty = 2)
abline(h = c(min(d_fwl$y), 0, max(d_fwl$y)), col = "darkgray", lty = 2)
abline(lm(y ~ X - 1, data = d_fwl), lwd = 2, col = "#800000")
title(ylab = "Residualized outcome", line = 3, font.lab = 2)
title(xlab = "Marginalized instrument", line = 2, font.lab = 2)
plot.new()
plot.window(xlim = range(d_fwl$X) + c(-.01, .01), ylim = range(d_fwl$y) + c(-.01, .01))
points(x = r_fwl$X, y = r_fwl$y, pch = 1, cex = 1.25)
axis(1, at = round(c(min(d_fwl$X), 0, max(d_fwl$X)), 2))
axis(2, at = round(c(min(d_fwl$y), 0, max(d_fwl$y)), 2), las = 1)
abline(v = c(min(d_fwl$X), 0, max(d_fwl$X)), col = "darkgray", lty = 2)
abline(h = c(min(d_fwl$y), 0, max(d_fwl$y)), col = "darkgray", lty = 2)
abline(lm(y ~ X - 1, data = r_fwl), lwd = 2, col = "#008080")
title(xlab = "Marginalized instrument", line = 2, font.lab = 2)
dev.off()
