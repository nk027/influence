
# Poverty convergence ---

library("dplyr")
devtools::load_all()

mdl_rav <- readRDS("paper/data/ravallion2012.rds")

# Check sensitivity
sens_rav <- sens.lm(mdl_rav, lambda = set_lambda("tstat", pos = 2, sign = 1),
  options = set_options(n_max = 20L, tstat = TRUE, cluster = TRUE))
# plot.sensitivity(sens)
summary.sensitivity(sens_rav)

# Redo with removals
ids <- sens_rav$influence$id[1:4]
data <- get_data.lm(mdl_rav)
m_rm <- lm(y[-ids] ~ X[-ids, ] - 1, data = data)

x <- data$X[, 2]
y <- data$y
countries <- c(
  "Albania", "Algeria", "Argentina", "Armenia", "Bangladesh", "Belarus",
  "Bolivia", "Botswana", "Brazil", "Burkina Faso", "Burundi", "Cambodia",
  "Cameroon", "Central African Republic", "Chile", "China", "Colombia",
  "Costa Rica", "Côte d'Ivoire", "Djibouti", "Dominican Republic", "Ecuador",
  "Egypt, Arab Rep.", "El Salvador", "Estonia", "Ethiopia", "Gambia, The",
  "Georgia", "Ghana", "Guatemala", "Guinea", "Guinea-Bissau", "Guyana",
  "Honduras", "India", "Indonesia", "Iran, Islamic Rep.", "Jamaica", "Jordan",
  "Kazakhstan", "Kenya", "Lao PDR", "Latvia", "Lesotho", "Lithuania",
  "Macedonia, FYR", "Madagascar", "Malawi", "Malaysia", "Mali", "Mauritania",
  "Mexico", "Moldova", "Mongolia", "Morocco", "Mozambique", "Nepal",
  "Nicaragua", "Niger", "Nigeria", "Pakistan", "Panama", "Paraguay", "Peru",
  "Philippines", "Poland", "Romania", "Russian Federation", "Rwanda", "Senegal",
  "Sierra Leone", "South Africa", "Sri Lanka", "Swaziland", "Tajikistan",
  "Tanzania", "Thailand", "Trinidad and Tobago", "Tunisia", "Turkey",
  "Turkmenistan", "Uganda", "Ukraine", "Uruguay  ", "Uzbekistan",
  "Venezuela, RB", "Vietnam", "Yemen, Rep.", "Zambia")
countries[ids]

# Check sensitivity with EE interaction
ee <- countries %in% c("Belarus", "Estonia", "Georgia", "Kazakhstan", "Latvia",
  "Lithuania", "Moldova", "Poland", "Romania", "Russian Federation", "Ukraine")
mdl_ee <- lm(y ~ x * ee)
sens_ee <- sens.lm(mdl_ee,
  lambda = set_lambda("tstat", pos = 2, sign = -1),
  options = set_options(n_max = 50L, tstat = TRUE, cluster = TRUE))
plot.sensitivity(sens_ee)
summary.sensitivity(sens_ee)

# Plot the sensitive set
cairo_pdf("paper/output/pov-convergence.pdf", height = 4.2, width = 4)
op <- par(mar = c(2, 3, 0, .5), bg = "transparent", family = "Noto Sans")
plot.new()
plot.window(xlim = range(x) + c(-.2, .01), ylim = range(y) + c(-.01, .01))
points(x, y, cex = 1.25)
axis(1, at = round(c(min(x), max(x)), 2))
axis(2, at = round(c(min(y), 0, max(y)), 2), las = 1)
# Add regression line
abline(coef(m), col = "darkgray", lty = 1, lwd = 2)
abline(coef(m_rm), col = "#008080", lty = 2, lwd = 2)
# Highlight sets
points(x, y, cex = 1.25)
points(x[ids], y[ids], cex = 1.5, col = "darkgray",
  pch = 21, bg = (viridisLite::inferno(n = length(ids), begin = 0.5)))
points(x[ids[1:2]], y[ids[1:2]], cex = 1.75, pch = 4)
points(x[ids[-1:-2]], y[ids[-1:-2]], cex = 1.75, pch = 3)
text(x = x[ids[1]] - .1, y = y[ids[1]] + .032,
  labels = countries[ids[1]], cex = 1.5)
text(x = x[ids[2]] + .2, y = y[ids[2]] + .032,
  labels = countries[ids[2]], cex = 1.5)
text(x = x[ids[3]] - .8, y = y[ids[3]] - .02,
  labels = countries[ids[3]], cex = 1.5)
text(x = x[ids[4]] + .75, y = y[ids[4]] + .02,
  labels = countries[ids[4]], cex = 1.5)
dev.off()


# With PovCalNet data ---

mdl_cre <- readRDS("paper/data/poverty.rds")

# Check sensitivity
sens_cre <- sens.lm(mdl_cre, lambda = set_lambda("tstat", pos = 2, sign = -1),
  options = set_options(n_max = 100, tstat = TRUE, cluster = TRUE))
plot.sensitivity(sens_cre)
summary.sensitivity(sens_cre)

id0 <- sens_cre$initial$id[1:26]
id2 <- sens_cre$influence$id[1:26]
data <- get_data.lm(mdl_cre)
m_rm0 <- lm(y[-id0] ~ X[-id0, ] - 1, data = data)
m_rm2 <- lm(y[-id2] ~ X[-id2, ] - 1, data = data)

x <- data$X[, 2]
y <- data$y

# Plot the influential sets of Alg2 and Alg0
cairo_pdf("paper/output/pov-convergence2.pdf", height = 4.2, width = 7.2)
op <- par(mar = c(2, 3, 0, .5), mfrow = c(2, 1),
  bg = "transparent", family = "Noto Sans")
# Algorithm 1
plot.new()
plot.window(xlim = range(x) + c(0, .01), ylim = range(y) + c(-0.005, 0.005))
points(x, y, cex = 1.25)
axis(1, at = round(c(min(x), max(x)), 2))
axis(2, at = round(c(min(y), 0, max(y)), 2), las = 1)
# Add regression line
abline(coef(m), col = "darkgray", lty = 1, lwd = 2)
abline(coef(m_rm0), col = "#408040", lty = 2, lwd = 2)
# Highlight sets
points(x, y, cex = 1.25)
points(x[id0], y[id0], cex = 1.5, col = "darkgray",
  pch = 21, bg = (viridisLite::inferno(n = length(id0), begin = 0.5)))
points(x[id0[1:13]], y[id0[1:13]], cex = 1.75, pch = 4)
points(x[id0[-1:-13]], y[id0[-1:-13]], cex = 1.75, pch = 3)
text(.12, -.034, labels = expression("Algorithm " * phantom("(1)")),
  col = "black", cex = 1.5, font = 1, family = "Merriweather")
text(.12, -.034, labels = expression(phantom("Algorithm ") * "(1)"),
  col = "#408040", cex = 1.5, font = 1, family = "Merriweather")
symbols(add = TRUE, 0.005, 0,
  circle = .035, inches = FALSE,
  fg = viridisLite::inferno(n = 1, begin = 0.8), lwd = 2.5)
# Algorithm 2
plot.new()
plot.window(xlim = range(x) + c(0, .01), ylim = range(y) + c(-0.005, 0.005))
points(x, y, cex = 1.25)
axis(1, at = round(c(min(x), max(x)), 2))
axis(2, at = round(c(min(y), 0, max(y)), 2), las = 1)
# Add regression line
abline(coef(m), col = "darkgray", lty = 1, lwd = 2)
abline(coef(m_rm2), col = "#008080", lty = 2, lwd = 2)
# Highlight sets
points(x, y, cex = 1.25)
points(x[id2], y[id2], cex = 1.5, col = "darkgray",
  pch = 21, bg = (viridisLite::inferno(n = length(id2), begin = 0.5)))
points(x[id2[1:13]], y[id2[1:13]], cex = 1.75, pch = 4)
points(x[id2[-1:-13]], y[id2[-1:-13]], cex = 1.75, pch = 3)
text(.12, -.034, labels = expression("Algorithm " * phantom("(2)")),
  col = "black", cex = 1.5, font = 1, family = "Merriweather")
text(.12, -.034, labels = expression(phantom("Algorithm ") * "(2)"),
  col = "#008080", cex = 1.5, font = 1, family = "Merriweather")
dev.off()


# Main table ---
tbl <- data.frame(v = c("Convergence$^\\dag$", "", "Convergence, Eastern Europe",
  "", "Thresholds$^\\dag$", "Observations", "$R^2$"))
tbl$ravallion <- c(round(coef(mdl_rav)[2], 3),
  paste0("(", round(summary(mdl_rav)$coefficients[2, 3], 2), ")"),
  "--", "--", "--[1]\\{4\\}",
  NROW(resid(mdl_rav)), round(summary(mdl_ee)$r.squared, 3))
tbl$ee <- c(round(coef(mdl_ee)[2], 3),
  paste0("(", round(summary(mdl_ee)$coefficients[2, 3], 2), ")"),
  round(coef(mdl_ee)[4], 3),
  paste0("(", round(summary(mdl_ee)$coefficients[4, 3], 2), ")"),
  "3[10]\\{24\\}",
  NROW(resid(mdl_ee)), round(summary(mdl_ee)$r.squared, 3))
tbl$cre <- c(round(coef(mdl_cre)[2], 3),
  paste0("(", round(summary(mdl_cre)$coefficients[2, 3], 2), ")"),
  "--", "--", "26[32]\\{42\\}",
  NROW(resid(mdl_cre)), round(summary(mdl_cre)$r.squared, 3))
tbl

kableExtra::kbl(tbl, format = "latex", booktabs = TRUE, align = "lccccccc",
  escape = FALSE, label = "poverty", position = "th", row.names = FALSE,
  linesep = c(rep("", 3), "\\addlinespace", rep("", 2)),
  caption = "Sensitivity of poverty convergence")
