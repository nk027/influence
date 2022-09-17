
# Ruggedness ---------

library("dplyr")
library("sandwich")
library("lmtest")
library("MASS") # M estimation
devtools::load_all()


# Prepare ---

# Load data and add the diamonds variable
x <- read.csv("paper/data/rugged_data.csv") |>
  mutate(diamonds = gemstones / (land_area / 100))

# Short-hand function to link observations to their ISO
get_iso <- \(m, id) {x[rownames(m$x)[id], "isocode"]}

# Short-hand to target influence on t-values
get_sens <- \(m, pos) {
  sens.lm(m,
    lambda = set_lambda("tstat", pos = pos, sign = sign(coef(m)[pos])),
    options = set_options(n_max = 80, tstat = TRUE, cluster = TRUE),
    cluster = seq_along(resid(m)))
}


# Reproduce results of Nunn and Puga (2012) -----

# Our baseline (their baseline + controls) ---
mdl_base <- lm(log(rgdppc_2000) ~ rugged * cont_africa +
  diamonds * cont_africa + soil * cont_africa +
  tropical * cont_africa + dist_coast * cont_africa, data = x)
coeftest(mdl_base, vcov = vcovHC(mdl_base, "HC1"))

# Find influential sets
sens_base <- get_sens(mdl_base, 8) # Interacted ruggedness at position 8
plot.sensitivity(sens_base)
summary.sensitivity(sens_base) # 2/5/11 (Alg0 needs 2/20/-)

# These five should be particularly bad
infl_obs <- c("SYC", "LSO", "RWA", "SWZ", "COM")
get_iso(mdl_base, sens_base$influence$id[1:5]) == infl_obs

# Compare the overlap of adaptive and initial sets in percent
overlap <- vector("numeric", 20)
for(i in seq_along(overlap)) {
  overlap[i] <- sum(sens_base$influence$id[seq(i)] %in%
    sens_base$initial$id[seq(i)]) / i
}
plot.ts(print(overlap)) # The initial set performs rather poorly

# Check their baseline specification (i.e. only coast controls) ---
mdl_plain <- lm(log(rgdppc_2000) ~ rugged * cont_africa +
  dist_coast * cont_africa, data = x)
coeftest(mdl_plain, vcov = vcovHC(mdl_plain, "HC1"))

# Find influential sets
sens_plain <- get_sens(mdl_plain, 5) # Interacted ruggedness at position 5
plot.sensitivity(sens_plain)
summary.sensitivity(sens_plain) # 2/7/16 (Alg0 needs 2/22/-)

# The top observations are similar
get_iso(mdl_plain, sens_plain$influence$id[1:7])

# Check the results with M-estimation ---
mdl_r <- MASS::rlm(log(rgdppc_2000) ~ rugged * cont_africa +
  diamonds * cont_africa + soil * cont_africa +
  tropical * cont_africa + dist_coast * cont_africa,
  data = x, method = "M", psi = psi.huber)
coeftest(mdl_r, vcov = vcovHC(mdl_r, "HC1"))

# Note that only the differential effect is significant ---
mdl_afr <- lm(log(rgdppc_2000) ~ rugged + diamonds + soil + tropical +
  dist_coast, data = x |> filter(cont_africa == 1))
coeftest(mdl_afr, vcov = vcovHC(mdl_afr, "HC1"))
mdl_row <- lm(log(rgdppc_2000) ~ rugged + diamonds + soil + tropical +
  dist_coast, data = x |> filter(cont_africa == 0))
coeftest(mdl_row, vcov = vcovHC(mdl_row, "HC1"))


# Alternative specifications -----

# Note that population in 1400 is:
#   - missing for URY and GUY
#   - zero for COM, CPV, MUS, SYC

# + log population in 1400 + 1 (and a dummy for zeros) removes significance
mdl_pop <- lm(log(rgdppc_2000) ~ rugged * cont_africa + diamonds * cont_africa +
  soil * cont_africa + tropical * cont_africa + dist_coast * cont_africa +
  log(1 + pop_1400) * cont_africa + I(pop_1400 == 0), data = x)
coeftest(mdl_pop, vcov = vcovHC(mdl_pop, "HC1"))

# + log land area removes significance
mdl_area <- lm(log(rgdppc_2000) ~ rugged * cont_africa + diamonds * cont_africa +
  soil * cont_africa + tropical * cont_africa + dist_coast * cont_africa +
  log(land_area) * cont_africa, data = x)
coeftest(mdl_area, vcov = vcovHC(mdl_area, "HC1"))

# Weighted by population in 1400 removes significance
mdl_wt <- lm(log(rgdppc_2000) ~ rugged * cont_africa +
  diamonds * cont_africa + soil * cont_africa +
  tropical * cont_africa + dist_coast * cont_africa,
  weights = pop_1400, data = x)
coeftest(mdl_wt, vcov = vcovHC(mdl_wt, "HC1"))

# Interestingly, a pooled model removes significance
mdl_pool <- lm(log(rgdppc_2000) ~ rugged + diamonds + soil +
  tropical + dist_coast, data = x)
coeftest(mdl_pool, vcov = vcovHC(mdl_pool, "HC1"))
# With either alternative (population or land area) it's back
mdl_pool_pop <- lm(log(rgdppc_2000) ~ rugged + diamonds + soil +
  tropical + dist_coast + log(1 + pop_1400) + I(pop_1400 == 0), data = x)
coeftest(mdl_pool_pop, vcov = vcovHC(mdl_pool_pop, "HC1"))

mdl_pool_area <- lm(log(rgdppc_2000) ~ rugged + diamonds + soil +
  tropical + dist_coast + log(land_area), data = x)
coeftest(mdl_pool_area, vcov = vcovHC(mdl_pool_area, "HC1"))

# Investigate the jointly influential set ---
mdl_rm5 <- lm(log(rgdppc_2000) ~ rugged * cont_africa + diamonds * cont_africa +
  soil * cont_africa + tropical * cont_africa + dist_coast * cont_africa,
  data = x |> filter(!isocode %in% infl_obs), x = TRUE, y = TRUE)
coeftest(mdl_rm5, vcov = vcovHC(mdl_rm5, "HC1"))

# Removing SYC and any one of the other four (e.g. #5) removes significance
mdl_rm2 <- lm(log(rgdppc_2000) ~ rugged * cont_africa + diamonds * cont_africa +
  soil * cont_africa + tropical * cont_africa + dist_coast * cont_africa,
  data = x |> filter(!isocode %in% infl_obs[c(1, 5)]), x = TRUE, y = TRUE)
coeftest(mdl_rm2, vcov = vcovHC(mdl_rm2, "HC1"))
# It seems fair to speak of a jointly influential set


# Prepare outputs -----

# Main table ---
tbl <- data.frame( # Extract basics
  v = c("Ruggedness, Africa$^\\dag$", "", "Ruggedness", "", "Observations", "$R^2$"),
  sapply(list(mdl_base, mdl_plain, mdl_pop, mdl_area, mdl_afr, mdl_r), \(m) {
    ct <- coeftest(m, vcov = vcovHC(m, "HC1"))
    vars <- if("rugged:cont_africa" %in% rownames(ct)) {
      ct[c("rugged:cont_africa", "rugged"), c(1, 3)]
    } else {rbind(ct["rugged", c(1, 3)], NA)}
    vars[, 1] <- round(vars[, 1], 3)
    vars[, 2] <- paste0("(", round(vars[, 2], 2), ")")
    # beta, t, beta, t, n, r²
    as.character(c(t(vars), nrow(m$model), round(summary(m)$r.squared, 3)))
  })
)

# Fill in the estimate for the separated model
vals <- coeftest(mdl_row, vcov = vcovHC(mdl_row, "HC1"))[2, c(1, 3)]
tbl[3:4, 6] <- paste0(c(round(vals[1], 3), paste0("(", round(vals[2], 2), ")")))
# Fill in R² for the robust model
tbl[6, 7] <- round(cor(mdl_r$model[, 1], mdl_r$fitted.values, method = "pearson")^2, 3)

# Add textual information
tbl <- rbind(tbl, c("Controls", rep("Yes", 6)))
tbl <- rbind(tbl, c("Population in 1400", c(rep("--", 2), "Yes", rep("--", 3))))
tbl <- rbind(tbl, c("Land area", c(rep("--", 3), "Yes", rep("--", 2))))
tbl <- rbind(tbl, c("Other controls", "Yes", "--", rep("Yes", 4)))
# These thresholds are based on the summary(sensitivity(mdl))
tbl <- rbind(tbl, c("Thresholds$^\\dag$",
  c("2[5]\\{11\\}", "2[7]\\{16\\}", "--[3]\\{6\\}", "--[4]\\{8\\}", "--[1]\\{6\\}", "--")))

# Move the textual information on controls below the estimates
tbl <- tbl[c(1:4, 7:10, 11, 5, 6), ]
# Move robust estimates, and kick the plain model (and other controls)
tbl <- tbl[-8, c(1, 2, 7, 4:6)]

kableExtra::kbl(tbl, format = "latex", booktabs = TRUE, align = "lcccccc",
  escape = FALSE, label = "rugged", position = "ht", row.names = FALSE,
  col.names = c("", "Baseline", "Robust", "Population", "Land area", "Separated"),
  linesep = c(rep("", 3), "\\addlinespace", rep("", 3), "\\addlinespace"),
  caption = "The differential effect of ruggedness in Africa")
# There's some Latex magic afterwards


# Plot -----

library("sf")
library("tmap")

# Panel A ---
data("World")
map <- rworldmap::getMap("low") |>
  st_as_sf() |> st_transform("+proj=igh") |>
  dplyr::select(isocode = ADM0_A3, name = NAME) |>
  left_join(data.frame(isocode = infl_obs, id = seq(-1, -length(infl_obs)))) |>
  left_join(x |> transmute(isocode, Ruggedness = rugged)) |>
  mutate(lab = paste0(abs(id), ". ", name), Influence = abs(id)) |>
  mutate(lab_a = ifelse(id < 0, lab, NA), lab_b = ifelse(id > 0, lab, NA))

# Viridis inferno colour scheme
cols <- viridisLite::inferno(length(infl_obs), begin = 0.5)

# Zooming in on Africa
bbox_infl <- st_bbox(map[map$id < 0 & !is.na(map$id), ]) + c(0, 0, 0, 50000)
bbox_infl <- c(2908504-5600000, -3411087-1250000, 6169455+1275000, -68337+4100000)

p <- map |>
  tm_shape(bbox = bbox_infl) +
  # Countries filled according to their ruggedness
  tm_fill("Ruggedness", palette = "Blues", style = "cont",
    legend.show = TRUE, showNA = FALSE, breaks = c(0, 7/3, 7/3*2, 7),
    labels = c("low", "", "", "high")) +
  tm_borders(col = "gray") +
  # Extra borders around Benin, Ghana, and Nigeria (see mistrust application)
  tm_shape(map |> filter(isocode %in% c("BEN", "GHA", "NGA"))) +
  tm_borders(col = "#333333", lwd = 1.5) +
  # Add lines for the coordinate system
  tm_graticules(labels.size = 1.2, col = "#333333", lwd = 1.5,
    labels.inside.frame = TRUE, y = c(0, 20), x = c(20)) +
  # Highlight the top 5 influential set
  tm_shape(map |> mutate(Influence = ifelse(id < 0, Influence, NA))) +
  tm_symbols(shape = 1, col = "Influence", border.lwd = 12,
    palette = "inferno", contrast = c(.5, 1),
    scale = 2.5, colorNA = NULL, shapeNA = NA, shape.showNA = FALSE,
    legend.col.show = FALSE) +
  tm_shape(map) + # Add labels
  tm_text("lab_a", size = 1.4, col = "#ffffff", bg.color = "#333333",
    bg.alpha = .75, auto.placement = FALSE, xmod = -4, ymod = 1) +
  # Add legend for the influence
  tm_add_legend(type = "fill", size = 3, title = "Influence",
    col = viridisLite::inferno(length(infl_obs), begin = .5),
    labels = c("highest", "2nd", "3rd", "4th", "5th")) +
  tm_layout(frame = "gray", fontfamily = "Noto Sans",
    title = "Two nations drive the blessing of bad geography",
    title.color = "#ffffff", title.bg.color = "#333333", title.bg.alpha = .75,
    title.position = c("right", "top"), title.size = 2,
    title.fontfamily = "Merriweather",
    outer.bg.color = "transparent", bg.color = "white",
    inner.margins = 0, outer.margins = c(0, 0, 0, 0),
    legend.position = c(0, 0), legend.frame = FALSE,
    legend.text.size = 1.2, legend.title.size = 1.4,
    legend.title.fontface = "bold")

tmap_save(p, "paper/output/ruggedness/map.pdf", units = "in", device = cairo_pdf,
  width = 7, height = 6, bg = "transparent")
# There's a lot of Inkscape magic afterwards

# Panel B ---
iso_a2 <- get_iso(mdl_base, sens_base$influence$id)[seq(11)]
iso_a0 <- get_iso(mdl_base, sens_base$initial$id)[seq(11)]
# Get the values
tvals_a2 <- c(coeftest(mdl_base, vcov = vcovHC(mdl_base, "HC1"))[8, 3],
  sens_base$influence$lambda[1:11])
tvals_a0 <- init.sensitivity(sens_base)$initial[1:12]

# We just plot points at the relevant values and use that to insert text later
cairo_ps("paper/output/ruggedness/stack.ps", height = 10, width = 5, pointsize = 12)
op <- par(mar = c(1, 4, 2, 1))
plot.new()
plot.window(xlim = c(0, 3), ylim = c(-2.75, 2.75))
axis(3, at = c(1, 2), labels = c("A2", "A0"), lwd = 0, lwd.ticks = 1)
axis(2, lwd = 1, lwd.ticks = 1, las = 1, col = "gray",
  at = c(max(tvals_a2), 1, 0, -1, min(tvals_a2), -Inf),
  labels = c(sprintf("%0.2f", c(max(tvals_a2), 1, 0, -1, min(tvals_a2))), ""))
abline(h = c(max(tvals_a2), 1, 0, -1, min(tvals_a2)), lty = 3, col = "gray")
text(x = 0.5, y = tvals_a2[-12], adj = c(0.5, 1),
  label = iso_a2, cex = 1, family = "Inconsolata")
points(x = rep(1, 12), y = tvals_a2, pch = 20, cex = 0.01)
text(x = 2.5, y = tvals_a0[-12], adj = c(0.5, 1),
  label = iso_a0, cex = 1, family = "Inconsolata")
points(x = rep(2, 12), y = tvals_a0, pch = 20, cex = 0.01)
dev.off()
# Magic in Inkscape follows
