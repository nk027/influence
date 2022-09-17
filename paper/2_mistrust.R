
# Mistrust ---------

library("dplyr")
library("sandwich")
library("lmtest")
devtools::load_all()


# Prepare ---

# Load data and add a dummy to separate Eastern and Western Africa
x <- foreign::read.dta("paper/data/trust_data.dta") |>
  mutate(west = centroid_long < 20)

# Build the datasets and prepare the clustering for both dependents
variables <- c("ln_export_area", "age", "age2", "male", "urban_dum",
  "education", "occupation", "religion", "living_conditions",
  "district_ethnic_frac", "frac_ethnicity_in_district", "isocode")
# Trust of relatives ~
id_rel <- !apply(x[, c("trust_relatives", variables)], 1, \(x) any(is.na(x)))
x_rel <- x[id_rel, ]
cluster_rel <- x[id_rel, c("murdock_name", "district")]
# Trust of neighbours ~
id_nei <- !apply(x[, c("trust_neighbors", variables)], 1, \(x) any(is.na(x)))
x_nei <- x[id_nei, ]
cluster_nei <- x[id_nei, c("murdock_name", "district")]

# Short-hand function to link observations to their countries
get_nation <- \(m, s, max = 1000, each = 100) {
  m$model$id <- NA_integer_
  m$model$id[s$influence$id[seq(max)]] <- rep(seq(max / each), each = each)
  m$model$iso <- as.character(m$model[["factor(isocode)"]])
  isos <- m$model |> filter(!is.na(id)) |> group_by(iso) |>
    count() |> arrange(desc(n))
  rank <- m$model |> filter(!is.na(id)) |> group_by(iso, id) |>
    count() |> tidyr::pivot_wider(values_from = "n", names_from = "id")
  isos |> left_join(rank)
}

# Short-hand to target influence on t-values
get_sens <- \(m, cluster, pos, fwl = 1:6) {
  sens.lm(m,
    lambda = set_lambda("tstat", pos = pos, sign = sign(coef(m)[pos])),
    options = set_options(n_max = 1200, tstat = TRUE, fwl = fwl, fwl_re = 1),
    cluster = cluster)
}


# Reproduce results of Nunn and Wantchekon (2011) -----

# Trust of relatives ---
mdl_rel <- lm(trust_relatives ~
  ln_export_area + age + age2 + male + urban_dum + factor(education) +
  factor(occupation) + factor(religion) + factor(living_conditions) +
  factor(isocode) + district_ethnic_frac + frac_ethnicity_in_district,
  data = x_rel)
coeftest(mdl_rel, vcov = vcovCL(mdl_rel, "HC1", cluster = cluster_rel))
# -0.133 (0.03609) vs. -0.133 (0.037) --- there may be a rounding error

rmdl_rel <- MASS::rlm(trust_relatives ~
  ln_export_area + age + age2 + male + urban_dum + factor(education) +
  factor(occupation) + factor(religion) + factor(living_conditions) +
  factor(isocode) + district_ethnic_frac + frac_ethnicity_in_district,
  data = x_rel)
coeftest(rmdl_rel, vcov = vcovCL(rmdl_rel, "HC1", cluster = cluster_rel))

# Trust of Neighbours ---
mdl_nei <- lm(trust_neighbors ~
  ln_export_area + age + age2 + male + urban_dum + factor(education) +
  factor(occupation) + factor(religion) + factor(living_conditions) +
  factor(isocode) + district_ethnic_frac + frac_ethnicity_in_district,
  data = x_nei)
coeftest(mdl_nei, vcov = vcovCL(mdl_nei, "HC1", cluster = cluster_nei))
# -0.159 (0.03414) vs. -0.159 (0.034)

rmdl_nei <- MASS::rlm(trust_neighbors ~
  ln_export_area + age + age2 + male + urban_dum + factor(education) +
  factor(occupation) + factor(religion) + factor(living_conditions) +
  factor(isocode) + district_ethnic_frac + frac_ethnicity_in_district,
  data = x_nei)
coeftest(rmdl_nei, vcov = vcovCL(rmdl_nei, "HC1", cluster = cluster_nei))

# Find influential sets (runs a while) ---
sens_rel <- get_sens(mdl_rel, cluster_rel, 2)
sens_nei <- get_sens(mdl_nei, cluster_nei, 2)

# Note that Alg0 doesn't use the clustered standard errors
op <- par(mfrow = c(1, 2))
plot.sensitivity(sens_rel)
plot.sensitivity(sens_nei)
par(op)
summary.sensitivity(sens_rel, threshold = qnorm(.995)) # 105/380/656
summary.sensitivity(sens_nei, threshold = qnorm(.995)) # 161/425/768
# Elements of the sets are mostly located in BEN, NGA, GHA
get_nation(mdl_rel, sens_rel)
get_nation(mdl_nei, sens_nei)


# Alternative specifications -----

# Use a single interaction term ---

mdl_rel_i <- lm(trust_relatives ~
  ln_export_area * west + age + age2 + male + urban_dum + factor(education) +
  factor(occupation) + factor(religion) + factor(living_conditions) +
  factor(isocode) + district_ethnic_frac + frac_ethnicity_in_district,
  data = x_rel)
coeftest(mdl_rel_i, vcov = vcovCL(mdl_rel_i, "HC1", cluster = cluster_rel))
# Only the west is impacted

mdl_nei_i <- lm(trust_neighbors ~
  ln_export_area * west + age + age2 + male + urban_dum + factor(education) +
  factor(occupation) + factor(religion) + factor(living_conditions) +
  factor(isocode) + district_ethnic_frac + frac_ethnicity_in_district,
  data = x_nei)
coeftest(mdl_nei_i, vcov = vcovCL(mdl_nei_i, "HC1", cluster = cluster_nei))
# Only the west is impacted

# Sensitivity of these results ---
sens_rel_i <- get_sens(mdl_rel_i, cluster_rel, 8, c(1:7, 80))
sens_nei_i <- get_sens(mdl_nei_i, cluster_nei, 8, c(1:7, 80))
summary.sensitivity(sens_rel_i)
summary.sensitivity(sens_nei_i)
# Here we get MWI, MOZ, and some TZA --- East and West move closer together
get_nation(mdl_rel_i, sens_rel_i, threshold = qnorm(.995))
get_nation(mdl_nei_i, sens_nel_i, threshold = qnorm(.995))

# Separate models for Eastern and Western Africa ---
mdl_rel_w <- lm(trust_relatives ~
  ln_export_area + age + age2 + male + urban_dum + factor(education) +
  factor(occupation) + factor(religion) + factor(living_conditions) +
  factor(isocode) + district_ethnic_frac + frac_ethnicity_in_district,
  data = x_rel |> filter(west))
coeftest(mdl_rel_w, vcov = vcovCL(mdl_rel_w, "HC1", cluster = cluster_rel[x_rel$west, ]))
# -0.145 (0.038)
mdl_rel_e <- lm(trust_relatives ~
  ln_export_area + age + age2 + male + urban_dum + factor(education) +
  factor(occupation) + factor(religion) + factor(living_conditions) +
  factor(isocode) + district_ethnic_frac + frac_ethnicity_in_district,
  data = x_rel |> filter(!west))
coeftest(mdl_rel_e, vcov = vcovCL(mdl_rel_e, "HC1", cluster = cluster_rel[!x_rel$west, ]))
# 0.053 (0.055)
mdl_nei_w <- lm(trust_neighbors ~
  ln_export_area + age + age2 + male + urban_dum + factor(education) +
  factor(occupation) + factor(religion) + factor(living_conditions) +
  factor(isocode) + district_ethnic_frac + frac_ethnicity_in_district,
  data = x_nei |> filter(west))
coeftest(mdl_nei_w, vcov = vcovCL(mdl_nei_w, "HC1", cluster = cluster_nei[x_nei$west, ]))
# -0.168 (0.038)
mdl_nei_e <- lm(trust_neighbors ~
  ln_export_area + age + age2 + male + urban_dum + factor(education) +
  factor(occupation) + factor(religion) + factor(living_conditions) +
  factor(isocode) + district_ethnic_frac + frac_ethnicity_in_district,
  data = x_nei |> filter(!west))
coeftest(mdl_nei_e, vcov = vcovCL(mdl_nei_e, "HC1", cluster = cluster_nei[!x_nei$west, ]))
# 0.023 (0.073)

# Sensitivity of these results ---
sens_rel_w <- get_sens(mdl_rel_w, cluster_rel[x_rel$west, ], 2)
sens_nei_w <- get_sens(mdl_nei_w, cluster_nei[x_nei$west, ], 2)
summary.sensitivity(sens_rel_w, threshold = qnorm(.995)) # 78/301/532
summary.sensitivity(sens_nei_w, threshold = qnorm(.995)) # 133/323/527


# Prepare outputs -----

# Table with main results ---
tbl <- data.frame(
  v = c("Exports/Area$^\\dag$", "", "Exports/Area, East", "",
    "Individual controls", "District controls", "Country fixed effects",
    "Thresholds$^\\dag$", "Observations", "Ethnicity clusters", "District clusters", "$R^2$"),
  relatives_re = c(
    coeftest(mdl_rel, vcov = vcovCL(mdl_rel,
      "HC1", cluster = cluster_rel))[2, c(1, 3)] |>
      {\(x) c(round(x[1], 3), paste0("(", round(x[2], 2), ")"))}(), "", "",
    rep("Yes", 3),
    "105[380]\\{656\\}", nrow(x_rel),
    length(unique(cluster_rel[, 1])), length(unique(cluster_rel[, 2])),
    round(summary(mdl_rel)$r.squared, 3)
  ),
  relatives_split = c(
    coeftest(mdl_rel_w, vcov = vcovCL(mdl_rel_w,
      "HC1", cluster = cluster_rel[x_rel$west, ]))[2, c(1, 3)] |>
      {\(x) c(round(x[1], 3), paste0("(", round(x[2], 2), ")"))}(),
    coeftest(mdl_rel_e, vcov = vcovCL(mdl_rel_e,
      "HC1", cluster = cluster_rel[!x_rel$west, ]))[2, c(1, 3)] |>
      {\(x) c(round(x[1], 3), paste0("(", round(x[2], 2), ")"))}(),
    rep("Yes", 3),
    "78[301]\\{532\\}", paste0(sum(x_rel$west), "|", sum(!x_rel$west), ""),
    paste0(length(unique(cluster_rel[x_rel$west, 1])), "|",
      length(unique(cluster_rel[!x_rel$west, 1])), ""),
    paste0(length(unique(cluster_rel[x_rel$west, 2])), "|",
      length(unique(cluster_rel[!x_rel$west, 2])), ""),
    paste0(round(summary(mdl_rel_w)$r.squared, 3), "|",
      round(summary(mdl_rel_e)$r.squared, 3), "")
  ),
  neighbours_re = c(
    coeftest(mdl_nei, vcov = vcovCL(mdl_nei,
      "HC1", cluster = cluster_nei))[2, c(1, 3)] |>
      {\(x) c(round(x[1], 3), paste0("(", round(x[2], 2), ")"))}(), "", "",
    rep("Yes", 3),
    "161[425]\\{768\\}", nrow(x_nei),
    length(unique(cluster_nei[, 1])), length(unique(cluster_nei[, 2])),
    round(summary(mdl_nei)$r.squared, 3)
  ),
  neighbours_split = c(
    coeftest(mdl_nei_w, vcov = vcovCL(mdl_nei_w,
      "HC1", cluster = cluster_nei[x_nei$west, ]))[2, c(1, 3)] |>
      {\(x) c(round(x[1], 3), paste0("(", round(x[2], 2), ")"))}(),
    coeftest(mdl_nei_e, vcov = vcovCL(mdl_nei_e,
      "HC1", cluster = cluster_nei[!x_nei$west, ]))[2, c(1, 3)] |>
      {\(x) c(round(x[1], 3), paste0("(", round(x[2], 2), ")"))}(),
    rep("Yes", 3),
    "133[323]\\{527\\}", paste0(sum(x_nei$west), "|", sum(!x_nei$west), ""),
    paste0(length(unique(cluster_nei[x_nei$west, 1])), "|",
      length(unique(cluster_nei[!x_nei$west, 1])), ""),
    paste0(length(unique(cluster_nei[x_nei$west, 2])), "|",
      length(unique(cluster_nei[!x_nei$west, 2])), ""),
    paste0(round(summary(mdl_nei_w)$r.squared, 3), "|",
      round(summary(mdl_nei_e)$r.squared, 3), "")
  )
)

# Add a header differentiating pooled and split models
tbl <- rbind(c("", rep(c("Pooled", "West|East"), 2)), tbl)

kableExtra::kbl(tbl, format = "latex", booktabs = TRUE, align = "lcccc",
  escape = FALSE, label = "mistrust", position = "ht", row.names = FALSE,
  col.names = c("", "trust of relatives $\\sim$", "", "trust of neighbours $\\sim$", ""),
  linesep = c(rep("", 4), "\\addlinespace", rep("", 2), "\\addlinespace", rep("", 3)),
  caption = "The origins of mistrust")


# Information on removals -----

library("ggplot2")

# Stacked barplot ---
n_max <- 600
n_split <- 100
df <- get_nation(mdl_rel, sens_rel, n_max, n_split) |>
  tidyr::pivot_longer(cols = -1:-2) |>
  mutate(value = ifelse(is.na(value), 0, value))

for(i in seq(2, n_max / n_split)) {
  df[df$name == i, "value"] <- df[df$name == i, "value"] +
    df[df$name == (i - 1), "value"]
}

# Group together all but the Top 3
df |> mutate(name = paste0("top ", as.integer(name) * n_split)) |>
  mutate(name = factor(name,
    levels = paste0("top ", seq(n_split, n_max, n_split)))) |>
  mutate(iso = ifelse(iso %in% c("BEN", "NGA", "GHA"), iso, "OTH")) |>
  mutate(iso = factor(iso, levels = rev(c("BEN", "NGA", "GHA", "OTH",
    "MWI", "KEN", "MLI", "MOZ", "SEN", "TZA", "UGA", "ZMB")))) |>
  ggplot() +
  geom_bar(stat = "identity", aes(x = name, fill = iso, y = value)) +
  scale_fill_viridis_d()

# Simple table ---
tbl <- get_nation(mdl_rel, sens_rel, 600, 100) # Top 600, binned in hundreds
t(rbind(tbl[1:3, ], "Other" = colSums(tbl[-1:-3, -1], na.rm = TRUE))) |>
  kableExtra::kbl(format = "latex", booktabs = TRUE, align = "lrrrr",
    escape = FALSE, label = "mistrust", position = "ht", row.names = TRUE,
    caption = "Distribution of the top 600 most influential observations on trust in relatives")
