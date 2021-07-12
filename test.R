
sapply(list.files("R/", ".R$"), function(x) source(paste0("R/", x)))

N <- 100
K <- 5
M <- 5
Z <- matrix(rnorm(N * M), N)
X <- Z %*% (b0 <- matrix(rnorm(M * K), M)) + (e0 <- rnorm(N * K))
y <- X %*% (b1 <- rnorm(K)) + (e1 <- X[, 1] * 2 * (u1 <- rnorm(N)))

sensitivity(lm(y ~ X))

require("pwt10")
require("dplyr")

# load pwt10
df <- pwt10::pwt10.0

# compute log gdppc
df$l_gdppc <- log(df$rgdpe / df$pop)

# subset to 5 year growth periods
df <- subset(df, year %in% seq(1949, 2019, 5))

# create necessary variables
df <- df %>%
  group_by(country) %>%
  mutate(lag_l_gdppc = lag(l_gdppc), g_gdppc = l_gdppc - lag(l_gdppc),
    g_pop = log(pop) - lag(log(pop)), lag_inv = lag(csh_i))

# subset and rename
df <- df[, c("country", "year", "g_gdppc", "lag_l_gdppc", "l_gdppc", "g_pop", "lag_inv")]
colnames(df) <- c("country", "year", "growth", "lag_gdppc", "l_gdppc", "growth_pop", "lag_inv")

# annual growth rate
df$growth <- df$growth / 5

# solow growth regression
summary(mdl <- lm(data = df,
  growth ~ lag_gdppc + growth_pop + lag_inv + as.factor(country) + as.factor(year)))
# on average, 3.6% conditional convergence p.a.

s1 <- sens(mdl,
  lambda = set_lambda("beta", position = 2L, sign = -1),
  options = set_options(fwl = 1L:4L, fwl_re = 10L))

s2 <- sens(mdl,
  lambda = set_lambda("beta", position = 2L, sign = -1),
  options = set_options(n_max = 100))

s1$influence[1:10, ]
s1$model[1:10, ]

s2$influence[1:10, ]
s2$model[1:10, 1:10]

plot(s1$model[, "beta-2"], type = "l", col = 1)
lines(s2$model[, "beta-2"], col = 2)
lines(cumsum(c(s1$model[1, "beta-2"],
  (-s1$influence[, "init_lambda"] - s1$model[1, "beta-2"]))),
  col = 3)


droller <- readRDS("~/repos/max_perturbation/data/droller2018.rds")

summary(mdl <- ivreg::ivreg(y ~ X - 1L | Z - 1L, data = droller$spec2,
  x = TRUE, y = TRUE))

s3 <- sens(mdl,
  lambda = set_lambda("beta", position = 2L))

s3$influence[1:10, ]
s3$model[1:10, 1:10]

s3$model <- na.omit(s3$model)

plot(s3$model[, "beta-2"], type = "l", col = 1)
lines(cumsum(c(s3$model[1, "beta-2"],
  (+s3$influence[, "init_lambda"] - s3$model[1, "beta-2"]))),
  col = 3)

summary(infl(mdl)$influence$hat)


summary(mdl <- ivreg::ivreg(y ~ X - 1L | Z - 1L, data = droller$spec2,
  x = TRUE, y = TRUE))

s3 <- sens(mdl,
  lambda = set_lambda("beta", position = 2L))

s3$influence[1:10, ]
s3$model[1:10, 1:10]

s3$model <- na.omit(s3$model)

plot(s3$model[, "beta-2"], type = "l", col = 1)
lines(cumsum(c(s3$model[1, "beta-2"],
  (+s3$influence[, "init_lambda"] - s3$model[1, "beta-2"]))),
  col = 3)
lines(influence_dfs$sign$neg$beta_est, col = 4)

summary(infl(mdl)$influence$hat)



library("zaminfluence")

reg_infl <- zaminfluence::ComputeModelInfluence(mdl)
grad_df <- GetTargetRegressorGrads(reg_infl, "Xpct_eu")
influence_dfs <- SortAndAccumulate(grad_df)
target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")
# PlotInfluence(influence_dfs$sign, "prop_removed", 0.01, target_change)

lines(influence_dfs$sign$neg$beta_est, col = 4)
