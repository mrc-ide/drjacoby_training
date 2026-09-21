# E1.Visualising_model_fit_SOLUTION.R
#
# Author: Bob Verity
# Date: 2024-06-20
#
# Inputs: Data/E1.mcmc_quadratic_regression.rds
#
# Outputs: (none)
#
# Purpose:
# "Solution" to the challenge in Extras/E1.Visualising_model_fit.R
#
# ------------------------------------------------------------------

# load packages
library(drjacoby)
library(tidyverse)

# read in quadratic fit model mcmc
mcmc <- readRDS("Data/E1.mcmc_quadratic_regression.rds")

# posterior density plots
mcmc$output |>
  filter(phase == "sampling") |>
  pivot_longer(cols = c(intercept, slope, quadratic, sigma)) |>
  ggplot(aes(x = value)) + theme_bw() +
  geom_density(fill = "blue", alpha = 0.5) +
  facet_wrap(~name, scales = "free")

# draw from the posterior
param_draws <- mcmc$output |>
  filter(phase == "sampling") |>
  select(intercept, slope, quadratic, sigma) |>
  sample_n(100) |>
  mutate(draw = row_number())

# set a continuous range of girth over which we want to predict
girth_predict <- seq(0, 25, l = 101)

# calculate our model predicted volume
trajectories <- param_draws |>
  expand_grid(girth_predict) |>
  mutate(vol_predict = intercept + slope * girth_predict + quadratic * girth_predict^2)

# plot the data and overlay our trajectories
ggplot(trees, aes(x = Girth, y = Volume)) + theme_bw() +
  geom_line(aes(x = girth_predict, y = vol_predict, group = draw),
            col = "dodgerblue", alpha = 0.2, data = trajectories) +
  geom_point()

# calculate and plot quantiles
trajectory_quantiles <- trajectories |>
  group_by(girth_predict) |>
  summarise(Q = list(quantile(vol_predict, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) |>
  unnest_wider(Q, names_sep = "_")

ggplot(trees) + theme_bw() +
  geom_ribbon(aes(x = girth_predict, ymin = `Q_2.5%`, ymax = `Q_97.5%`),
              fill = "dodgerblue", alpha = 0.2, data = trajectory_quantiles) +
  geom_ribbon(aes(x = girth_predict, ymin = `Q_25%`, ymax = `Q_75%`),
              fill = "dodgerblue", alpha = 0.2, data = trajectory_quantiles) +
  geom_point(aes(x = Girth, y = Volume))


# draw from posterior predictive
n_PPI_sims <- 100
PPI_sim_data <- trajectories |>
  group_by(girth_predict, draw) |>
  reframe(vol_predict = vol_predict,
          sigma = sigma,
          sim = 1:n_PPI_sims,
          vol_sim = rnorm(n_PPI_sims, mean = vol_predict, sd = sigma))

# get quantiles over simulated data
PPI_sim_quantiles <- PPI_sim_data |>
  group_by(girth_predict) |>
  summarise(Q = list(quantile(vol_sim, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) |>
  unnest_wider(Q, names_sep = "_")

ggplot(trees) + theme_bw() +
  geom_ribbon(aes(x = girth_predict, ymin = `Q_2.5%`, ymax = `Q_97.5%`),
              fill = "firebrick2", alpha = 0.2, data = PPI_sim_quantiles) +
  geom_ribbon(aes(x = girth_predict, ymin = `Q_25%`, ymax = `Q_75%`),
              fill = "firebrick2", alpha = 0.2, data = PPI_sim_quantiles) +
  geom_point(aes(x = Girth, y = Volume))

