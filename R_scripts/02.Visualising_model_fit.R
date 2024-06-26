# 02.Visualising_model_fit.R
#
# Author: Bob Verity
# Date: 2024-06-20
#
# Inputs: Data/02.mcmc_linear_regression.rds
#         Data/02.mcmc_quadratic_regression.rds
#
# Outputs: (none)
#
# Purpose:
# Demonstrates ways that we can explore the fit of a model to data, and how this
# differs from comparing between models.
#
# ------------------------------------------------------------------

# load packages
library(tidyverse)

# ------------------------------------------------------------------
### PART1 - Credible intervals

# start by loading in a successful MCMC run from the previous script on linear
# modelling of the trees dataset
mcmc <- readRDS("Data/02.mcmc_linear_regression.rds")

# we can make simple density plots of our parameters, but what does this
# actually tell us about our model fit? Not much really
plot_density(mcmc)

# how can we visualise our model fit against our data? Remember that we assumed
# a linear model for the trees data:
# 
# Volume = intercept + slope*Girth + epsilon
#
# where epsilon was error (Gaussian noise) with standard deviation sigma. We can
# divide this into two parts: the linear model and the error. Looking at just
# the linear model we have:
#
# Volume_expected = intercept + slope*Girth
#
# This is the part of the model that describes the relationship we are trying to
# fit. We can construct a 95% CrI on Volume_expected to explore this
# relationship. Note that this is not one of the raw parameters of our model,
# i.e. it is not the slope, intercept, or sigma. Rather, it is a compound
# parameter that is made up of both the intercept and the slope.

# lets start by selecting 100 parameter draws from the posterior
param_draws <- sample_chains(mcmc, 100)

# set a continuous range of girth over which we want to predict
girth_predict <- seq(0, 25, l = 101)

# we can calculate our model predicted volume at each one of these girth points,
# and for every one of our posterior draws
trajectories <- param_draws |>
  expand_grid(girth_predict) |>
  mutate(vol_predict = intercept + slope * girth_predict)

# note that this trajectories object has a row for every value of girth_predict
# and every one of our posterior draws. It can get quite large if we use a large
# number of draws
dim(trajectories)

# let's plot the data and overlay our trajectories
ggplot(trees, aes(x = Girth, y = Volume)) + theme_bw() +
  geom_line(aes(x = girth_predict, y = vol_predict, group = sample),
            col = "dodgerblue", alpha = 0.2, data = trajectories) +
  geom_point()

# the first thing to note about this model fit is that it's absolute trash. A
# linear model is not a good choice for these data, which clearly have some sort
# of curve. Also our limitation that the intercept must be greater than zero has
# constrained the model so it barely captures the trend. Don't worry, we'll have
# a lovely example of a quadratic model fit shortly.

# another way to summarise this same information is to construct quantiles over
# the trajectories. Let's do that by first grouping by the girth_predict
# variable and then using the quantile function. We will explore the 95% CrI
# (going from 2.5% to 97.5%) and also the 50% CrI (going from 25% to 75%).
#
# To make tidyr happy we have to summarise as a list and then un-nest, but
# really we're not doing anything fancy here, and the same thing could be
# achieved through a loop
trajectory_quantiles <- trajectories |>
  group_by(girth_predict) |>
  summarise(Q = list(quantile(vol_predict, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) |>
  unnest_wider(Q, names_sep = "_")

# this is a much more compact representation of the fit
dim(trajectory_quantiles)

# let's plot against the data
ggplot(trees) + theme_bw() +
  geom_ribbon(aes(x = girth_predict, ymin = `Q_2.5%`, ymax = `Q_97.5%`),
              fill = "dodgerblue", alpha = 0.2, data = trajectory_quantiles) +
  geom_ribbon(aes(x = girth_predict, ymin = `Q_25%`, ymax = `Q_75%`),
              fill = "dodgerblue", alpha = 0.2, data = trajectory_quantiles) +
  geom_point(aes(x = Girth, y = Volume))

# QUESTION!
# what are some of the advantages and disadvantages of plotting trajectories vs.
# quantiles?

# ------------------------------------------------------------------
### PART2 - Predictive intervals

# the plots above show us our parameter uncertainty in a more meaningful way
# than simple posterior distributions of the raw parameters. However, they do
# not tell us how well we capture the spread in the data. For this, we need
# posterior predictive intervals. This is the interval that we expect *new* data
# to fall within, if we were to observe it.

# one way to construct this interval is to actually simulate new data from the
# model. We need to do this over the entire interval that we are interested in,
# so over every value of girth_predict.
# 
# We also want to do this over every one of our parameter draws. Doing so means
# we "marginalise out" these parameters, in other words we consider what we do
# and don't know about these parameters when making predictions.
n_PPI_sims <- 10
PPI_sim_data <- trajectories |>
  group_by(girth_predict, sample) |>
  reframe(sim = 1:n_PPI_sims,
          vol_predict = vol_predict,
          sigma = sigma,
          vol_sim = rnorm(n_PPI_sims, mean = vol_predict, sd = sigma))

# then we get quantiles over the simulated data, but *not* grouping by sample
PPI_sim_quantiles <- PPI_sim_data |>
  group_by(girth_predict) |>
  summarise(Q = list(quantile(vol_sim, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))) |>
  unnest_wider(Q, names_sep = "_")

# we can now plot this as a predictive interval
ggplot(trees) + theme_bw() +
  geom_ribbon(aes(x = girth_predict, ymin = `Q_2.5%`, ymax = `Q_97.5%`),
              fill = "firebrick2", alpha = 0.2, data = PPI_sim_quantiles) +
  geom_ribbon(aes(x = girth_predict, ymin = `Q_25%`, ymax = `Q_75%`),
              fill = "firebrick2", alpha = 0.2, data = PPI_sim_quantiles) +
  geom_point(aes(x = Girth, y = Volume))

# ----------------------------
#### CHALLENGES!
# You know when I said I have a lovely example of a quadratic model fit? The bad
# news is, I lied. The good news is, you get the pleasure of doing it yourself!
# I already ran the MCMC for you so you can just import the mcmc object
mcmc <- readRDS("Data/02.mcmc_quadratic_regression.rds")

# you'll see that this output has columns for intercept, slope, and quadratic.
head(mcmc$output)

# The formula assumed in this model is:
#
# Volume_expected = intercept + slope*Girth + quadratic*Girth^2
#
# your task is to produce the same plots as above, looking at trajectories,
# credible intervals, and predictive intervals.

# Visually, does this model fit any better? Can you think of any tweaks or
# improvements?

