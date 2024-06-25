# 03b.Multimodal_example.R
#
# Author: Bob Verity
# Date: 2024-06-25
#
# Inputs:
# Utilities/gaussian_process_utilities.R
# Data/rainfall.rds
#
# Outputs: (none)
#
# Purpose:
# This is an example of a model where the posterior distribution is multi-modal
# (contains multiple peaks). If we're not careful then our MCMC can mix very
# badly in these situations, or in extreme cases can fail to explore some peaks
# entirely. Your job is to read through the script and see if you can get the
# MCMC to mix.
#
# Our data consists of (made up) rainfall data, measured every 2 weeks for a
# period of 2 years. The pattern in the data is clearly highly seasonal, which
# is typical of some regions of e.g. northern Africa.
#
# Our plan is to model this data using a Gaussian process (GP). This method is useful
# for data that changes smoothly over time or space. You'll see
# in the likelihood function that we start by making a Gaussian kernel, that
# describes the autocorrelation in the data in terms of three parameters:
#
# 1) length_scale: The length_scale parameter determines how quickly the
#   correlation between points decreases with distance. A small length_scale means
#   that the function varies rapidly, while a large length_scale indicates that
#   the function varies slowly.
# 2) sigma_f: sigma_f represents the variance of the underlying function, or the
#   amplitude. It defines the average vertical variation of the function from its
#   mean. A larger sigma_f means that the function values can deviate more from
#   the mean.
# 3) sigma_N: sigma_n denotes the noise level in the observations (the standard
#   deviation of the noise). It accounts for the uncertainty or measurement noise
#   in the observed data. A larger sigma_n indicates more noise in the
#   observations, and the GP will be less confident about fitting the data
#   exactly.
#
# Once we have the Gaussian kernel, our last free parameter is mu, which defines
# the global mean. We will fit all four of these parameters to data.
# 
# ------------------------------------------------------------------
# load packages
#install.packages("mvtnorm")
library(mvtnorm)
#install.packages("GGally")
library(GGally)
library(tidyverse)

# we will make use of a few utility functions. You don't need to look deeply
# into these functions (although you can)! They allow us to do things like draw
# from the posterior mean under a Gaussian process
source("Utilities/gaussian_process_utilities.R")

# --------------------------

# read in rainfall data
rainfall_data <- readRDS("Data/rainfall.rds")

# plot rainfall data
ggplot(rainfall_data) + theme_bw() +
  geom_point(aes(x = time, y = rainfall)) +
  ylim(0, 130) + xlab("Time (months)") + ylab("Rainfall (mm)")

# define parameters
df_params <- define_params(name = "length_scale", min = 0, max = 10,
                           name = "sigma_f", min = 0, max = 100,
                           name = "sigma_n", min = 0, max = 100,
                           name = "mu", min = 0, max = 100)

# define the log-likelihood function using the Gaussian kernel and the dmvnorm
# function
loglike <- function(params, data, misc) {
  
  # unpack parameters
  length_scale <- params["length_scale"]
  sigma_f <- params["sigma_f"]
  sigma_n <- params["sigma_n"]
  mu <- params["mu"]
  
  n <- nrow(data)
  K <- Gaussian_kernel(data$time, data$time, length_scale, sigma_f) + sigma_n^2 * diag(n)
  dmvnorm(data$rainfall, mean = rep(mu, n), sigma = K, log = TRUE)
}

# define the logprior function to be uniform over all model parameters
logprior <- function(params, misc) {
  dunif(params["length_scale"], 0, 10, log = TRUE) + 
    dunif(params["sigma_f"], 0, 100, log = TRUE) + 
    dunif(params["sigma_n"], 0, 100, log = TRUE) +
    dunif(params["mu"], 0, 100, log = TRUE)
}

# run the MCMC
mcmc <- run_mcmc(data = dat,
                 df_params = df_params,
                 loglike = loglike,
                 logprior = logprior,
                 burnin = 1e2,
                 samples = 1e3,
                 chains = 5)

# oh no, based on the ESS our MCMC doesn't appear to have mixed well!
mcmc$diagnostics

# ----------------------------
#### CHALLENGES!
# - how can we detect poor mixing in other ways than the ESS?
# - can you get this MCMC to mix? There are multiple ways of doing this, can you think of more than one?
# - what is going on here? Why are we getting multiple peaks? (hint, check out the model fits code below)
# - although the GP is a reasonable modelling choice here, what are some other approaches that might be more suitable?

# ------------------------------------------------------------------
# Check model fits
# NB, you can use the code below to check the model fits to the data

# draw some parameters from the posterior
param_draws <- mcmc$output |>
  filter(phase == "sampling") |>
  select(length_scale, sigma_f, sigma_n, mu) |>
  sample_n(100) |>
  mutate(draw = row_number())

# set the time period that we want to predict over
time_pred <- seq(0, 25, l = 101)

# draw a series of trajectories over this time period. Each of these represents
# a draw from the posterior mean under our fitted model
trajectories <- param_draws |>
  group_by(draw) |>
  reframe(time_pred = time_pred,
          mean_draw = draw_posterior_mean(x = dat$time,
                                          y = dat$rainfall,
                                          x_pred = time_pred,
                                          length_scale = length_scale,
                                          sigma_f = sigma_f,
                                          sigma_n = sigma_n,
                                          mu = mu))

# plot trajectories against the data
trajectories |>
  ggplot() + theme_bw() +
  geom_line(aes(x = time_pred, y = mean_draw, group = draw),
            col = "dodgerblue", alpha = 0.25) +
  ylim(c(0.5*min(dat$rainfall), 1.5*max(dat$rainfall))) +
  geom_point(aes(x = time, y = rainfall), data = dat)
