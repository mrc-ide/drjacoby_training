# 01.Basic_linear_regression.R
#
# Author: Bob Verity
# Date: 2024-06-20
#
# Inputs: (none)
#
# Outputs: (none)
#
# Purpose:
# This first example demonstrates how to correctly format the inputs that
# drjacoby expects. We will work with the "trees" dataset that comes inbuilt
# with R. Our objective will be to fit a model to predict the volume of a tree
# based on the girth. Once complete, we will look at the MCMC output format
# produced by drjacoby, and explore some plotting functions and MCMC
# diagnostics. Finally, we will plot the fitted model against the data to see
# whether it is actually any good.
#
# PLEASE MESS WITH THIS FILE!
# The best way to get used to drjacoby is to fiddle with it. Try changing
# things, even try breaking things if you want. Try copying code into another
# script and running on your own data, don't feel limited to this example.
#
# ------------------------------------------------------------------

# the simplest way to install drjacoby is via the runiverse. Choose "Y" when
# asked to install from source.
#install.packages('drjacoby', repos = 'https://mrc-ide.r-universe.dev')

# if that fails, you should follow the installation instructions here:
# https://mrc-ide.github.io/drjacoby/articles/installation.html

#devtools::install_github("mrc-ide/drjacoby")
library(drjacoby)

# load other packages
library(tidyverse)

# ----------------------------
#### STEP 1 - data
# format your data as either a named list or a data.frame (which is a special
# type of list)
trees
is.list(trees)

# ----------------------------
#### STEP 2 - parameters
# define parameters data.frame. At a minimum, this needs to have a "name", "min"
# and "max" column. You can make this any way you like, but we have a simple
# utility function in the package to make life easier:
df_params <- define_params(name = "slope", min = -10, max = 10,
                           name = "intercept", min = 0, max = 100,
                           name = "sigma", min = 0, max = Inf)
df_params

# ----------------------------
#### STEP 3 - loglikelihood and logprior functions
# these functions must follow this exact convention in terms of input arguments.
# Inside the function you can do whatever you want. It should return the result
# IN LOG SPACE!

# loglikelihood function
loglike <- function(params, data, misc) {
  
  # unpack parameters
  slope <- params["slope"]
  intercept <- params["intercept"]
  sigma <- params["sigma"]
  
  # our linear model of volume as a function of girth
  vol_predict <- intercept + slope * data$Girth
  
  # get probability of each data point (aka the likelihood) IN LOG SPACE
  loglike_point <- dnorm(data$Volume, mean = vol_predict, sd = sigma, log = TRUE)
  
  # sum these up to get the overall log-likelihood and return
  sum(loglike_point)
}

# logprior function
logprior <- function(params, misc) {
  
  # unpack parameters
  slope <- params["slope"]
  intercept <- params["intercept"]
  sigma <- params["sigma"]
  
  # apply priors IN LOG SPACE
  dunif(slope, -10, 10, log = TRUE) + 
    dunif(intercept, 0, 100, log = TRUE) + 
    dlnorm(sigma, meanlog = 0, sdlog = 1, log = TRUE)
}

# ----------------------------
#### STEP 4 - run the MCMC
# bring in the objects and functions defined above
mcmc <- run_mcmc(data = trees,
                 df_params = df_params,
                 loglike = loglike,
                 logprior = logprior,
                 burnin = 5e2,
                 samples = 1e3,
                 chains = 5)

# ----------------------------
#### STEP 5 - explore output and check behaviour

# printing the raw object gives an overview
mcmc

# looking under the hood, the output is a list with four elements
names(mcmc)

# the parameter draws are in the output element
head(mcmc$output)

# things to do with parallel tempering are in the pt element. We can ignore this
# most of the time
head(mcmc$pt)

# MCMC diagnostics are in the diagnostics element

# the rhat element applies to the burn-in phase and tells us how close each
# parameter is to convergence. As a rule of thumb we aim for values below 1.1.
# Note that you will only get this output when running multiple chains.

# the ess element applies to the sampling phase and tells us how many
# effectively independent samples we have from the posterior. There are no
# strict rules on what is "enough" samples, but we want to produce smooth
# posterior distributions, precise summary statistics, and we want to be
# confident we've explored all the major peaks in the posterior. I would be
# ashamed of anything less than 100.
mcmc$diagnostics

# the parameters element stores a record of everything that went into the MCMC.
# This includes keeping a copy of the data - in some cases we might not want
# this, in which case you should run with run_mcmc(save_data = FALSE)
mcmc$parameters

# in addition to looking at the diagnostics, we can get a feel for how the MCMC
# behaved by looking at plots

# we can browse through the trace plots for each parameter in turn. This plots
# the sampling phase only
plot_trace(mcmc)

# we can also be more specific with this plot, for example setting the
# parameter(s) we are interested in and which phases to plot. Let's plot both
# the burn-in and sampling phases together
plot_trace(mcmc, show = "slope", phase = "both")

# what we are looking for:
# - chains started far apart but converged over time. This is why it is important
#   to use multiple chains (5 or more) and to deliberately start them in diverse
#   locations. If we start them all within the peak then it's harder to see convergence
# - convergence happened within our burn-in iterations, and did not creep into the
#   sampling phase
# - trace plot looks like a "hairy caterpillar". No long stretches where chains disagree
# - lag plot shows how many iterations apart samples can be considered independent
#   (when autocorrelation falls to around zero). Ideally we want a fast fall-off in
#   autocorrelation, but it's OK if you don't get this as long as you can make up for it
#   with enough iterations. The ESS is your best guide here as it takes autocorrelation
#   into account

# we can get a quick summary over multiple parameters by plotting the 95% CrI
plot_credible(mcmc) +
  coord_flip()

# we can explore the relationship between two parameters via a bivariate
# scatterplot
plot_scatter(mcmc, "slope", "intercept")

# we can also produce scatterplots of every pair of parameters
plot_pairs(mcmc)

# we may want to know which parameters are strongly correlated. We can explore
# this by plotting the correlation matrix
plot_cor_mat(mcmc)

# ----------------------------
#### STEP 6 - does the model actually fit the data?
# everything in step 5 tells us whether the MCMC behaved. None of it tells us
# whether the model is any good! For that we need to look at the fitted model
# against the data.

# start by taking 100 parameter draws from the posterior
param_draws <- sample_chains(mcmc, 100)
head(param_draws)

# each of these draws defines a straight line. Let's calculate the predicted
# volume along each line, over a continuous range of girths
girth_predict <- seq(0, 25, l = 101)

trajectories <- param_draws |>
  expand_grid(girth_predict) |>
  mutate(vol_predict = intercept + slope * girth_predict)

# now overlay these lines on the data. The spread of the lines shows our
# posterior uncertainty in the fitted relationship
ggplot(trees, aes(x = Girth, y = Volume)) + theme_bw() +
  geom_line(aes(x = girth_predict, y = vol_predict, group = sample),
            col = "dodgerblue", alpha = 0.2, data = trajectories) +
  geom_point()

# the first thing to note about this model fit is that it's absolute trash. A
# linear model is not a good choice for these data, which clearly have some sort
# of curve. Also our limitation that the intercept must be greater than zero has
# constrained the model so it barely captures the trend.
#
# this is worth dwelling on: our MCMC diagnostics were all excellent! Good
# diagnostics tell you that you have found the posterior. They tell you nothing
# about whether the model describes the data. Always plot the fit.
#
# the recipe we used here (draw from the posterior, push each draw through the
# model, plot against the data) works for any model, and we will use it again
# later on an epidemic model. If you want to take it further, the Extras folder
# has self-study scripts on credible intervals, posterior predictive intervals
# and posterior predictive checks.

# ----------------------------
#### CHALLENGES!
# - try adding a quadratic term to get a better fit. Keep an eye on the
#   diagnostics when you do - we will pick this up in the next script
# - try fiddling with the priors. What happens if you give sigma a normal prior going from -Inf to Inf?
# - try making a linear model in terms of both Girth and Height




