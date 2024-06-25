# 03.Parallel_tempering.R
#
# Author: Bob Verity
# Date: 2024-06-20
#
# Inputs: (none)
#
# Outputs: (none)
#
# Purpose:
# We use the example of a quadratic model fit to the trees dataset to
# demonstrate a poorly mixing MCMC. After seeing what this looks like, we
# activate parallel tempering (PT) to solve the problem. We go through the
# various steps involved in tuning PT to get best results.
#
# ------------------------------------------------------------------
# load packages. For this script you will need the GGally package
#install.packages("GGally")
library(GGally)
library(tidyverse)

# ------------------------------------------------------------------
#### PART1 - A POORLY MIXING MCMC
# let's have a look at what a bad MCMC looks like, so we can recognise the
# features and try to remedy them

# this time our parameters data.frame will contain a "quadratic" term, along
# with the intercept and slope
df_params <- define_params(name = "intercept", min = 0, max = 100,
                           name = "slope", min = -10, max = 10,
                           name = "quadratic", min = 0, max = 1,
                           name = "sigma", min = 0, max = Inf)


# the loglikelihood function is extended to include the quadratic term
loglike <- function(params, data, misc) {
  
  # unpack parameters
  intercept <- params["intercept"]
  slope <- params["slope"]
  quadratic <- params["quadratic"]
  sigma <- params["sigma"]
  
  # our quadratic model of volume as a function of girth
  vol_predict <- intercept + slope * data$Girth + quadratic * data$Girth^2
  
  sum(dnorm(data$Volume, mean = vol_predict, sd = sigma, log = TRUE))
}

# logprior function also includes the quadratic term
logprior <- function(params, misc) {
  dunif(params["intercept"], 0, 100, log = TRUE) + 
    dunif(params["slope"], -10, 10, log = TRUE) + 
    dunif(params["quadratic"], 0, 1, log = TRUE) + 
    dlnorm(params["sigma"], meanlog = 0, sdlog = 1, log = TRUE)
}

# run MCMC
mcmc <- run_mcmc(data = trees,
                 df_params = df_params,
                 loglike = loglike,
                 logprior = logprior,
                 burnin = 1e3,
                 samples = 1e3,
                 chains = 5)

# when we look at diagnostics, we find poor convergence (rhat) an extremely low
# ESS. This is a sign that our MCMC is not mixing well
mcmc$diagnostics

# when you browse trace plots you should see chains failing to converge
plot_par(mcmc, phase = "both")

# when we produce bivariate scatterplots we can see the reason - strong
# correlations making it difficult for the sampler to move around freely
plot_cor(mcmc, "intercept", "slope")

# we can use the ggpairs function from the GGally package to produce a similar
# plot comparing all variables. All three coefficients are strongly correlated,
# although sigma appears relatively independent as we might expect
mcmc$output |>
  filter(phase == "sampling") |>
  sample_n(1e3) |>
  ggpairs(aes(col = as.factor(chain)),
          columns = c("intercept", "slope", "quadratic", "sigma"),
          lower = list(continuous = wrap("points", size = 0.5))) + theme_bw()

# we should NOT go ahead with the MCMC output! Failing to converge or mix means
# we do not have a reliable description of the posterior distribution. Our
# parameter estimates will likely be biased, and there may even be regions of
# the posterior that are completely unexplored. We need to fix this.

# ------------------------------------------------------------------
#### PART2 - TUNING PARALLEL TEMPERING
# parallel tempering (PT) can get us out of this quandry, but it also brings
# with it some extra tuning steps and diagnostic checks.

# let's re-run the MCMC, this time using 5 temperature rungs. This will take
# around 5 times longer to run. I usually decrease the number of samples when
# running these sorts of test runs to make things faster
mcmc <- run_mcmc(data = trees,
                 df_params = df_params,
                 loglike = loglike,
                 logprior = logprior,
                 burnin = 1e3,
                 samples = 1e3,
                 chains = 5,
                 rungs = 5)

# before we do any of our normal checks, we should check that PT did what we
# expected. Let's produce a "Metropolis coupling" plot
plot_mc_acceptance(mcmc)

# the x-axis in this plot shows the different temperature rungs (we used 5).
# Notice that the points on the plot fall *between* the rungs. That is because
# these points represent the percentage of time that information was swapped
# between adjacent rungs. We need to see positive values all the way along this
# plot, and ideally somewhere around 25%. For this initial MCMC run we have zero
# acceptance between rungs 1 and 2 (the "hottest" rungs). This means whatever
# rung1 was doing, it wasn't being passed up to rung2 and so was wasted effort.
# How can we fix this? First, we can try increasing the number of rungs:
mcmc <- run_mcmc(data = trees,
                 df_params = df_params,
                 loglike = loglike,
                 logprior = logprior,
                 burnin = 1e3,
                 samples = 1e3,
                 chains = 5,
                 rungs = 10)

plot_mc_acceptance(mcmc)

# this has increased the coupling rate between almost all rungs, but the rate
# between rung1 and rung2 is still stubbornly at zero. We could further increase
# the rungs, but that is starting to cost us in terms of efficiency. A better
# option is to redistribute the existing rungs. Currently they are equally
# spaced between 0 and 1, but we want to concentrate them closer to zero. For
# example:
z <- seq(0, 1, l = 10)^5
print(z)

# we can use the beta_manual argument to input these values (we don't need to
# specify the number of rungs when we do this)
mcmc <- run_mcmc(data = trees,
                 df_params = df_params,
                 loglike = loglike,
                 logprior = logprior,
                 burnin = 1e3,
                 samples = 1e3,
                 chains = 5,
                 beta_manual = z)

plot_mc_acceptance(mcmc)

# we now have positive swap rates over the entire temperature ladder, meaning we
# can be confident that PT is working as expected

# returning to our usual diagnostics, we find that things are looking much better
mcmc$diagnostics

plot_par(mcmc)

mcmc$output |>
  filter(phase == "sampling") |>
  sample_n(1e3) |>
  ggpairs(aes(col = as.factor(chain)),
          columns = c("intercept", "slope", "quadratic", "sigma"),
          lower = list(continuous = wrap("points", size = 0.5))) + theme_bw()

# at this point I would run a "final" MCMC using more samples to get the ESS up
# a bit higher
mcmc <- run_mcmc(data = trees,
                 df_params = df_params,
                 loglike = loglike,
                 logprior = logprior,
                 burnin = 1e3,
                 samples = 1e4,
                 chains = 5,
                 beta_manual = z)

# everything should now look much healthier
plot_mc_acceptance(mcmc)

mcmc$diagnostics

plot_par(mcmc)

mcmc$output |>
  filter(phase == "sampling") |>
  sample_n(1e3) |>
  ggpairs(aes(col = as.factor(chain)),
          columns = c("intercept", "slope", "quadratic", "sigma"),
          lower = list(continuous = wrap("points", size = 0.5))) + theme_bw()

# in summary, there is some extra tuning to be done when using PT, which is a
# bit annoying (the next version of drjacoby should do this tuning for you in a
# cleverer way). But once it's going, it provides a really simple way of
# improving mixing.
