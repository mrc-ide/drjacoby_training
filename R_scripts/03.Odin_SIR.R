# 03.Odin_SIR.R
#
# Author: Bob Verity
# Date: 2024-06-25
#
# Inputs:
# Utilities/discrete_deterministic_sir.R
# Data/03.prev_truth.rds
# Data/03.prev_survey.rds
# Data/03.prev_survey_enhanced.rds
#
# Outputs: (none)
#
# Purpose:
# Often, we want to use drjacoby alongside other packages. This is extremely
# easy, we just load the package the ordinary way and any package that you have
# in your workspace is available to use within the likelihood and prior
# functions.
# 
# We will demonstrate this through the Odin package, which can be used to solve
# systems of differential equations. We will aim to fit the classic SIR
# compartmental model to some made-up influenza data. Our data consists of a
# prevalence survey at various time points during an outbreak. Individuals were
# chosen from the population at random and tested for flu, giving us a numerator
# (number positive) and denominator (number tested). We will use this
# prevalence dataset to try and fit the three parameters of the SIR model:
#
# 1) beta: the transmission rate
# 2) gamma: the recovery rate
# 3) I_init: the initial number of infected individuals at the start of the
#           simulation
#
# The first MCMC run in this script is deliberately feeble. Part way through
# you will be asked to fix it yourself, using what you learned in the previous
# scripts.
#
# ------------------------------------------------------------------

# install the Odin package. Complete installation instructions can be found
# here: https://github.com/mrc-ide/odin/tree/master?tab=readme-ov-file#installation
#install.packages("odin")
library(odin)

library(drjacoby)
library(tidyverse)

# --------------------------
# LOAD DATA

# our data consist of a prevalence survey
prev_survey <- readRDS("Data/03.prev_survey.rds")
prev_survey

# Notice that the number tested was small initially. This is typical of an
# outbreak situation, where data quality tends to improve over time

# we can plot the 95% confidence intervals on prevalence to get an idea of the
# precision of our survey
ggplot(prev_survey) + theme_bw() +
  geom_pointrange(aes(x = time, y = proportion, ymin = lower, ymax = upper)) +
  xlim(c(0, 200)) + ylim(c(0, 0.4)) +
  xlab("Time (days)") + ylab("Prevalence (proportion)")

# it looks like infections peaked somewhere around 60 days into our window, but
# there is considerable uncertainty at all time points, and especially near the
# start

# --------------------------
# LOAD THE SIR MODEL

# now let's load the SIR model from file. Please have a look into this file and
# try to see how it's working. Note that it is *not* R code, even though it might
# look like R code. Full details of the Odin syntax can be found in the package
# help
sir_generator <- odin("Utilities/discrete_deterministic_sir.R")

# NB, the line above compiles the model into C code, which requires a working C
# compiler. If it fails for you, use the following line instead. This produces
# the same model in plain R. It will run around 6 times slower, which is fine
# for our purposes
#sir_generator <- odin("Utilities/discrete_deterministic_sir.R", target = "r")

# let's explore the SIR model. We can create a new instance of this model with
# some parameters defined. We'll use arbitrary values for our three parameters
# of interest - feel free to fiddle with these numbers and see how they change
# the results. Note that values in this model are defined as *proportions* of
# the entire population, meaning the maximum value any state can take is 1.
sir_model <- sir_generator$new(beta = 0.3, gamma = 0.05, I_init = 0.01)

# now we have to run the model, which solves the deterministic system of
# equations
sir_output <- sir_model$run(step = 0:200)
head(sir_output)

# lets look at a plot of all model states over time. You should see the
# proportion of infected individuals increasing, and then peaking and crashing.
# It is this curve that we want to fit to our observed prevalence data.
sir_output |>
  as.data.frame() |>
  pivot_longer(cols = -step, names_to = "State") |>
  ggplot() + theme_bw() +
  geom_line(aes(x = step, y = value, col = State)) +
  xlab("Time (days)") + ylab("Proportion")

# --------------------------
# RUN MCMC
# by now the steps involved in running the MCMC should be quite familiar. BUT,
# we are going to do two new things here...

# define parameters
df_params <- define_params(name = "beta", min = 0, max = 1,
                           name = "gamma", min = 0, max = 1,
                           name = "I_init", min = 0, max = 1)

# NEW THING 1: keeping the likelihood fast. The likelihood is evaluated at every
# iteration of every chain, so anything slow in here gets multiplied up many
# thousands of times. Creating a new model with sir_generator$new() is
# relatively slow, so rather than doing this inside the likelihood we reuse the
# sir_model object that we already created above, and simply update its
# parameters with set_user(). For the same reason we avoid tidyverse-style data
# wrangling in here, and index the output matrix directly. Together these make
# the likelihood around 8 times faster, with identical results.
#
# NEW THING 2: the misc argument. Inside the likelihood we will assume we have
# access to an element inside misc called misc$return_fit. This is a flag that
# we can use to switch how this function works. If FALSE then the likelihood
# function operates as expected, returning the log-likelihood. But if TRUE then
# this function returns the model fit. Although this might seem strange, and it
# might be more natural to create a seperate function to return the model fit,
# having it in one place has two advantages: 1) avoids writing the same code
# twice, 2) means the model fit and likelihood used in MCMC never diverge from
# each other as we update the model.

loglike <- function(params, data, misc) {
  
  # update the parameters of our existing model
  sir_model$set_user(beta = params[["beta"]],
                     gamma = params[["gamma"]],
                     I_init = params[["I_init"]])
  
  # solve the model over days 0 to 200. The output is a matrix with one row per
  # day, meaning day t is found in row t + 1. Extract the proportion in state I
  # at the sampling times
  sir_output <- sir_model$run(step = 0:200)
  sir_prev <- sir_output[data$time + 1, "I"]
  
  # option to return model fit rather than log-likelihood
  if (misc$return_fit) {
    return(sir_prev)
  }
  
  # binomial probability of the data at each observation time
  sum(dbinom(data$n_pos, size = data$n_tested, prob = sir_prev, log = TRUE))
}

# define the logprior function to be uniform over all model parameters
logprior <- function(params, misc) {
  dunif(params["beta"], 0, 1, log = TRUE) + 
    dunif(params["gamma"], 0, 1, log = TRUE) +
    dunif(params["I_init"], 0, 1, log = TRUE)
}

# run the MCMC. We have to set the misc list so our log-likelihood function
# knows what to do
mcmc <- run_mcmc(data = prev_survey,
                 df_params = df_params,
                 misc = list(return_fit = FALSE),
                 loglike = loglike,
                 logprior = logprior,
                 burnin = 1e2,
                 samples = 1e3,
                 chains = 5)

# did the MCMC converge? Did it mix?
mcmc$diagnostics

plot_trace(mcmc, phase = "both")

plot_pairs(mcmc)

# --------------------------
# YOUR TURN - GET THE MCMC TO MIX
# that run was deliberately feeble, and the diagnostics should reflect that. We
# should not go any further until it is fixed, and fixing it is your job! There
# are (at least) two routes:
#
# 1) brute force. Increase burnin and samples until the diagnostics look OK
# 2) parallel tempering. Add rungs, then use plot_mc_acceptance() to tune alpha,
#    exactly as in the previous script
#
# Try both, by editing the run below. As a guide, no attempt should take much
# more than a minute, and you are aiming for an rhat below 1.1 and an ESS in the
# high hundreds for every parameter.
#
# QUESTION!
# which route gives you the most effective samples (ESS) per second of run-time?
# Is that what you expected? What does it tell you about when PT is worth it?

t0 <- Sys.time()
mcmc <- run_mcmc(data = prev_survey,
                 df_params = df_params,
                 misc = list(return_fit = FALSE),
                 loglike = loglike,
                 logprior = logprior,
                 burnin = 1e2,    # <- change me
                 samples = 1e3,   # <- change me
                 chains = 5)      # <- try adding rungs and alpha here
Sys.time() - t0

mcmc$diagnostics$rhat
mcmc$diagnostics$ess

plot_pairs(mcmc)

# notice the shape of the posterior in these scatterplots. This is why a simple
# sampler finds this problem hard work

# --------------------------
# VIEW MODEL FITS

# draw some parameters from the posterior
param_draws <- sample_chains(mcmc, sample_n = 100)

# set the time period over which we want to predict
time_pred <- 0:200

# calculate the model fit for each row in param_draws by using apply. This means
# every row of param_draws effectively becomes the params vector in the loglike
# function. The data element sets the times that we want the fit to be
# calculated. Finally, we must use the misc return_fit flag so we get the model
# fit and not the loglikelihood.
trajectory_mat <- apply(param_draws, 1, loglike,
                        data = data.frame(time = time_pred),
                        misc = list(return_fit = TRUE))

# the major advantage of this approach is that *we did not have to write the
# model fit twice*, once when computing the likelihood and once when exploring
# the fit. Writing the same model twice is a waste of effort, and also
# increases the chance of bugs being introduced. Doing it all in one place means
# we are definitely summarising the same model that was used in the fitting
# stage.

# let's get the trajectories into long format to make ggplot happy
trajectory_df <- trajectory_mat |>
  as.data.frame() |>
  mutate(time = time_pred) |>
  pivot_longer(cols = -time)

# at this point, we will also load in the "true" curve of infecteds over time.
# This is the actual curve that I used to generate the data - obviously we would
# not know this in reality, but it's useful here as it allows us to check how
# well our inference did.
prev_truth <- readRDS("Data/03.prev_truth.rds")

# let's plot everything together
trajectory_df |>
  ggplot() + theme_bw() +
  geom_line(aes(x = time, y = value, group = name), col = "dodgerblue", alpha = 0.25) +
  geom_line(aes(x = time, y = prev), col = "red", data = prev_truth) +
  geom_pointrange(aes(x = time, y = proportion, ymin = lower, ymax = upper), data = prev_survey)

# ----------------------------
#### CHALLENGES!
# - what's going on with the model fits at early timepoints? Why is it sometimes
#   starting from such a high prevalence?
# - the basic reproduction number, R0, under this model is given by beta/gamma. What
#   does the posterior distribution of R0 look like? Is it more or less certain than
#   the posteriors on beta and gamma?
# - what are some weaknesses of fitting a deterministic model to stochastic data? What
#   approaches could we use to fit a stochastic SIR model?
# - imagine a parallel universe in which we had deployed an enhanced prevalence survey
#   with larger sample sizes. How would our inference differ? If you want to jump into 
#   this universe, here is the enhanced survey data:
readRDS("Data/03.prev_survey_enhanced.rds")

