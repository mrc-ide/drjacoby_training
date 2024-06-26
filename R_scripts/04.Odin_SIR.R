# 04.Odin_SIR.R
#
# Author: Bob Verity
# Date: 2024-06-25
#
# Inputs:
# Utilities/discrete_deterministic_sir.R
# Data/04.prev_truth.rds
# Data/04.prev_survey.rds
# Data/04.prev_survey_enhanced.rds
#
# Outputs: (none)
#
# Purpose:
# Often we want to use drjacoby alongside other packages. This is extremely
# easy, we just load the package the ordinary way and any package that you have
# in your workspace is available to use within the likelihood and prior
# functions.
# 
# We will demonstrate this through the Odin package, which can be used to solve
# systems of differential equations. We will aim to fit the classic SIR
# compartmental model to some made-up influenza data. Our data consists of a
# prevalence survey at various time points during an outbreak. Individuals were
# chosen from the population at random and tested for flu, giving us a numerator
# (number positive) and denominator (number tested). We will use these
# prevalence data to try and fit the three parameters of the SIR model:
#
# 1) beta: the transmission rate
# 2) gamma: the recovery rate
# 3) I_ini: the initial number of infected individuals at the start of the
#           simulation
#
# ------------------------------------------------------------------

# install the Odin package. Complete installation instructions can be found
# here: https://github.com/mrc-ide/odin/tree/master?tab=readme-ov-file#installation
#install.packages("odin")
library(odin)

library(tidyverse)

# --------------------------
# LOAD DATA

# our data consist of a prevalence survey
prev_survey <- readRDS("Data/04.prev_survey.rds")
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
# try to see how its working. Note that it is *not* R code, even though it might
# look like R code. Full details of the Odin syntax can be found in the package
# help
sir_generator <- odin("Utilities/discrete_deterministic_sir.R")

# let's explore the SIR model. We can create a new instance of this model with
# some parameters defined. We'll use arbitrary values for our three parameters
# of interest - feel free to fiddle with these numbers and see how they change
# the results. Note that values in this model are defined as *proportions* of
# the entire population, meaning the maximum value any state can take is 1.
sir_model <- sir_generator$new(beta = 0.3, gamma = 0.05, I_ini = 0.01)

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
# we are going to do one new thing here, using the misc argument to the
# likelihood function...

# define parameters
df_params <- define_params(name = "beta", min = 0, max = 1,
                           name = "gamma", min = 0, max = 1,
                           name = "I_ini", min = 0, max = 1)

# define the log-likelihood function. We are going to assume that the misc
# argument is a list that contains the element "return_fit". If this is TRUE
# then the function actually returns the model fit, and not the log-likelihood.
loglike <- function(params, data, misc) {
  
  # unpack parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  I_ini <- params["I_ini"]
  
  # solve ODE to calculate proportion in state I at the sampling times
  sir_model <- sir_generator$new(beta = beta, gamma = gamma, I_ini = I_ini)
  sir_prev <- sir_model$run(step = 200) |>
    as.data.frame() |>
    filter(step %in% data$time) |>
    pull(I)
  
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
    dunif(params["I_ini"], 0, 1, log = TRUE)
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

# did the MCMC converge? Did it mix? It's up to you to fix if not!
mcmc$diagnostics

plot_trace(mcmc)

plot_pairs(mcmc)

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
# the fit. Writing the same model twice creates is a waste of effort, and also
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
prev_truth <- readRDS("Data/04.prev_truth.rds")

# let's plot everything together
trajectory_df |>
  ggplot() + theme_bw() +
  geom_line(aes(x = time, y = value, group = name), col = "dodgerblue", alpha = 0.25) +
  geom_line(aes(x = time, y = prev), col = "red", data = prev_truth) +
  geom_pointrange(aes(x = time, y = proportion, ymin = lower, ymax = upper), data = prev_survey)

# ----------------------------
#### CHALLENGES!
# - first, did you manage to get the MCMC to mix!?
# - what's going on with the model fits at early timepoints? Why is it sometimes
#   starting from such a high prevalence
# - the basic reproduction number, R0, under this model is given by beta/gamma. What
#   does the posterior distribution of R0 look like? Is it more or less certain than
#   the posteriors on beta and gamma?
# - what are some weaknesses of fitting a deterministic model to stochastic data? What
#   approaches could we use to fit a stochastic SIR model?
# - imagine a parallel universe in which we had deployed an enhanced prevalence survey
#   with larger sample sizes. How would our inference differ? Here is the enhanced
#   survey data:
readRDS("Data/04.prev_survey_enhanced.rds")

