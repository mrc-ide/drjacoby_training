# 05.Cpp_and_parallel.R
#
# Author: Bob Verity
# Date: 2024-06-26
#
# Inputs:
# Utilities/trees_quadratic.cpp
#
# Outputs: (none)
#
# Purpose:
# Demonstrates how we can speed things up by writing our log-likelihood and
# log-prior functions in C++, or by running chains in parallel. Uses the example
# of a quadratic regression of the "trees" data, as explored in a previous
# script.
#
# ------------------------------------------------------------------

#install.packages("parallel")
library(parallel)
library(tidyverse)

# -----------------------------
# R FUNCTIONS

# as in an earlier example, we want to model the "trees" dataset in order to
# predict Volume as a function of Girth
head(trees)

# by now you are very familiar with the inputs to the run_mcmc() function, so
# I'll write them here in a compact form with less comments than normal
df_params <- define_params(name = "intercept", min = 0, max = 100,
                           name = "slope", min = -10, max = 10,
                           name = "quadratic", min = 0, max = 1,
                           name = "sigma", min = 0, max = Inf)


r_loglike <- function(params, data, misc) {
  vol_predict <- params["intercept"] + params["slope"] * data$Girth + params["quadratic"] * data$Girth^2
  sum(dnorm(data$Volume, mean = vol_predict, sd = params["sigma"], log = TRUE))
}

r_logprior <- function(params, misc) {
  dunif(params["intercept"], 0, 100, log = TRUE) + 
    dunif(params["slope"], -10, 10, log = TRUE) + 
    dunif(params["quadratic"], 0, 1, log = TRUE) + 
    dlnorm(params["sigma"], meanlog = 0, sdlog = 1, log = TRUE)
}

# run the MCMC using these R functions. We need parallel tempering turned on
# here in order to help with mixing, which unfortunately slows the MCMC down
mcmc <- run_mcmc(data = trees,
                 df_params = df_params,
                 loglike = r_loglike,
                 logprior = r_logprior,
                 burnin = 1e3,
                 samples = 1e4,
                 chains = 5,
                 rungs = 10,
                 alpha = 5)

# make a note of how long it took to run this MCMC, we'll see if we can improve
# on this

# -----------------------------
# C++ FUNCTIONS

# now we'll run the same analysis but using C++ version of the log-likelihood
# and log-prior. Note that for this section to work you will need to have a C++
# compiler installed on your computer. If you installed drjacoby via the
# R-universe then this might not be the case. Have a look at the installation
# instructions here (https://mrc-ide.github.io/drjacoby/articles/installation.html)
# and make sure you have Rpp installed and working.

# we will write C++ functions in a separate file, with heading .cpp, then we
# will "source" them back into this script. By far the easiest way to create
# this file in the correct format is to use the cpp_template() function. Try
# uncommenting and running this line of code, which will create a new foo.cpp
# file in your current working directory. Then you can have a look inside this
# file to see what we're dealing with
#cpp_template("foo.cpp")

# if you look inside the Utilities folder of this repository you will find a
# completed version of this file called trees_quadratic.cpp. Have a look at this
# file. Even if you're not a C++ coder, try to work through the rough logic of
# the loglike and logprior functions compared with the R versions above. You'll
# notice that R is a much more expressive language than C++, meaning we can
# achieve the same result in fewer lines of code, but on the other hand C++
# gives deeper control over the operations we are performing which can make it
# faster.

# once you're happy with the C++ functions, you can source them into this file
# as follows
Rcpp::sourceCpp("Utilities/trees_quadratic.cpp")

# then we include them in the run_mcmc() function. Note that when running C++
# functions you should put the names of the functions in quotes, i.e. "loglike",
# "logprior". This tells the MCMC that you want to use C++ style functions
mcmc <- run_mcmc(data = trees,
                 df_params = df_params,
                 loglike = "loglike",
                 logprior = "logprior",
                 burnin = 1e3,
                 samples = 1e4,
                 chains = 5,
                 rungs = 10,
                 alpha = 5)

plot_mc_acceptance(mcmc)

mcmc$diagnostics$rhat
mcmc$diagnostics$ess

plot_pairs(mcmc)

# in this example I get around a 5-times speed up using the C++ functions
# compared to the R functions. The speed up will vary depending on many factors,
# for example if you find yourself doing lots of tight nested loops then C++ may
# be much more efficient than R. But in general it's impossible to give clear
# recommendations on when you will see significant speed ups as it depends on so
# many factors. You should also keep in mind that C++ is (arguably) longer to
# write, more liable to introduce mistakes, and harder to debug. Therefore you
# should always consider the time it takes you to write code alongside the time
# it takes to run code when deciding the best course of action. If you have slow
# R code but you only plan on running the MCMC once then I would just set it
# going and go get a coffee, rather than recoding in C++.

# -----------------------------
# PARALLEL

# every chain in our MCMC is completely independent of every other. This makes
# it an "embarassingly parallel" problem, which is another way of saying there
# is no good reason *not* to be running in parallel!

# start by finding out how many cores you have to play with
cores <- parallel::detectCores()
cores

# it's up to you how many of these cores you want to allocate to drjacoby, vs
# some that you want to leave open. For simplicity we will assume that we are
# using all cores. We start by specifying a cluster object over these cores
cl <- parallel::makeCluster(cores)

# now we run the MCMC as before, but using the cluster argument. Let's start by
# using the R versions of the log-likelihood and log-prior. Note that the
# run-time reported by drjacoby is the sum of the time it took to run each
# chain, irrespective of if they were on the same or different cores. What we
# want to know is actually how long it took in seconds to run, so we'll use the
# Sys.time() function to measure this directly
t0 <- Sys.time()
mcmc <- run_mcmc(data = trees,
                 df_params = df_params,
                 loglike = r_loglike,
                 logprior = r_logprior,
                 burnin = 1e3,
                 samples = 1e4,
                 chains = 5,
                 rungs = 10,
                 alpha = 5,
                 cluster = cl)
Sys.time() - t0

# I get about a 4-times speed up on my machine, despite the fact that I have 12
# cores. This is because we are only running 5 chains, so not making the most of
# my resources. In my case, I get almost the same evaluation time if I run 10
# chains, which is why it's worth knowing how many cores you have and thinking
# about this!

# if we want to run our C++ functions over multiple cores then we have an extra
# step. After running makeCluster(), we need to give every one of our cores
# access to the sourced C++ functions
cluster_notes <- parallel::clusterEvalQ(cl, Rcpp::sourceCpp("Utilities/trees_quadratic.cpp"))

# now we run as before, using the function names in quotes and bringing in the
# cluster object
t0 <- Sys.time()
mcmc <- run_mcmc(data = trees,
                 df_params = df_params,
                 loglike = "loglike",
                 logprior = "logprior",
                 burnin = 1e3,
                 samples = 1e5,
                 chains = 10,
                 rungs = 10,
                 alpha = 5,
                 cluster = cl)
Sys.time() - t0

# I get around a 14-times speed up compared with our original run. If I use more
# chains and more samples then I can get this up to 20 times. To put that in
# context, an MCMC run that would have taken a week now becomes an overnight
# job. In some cases this can be the difference between an analysis that is
# worth doing and one that would cause you to pull your hair out!

# it is good practice to shut down the workers once we are finished
parallel::stopCluster(cl)

# ----------------------------
#### CHALLENGES!
# - try editing the C++ code in some way. For example, you could use a lognormal
#   likelihood rather than normal, or you could change the model from quadratic to something else

