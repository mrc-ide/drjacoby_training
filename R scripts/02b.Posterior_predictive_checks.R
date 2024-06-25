# 02b.Posterior_predictive_checks.R
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
# Perform posterior predictive checks on our MCMC output for the trees data analysis.
#
# ------------------------------------------------------------------

# load packages
library(tidyverse)

# we'll start by reading in completed MCMC output from the previous sessions.
# Let's focus first on the linear model
mcmc <- readRDS("Data/02.mcmc_linear_regression.rds")

# take a bunch of parameter draws from the posterior
param_draws <- mcmc$output |>
  filter(phase == "sampling") |>
  sample_n(1e4) |>
  mutate(draw = row_number())

# we want to simulate a new "synthetic" dataset for each value of our posterior
# parameter draws. This dataset should have exactly the same dimensions as the
# real data
sim_data <- param_draws |>
  group_by(draw) |>
  reframe(Girth = trees$Girth,
          Volume_sim = rnorm(nrow(trees), mean = intercept + slope*trees$Girth, sd = sigma))

# sim_data shows the sorts of values that we expect to see in a dataset with the
# same dimensions as the one we observed. If the model describes the data well,
# then summary statistics calculated on these synthetic datasets should broadly
# agree with the same summaries on the real data

# let's choose some summary statistics. The nice thing about posterior
# predictive checks is that you can choose any summary you want, but ideally
# these should be driven by sensible arguments about the sorts of things we care
# about getting right. You can think of them as sanity checks. For example, it
# may be that the maximum volume has some special importance to us - perhaps we
# are a logging company and we price large volume trees differently. In this
# case, we would want to ensure max(Volume) was on our list of summary
# statistics. We should also try and choose statistics that explore different
# qualities of the data, and are not redundant. Here are just a few basic
# examples:
sim_summaries <- sim_data |>
  group_by(draw) |>
  summarise(mean = mean(Volume_sim),
            variance = var(Volume_sim),
            min = min(Volume_sim),
            max = max(Volume_sim))

# lets be clear about what this data.frame holds: it gives summary statistics
# calculated on a range of simulated datasets drawn from the model, each using a
# different posterior draw of our model parameters
head(sim_summaries)

# now we need to calculate the same summaries on our real data
real_summaries <- data.frame(mean = mean(trees$Volume),
                             variance = var(trees$Volume),
                             min = min(trees$Volume),
                             max = max(trees$Volume)) |>
  pivot_longer(everything())

# finally, let's compare the posterior predictive distributions against the true values
sim_summaries |>
  pivot_longer(cols = -draw) |>
  ggplot() + theme_bw() +
  geom_histogram(aes(x = value), bins = 50) +
  facet_wrap(~name, scales = "free") +
  geom_vline(aes(xintercept = value), linetype = "dashed", data = real_summaries)

# for the linear model you should find that some of the true values are in the
# tails of the posterior predictive distributions. This indicates that the model
# does not fit the data well. We kind of knew this from looking directly at the
# credible intervals, but in some cases PPCs will pick up differences that
# weren't immediately obvious from looking at the model fit. They are also
# anchored in data-driven, real world measures that we care about, making them
# excellent sanity checks.

# ----------------------------
#### CHALLENGES!
# - repeat the process above, but using the quadratic model. Note that you'll have
# to change the formula for simulating synthetic data from a linear model to a
# quadratic model
# - try coming up with your own sanity check statistics
mcmc <- readRDS("Data/02.mcmc_quadratic_regression.rds")

