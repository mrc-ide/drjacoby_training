
<!-- README.md is generated from README.Rmd. Please edit README.Rmd and then knit -->

# Training workshop on *drjacoby* for Bayesian model fitting in R

This repository holds the material for a training workshop on *drjacoby*.

## What to do before the workshop

#### 1. Install *drjacoby*

Please install *drjacoby* in advance so we can hit the ground running. This will also help identify any issues with installation so we can iron them out in a less pressured environment!

You should be able to install directly from the R-universe using the following code. If it asks you to install from sources, type "Y" for Yes

``` r
install.packages("drjacoby", repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
```

The R-universe method is quite new, and so may not work. If it fails then please follow [these instructions](https://mrc-ide.github.io/drjacoby/articles/installation.html) instead.

This repository assumes you are using *drjacoby* **version 1.5.4**. If you use older or newer versions, there is no guarantee that all functions will work as described here.

#### 2. Install the other packages we will use

``` r
install.packages(c("tidyverse", "odin", "mvtnorm"))
```

#### 3. Check that everything works

Please run the following short test. If it runs through without errors and finishes by printing some numbers and drawing a plot then you are good to go. If not, please get in touch *before* the workshop so we don't lose time on the day.

``` r
library(drjacoby)
library(tidyverse)
library(odin)
library(mvtnorm)

packageVersion("drjacoby")  # should be 1.5.4

mcmc <- run_mcmc(data = list(x = rnorm(10)),
                 df_params = define_params(name = "mu", min = -10, max = 10),
                 loglike = function(params, data, misc) sum(dnorm(data$x, mean = params["mu"], log = TRUE)),
                 logprior = function(params, misc) dnorm(params["mu"], sd = 10, log = TRUE),
                 burnin = 1e2, samples = 1e3, chains = 2)
mcmc$diagnostics$ess
plot_trace(mcmc)

# finally, check that odin can compile a model
gen <- odin({
  update(x) <- x + 1
  initial(x) <- 0
})
gen$new()$run(step = 0:3)
```

The final step requires a working C compiler. If everything else runs but this step fails then don't worry, there is a workaround that we can use on the day.

#### 4. Download this repository

You will need a local copy of all the material in this repository. There are a couple of ways to get it. First, you could [clone the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository). Another option would be to click on the *Code* button above and download the whole lot as a zip.

## What the workshop will cover

We will cover the following topics:

- Getting *drjacoby* up and running with the correct input formats
- Producing basic plots of MCMC output and checking MCMC diagnostics
- Checking the fit of the model against the data
- Activating parallel tempering to improve mixing and deal with multi-modality
- Linking with other packages. Example of a compartmental model in *Odin*
- A brief look at running in parallel and using C++ functions to speed things up

## Structure of the material

All the lectures are contained in the Lectures folder, and each lecture has a corresponding numbered script in the R_scripts folder:

| Lecture                | Script                                |
|:-----------------------|:--------------------------------------|
| `0.Introduction`       | (none)                                |
| `1.Basic_example`      | `01.Basic_linear_regression.R`        |
| `2.Parallel_tempering` | `02.Parallel_tempering.R`             |
|                        | `02b.Multimodal_example.R` (optional) |
| `3.Linking_with_odin`  | `03.Odin_SIR.R`                       |
| `4.Speeding_up`        | `04.Cpp_and_parallel.R`               |

The main activity of the workshop will be to work through the scripts in the R_scripts folder. These are heavily commented with instructions and challenges directly in the script. You are encouraged (actually required) to edit these scripts, so please get used to modifying them! You should work from the *drjacoby_training.Rproj* to ensure all file paths work correctly.

There are also a few bits and pieces that we will bring in from the Data and Utilities folders.

Two of these scripts are not part of the main hands-on activity. 02b.Multimodal_example.R is an optional example that will be shown as a demonstration, which you can work through if you have time. The final script, on C++ and running in parallel, will also be demonstrated in the session and is designed to be worked through in your own time.

The Extras folder contains self-study scripts that we will not cover in the session. These go further into model fit, covering credible intervals, posterior predictive intervals and posterior predictive checks.
