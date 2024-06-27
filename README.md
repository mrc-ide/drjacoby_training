
# Training workshop on *drjacoby* for Bayesian model fitting in R

This repository holds the material for a training workshop on *drjacoby*.

## What to do before the workshop

#### 1. Install *drjacoby*

Please install *drjacoby* in advance so we can hit the ground running. This will also help identify any issues with installation so we can iron them out in a less pressured environment!

You should be able to install directly from the R-universe using the following code. If it asks you to install from sources, type "Y" for Yes
```{r}
install.packages('drjacoby', repos = 'https://mrc-ide.r-universe.dev')
```

The R-universe method is quite new, and so may not work. If it fails then please follow [these instructions](https://mrc-ide.github.io/drjacoby/articles/installation.html) instead.

This repository assumes you are using *drjacoby* **version 1.5.4**. If you use older or newer versions, there is no guarantee that all functions will work as described here.

#### 2. Download this repository

You will need a local copy of all the material in this repository. There are a couple of ways to get it. First, you could [clone the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository). Another option would be to click on the *Code* button above and download the whole lot as a zip.

## What the workshop will cover

We will cover the following topics:

- Getting *drjacoby* up and running with the correct input formats
- Producing basic plots of MCMC output and checking MCMC diagnostics
- Model fit and model comparison. Credible intervals and posterior predictive checks
- Activating parallel tempering to improve mixing and deal with multi-modality
- Linking with other packages. Example of a compartmental model in *Odin*
- Using C++ functions and running in parallel to speed things up

## Structure of the material

All the lectures are contained in the Lectures folder, and each lecture has a corresponding numbered script in the R_scripts folder.

The main activity of the workshop will be to work through the scripts in the R_scripts folder. These are heavily commented with instructions and challenges directly in the script. You are encouraged (actually required) to edit these scripts, so please get used to modifying them! You should work from the *drjacoby_training.Rproj* to ensure all file paths work correctly.

There are also a few bits and pieces that we will bring in from the Data and Utilities folders.


