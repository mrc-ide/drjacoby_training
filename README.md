
# Training workshop on *drjacoby* for Bayesian model fitting in R

This repository holds the material for a training workshop on *drjacoby*.

## What to do before the workshop

Please install *drjacoby* in advance so we can hit the ground running. This will also help identify any issues with installation so we can iron them out in a less pressured environment!

You should be able to install directly from the R-universe using the following code. If it asks you to install from sources type "Y" for Yes
```{r}
install.packages('drjacoby', repos = 'https://mrc-ide.r-universe.dev')
```

If the R-universe method fails then please follow [these instructions](https://mrc-ide.github.io/drjacoby/articles/installation.html) instead.

## What the workshop will cover

We will cover the following topics:

- Getting *drjacoby* up and running with the correct input formats
- Producing basic plots of MCMC output and checking MCMC diagnostics
- Model fit and model comparison. Credible intervals and posterior predictive checks
- Activating parallel tempering to improve mixing and deal with multi-modality
- Linking with other packages. Example of a compartmental model in *Odin*
- Using C++ functions to speed things up

## Structure of the material

- The main activity will be to work through the scripts in the R_scripts folder, which are heavily commented with instructions. You are encouraged (actually required) to edit these scripts, so please don't feel too precious about them!
- All the lectures are contained in the Lectures folder
- There are a few bits and pieces that we will bring in from the Data and Utilities folders
