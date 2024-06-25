# gaussian_process_utilities.R
#
# Author: Bob Verity
# Date: 2024-06-25
#
# Inputs: (none)
#
# Outputs: (none)
#
# Purpose:
# A series of functions that are useful in Gaussian process modelling.
#
# ------------------------------------------------------------------

#install.packages("mvtnorm")
library(mvtnorm)

# define the Gaussian kernel function
Gaussian_kernel <- function(x1, x2, length_scale, sigma_f) {
  sqdist <- outer(x1, x2, "-")^2
  sigma_f^2 * exp(-0.5*sqdist / length_scale^2)
}

# force a matrix to be positive semi-definite by manually adjusting the eigenvalues
make_positive_semidefinite <- function(mat, jitter = 1e-10) {
  eig <- eigen(mat)
  eig$values[eig$values < 0] <- jitter
  eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
}

# given observations at [x,y], a series of new prediction points x_pred, and the
# parameters of the GP, produce a single draw from the posterior distribution of
# the mean
draw_posterior_mean <- function(x, y, x_pred, length_scale, sigma_f, sigma_n, mu) {
  
  n <- length(x)
  n_pred <- length(x_pred)
  
  # compute the covariance matrices
  K <- Gaussian_kernel(x, x, length_scale, sigma_f) + sigma_n^2 * diag(n)
  K_s <- Gaussian_kernel(x, x_pred, length_scale, sigma_f)
  K_ss <- Gaussian_kernel(x_pred, x_pred, length_scale, sigma_f)
  
  # compute the Cholesky decomposition
  L <- chol(K)
  
  # solve for the mean of the posterior
  alpha <- backsolve(L, forwardsolve(t(L), y - mu))
  post_mean <- mu + t(K_s) %*% alpha
  
  # solve for the covariance of the posterior
  v <- solve(L, K_s)
  post_cov <- K_ss - t(v) %*% v
  
  # force posterior covariance matrix to be positive semi-definite
  post_cov_fixed <- make_positive_semidefinite(post_cov, jitter = 0)
  
  # draw from the posterior
  rmvnorm(1, mean = post_mean, sigma = post_cov_fixed)[1,]
}
