####################
### ZOO800: HW13 ###
####################

# Author: Rebekkah LaBlue
# Focus: Maximum Likelihood
# Due: December 2, 2025

library(tidyverse)
library(dplyr)
library(purrr)
library(magrittr)
library(car)
library(units)
library(stats)
library(ggplot2)
library(viridis)
library(ggfortify)
library(here)


### --- OBJECTIVE 1: ANALYTICAL --- ###

dragons_df <- read.csv("dragon_data.csv")

# yi = En,j=1* X,j * Bj
# ^B = (XTX)^-1 * X' * y

x <- dragons_df$size
y <- dragons_df$acres_on_fire

X <- cbind(1, x) # matrix

XtX <- t(X) %*% X # transpose, matrix math
Xty <- t(X) %*% y
beta_hat <- solve(XtX, Xty)

beta_hat

b0_hat <- beta_hat[1] # intercept [-1.38]
b1_hat <- beta_hat[2] # slope [1.35]



### --- OBJECTIVE 2: OLS --- ###

# A #

# Construct grid (search range includes actual analytical values, ie. +/- 5)
b0_vals <- seq(b0_hat - 5, b0_hat + 5, by = 0.1) # grid around analytical intercept
b1_vals <- seq(b1_hat - 5,  b1_hat + 5,  by = 0.1) # grid around analytical slope

SSE_grid <- matrix(NA, nrow = length(b0_vals), ncol = length(b1_vals)) # grid to store parameter pairings

# Grid search (loop over all plausible values)
for (i in seq_along(b0_vals)) {
  
  for (j in seq_along(b1_vals)) {
    
    b0 <- b0_vals[i]
    b1 <- b1_vals[j]
    
    # Sum of squared errors for  b0, b1
    y_pred <- b0 + b1 * x
    SSE <- sum((y - y_pred)^2)
    
    # store results
    SSE_grid[i, j] <- SSE
  }
}

# Extract minimized OLS values
min_index_ols <- which(SSE_grid == min(SSE_grid), arr.ind = TRUE) # return row, column indices in matrix
b0_grid_ols <- b0_vals[min_coord[1, "row"]] # return indices for combos that minimize least squares
b1_grid_ols <- b1_vals[min_coord[1, "col"]]

min_index_ols # 51, 51
b0_grid_ols # val = -1.375551
b1_grid_ols # val = 1.346694



# B # 

# Define object function to be used
Min_rss <- function(par, x, y) {
  
  b0 <- par[1] 
  b1 <- par[2]
  y_pred <- b0 + b1 * x
  sum((y- y_pred)^2)
}

initial_pars <- c(b0_grid, b1_grid) # start search at vals from grid search

ols_optim <- optim(par = initial_pars, fn = Min_rss, x = x, y = y) # run optimization

ols_optim
# Convergence


# C # 

# Run, store results from several optimizations with random starting values
set.seed(333)

initial_matrix <- matrix( # define number of random starting values to check
  c(
    runif(10, -5, 5), # random intercept guesses
    runif(10, -5, 5) # random slope guesses
  ),
  ncol = 2
)
colnames(initial_matrix) <- c("b0_init", "b1_init")

initial_matrix

# Automate optimizations
results <- data.frame() # empty df to store results

for (i in 1:nrow(initial_matrix)) { # loop over all randomized starting points

  # Apply structure of object function, parameters  
  res <- optim(
    par = initial_matrix[i, ], # starting vals
    fn = Min_rss,#object function
    x = x,
    y = y
  )
  
  # Store
  results <- rbind(results,
                   data.frame(
                     init_b0 = initial_matrix[i, 1],
                     init_b1 = initial_matrix[i, 2],
                     est_b0 = res$par[1],
                     est_b1 = res$par[2],
                     converged = res$convergence
                   ))
}

results
### All results achieved convergence, demonstrating the optimization process is robust to differences in starting values for OLS.



### --- OBJECTIVE 3: MLE --- ###

# A #

# Construct grid (search range includes actual analytical values, ie. +/- 5)
# Repeat structure from above

b0_vals # same from before
b1_vals
sigma_vals <- seq(0.1, 5, by = 0.1) # new param for ML, residual s.d.; positive values only 

NLL_grid <- array(NA, dim = c(length(b0_vals), length(b1_vals), length(sigma_vals))) # grid to store parameter sets

# Loop over all combinations
for (i in seq_along(b0_vals)) {
  for (j in seq_along(b1_vals)) {
    for (k in seq_along(sigma_vals)) {
      
      b0 <- b0_vals[i]
      b1 <- b1_vals[j]
      sigma <- sigma_vals[k]
      
      # Negative log-likelihood
      y_pred <- b0 + b1 * x
      NLL <- -sum(dnorm(x = y, mean = y_pred, sd = sigma, log = TRUE))
      
      # store
      NLL_grid[i, j, k] <- NLL
    }
  }
}

# Extract minimized NLL values
min_index_mle <- which(NLL_grid == min(NLL_grid), arr.ind = TRUE)

b0_grid_mle <- b0_vals[min_index[1]]
b1_grid_mle <- b1_vals[min_index[2]]
sigma_grid_mle <- sigma_vals[min_index[3]]

min_index_mle # 51, 51, 46
b0_grid_mle # val = -1.375551
b1_grid_mle # val = 1.1346694
sigma_grid_mle # val = 4.6


# B # 

Min_nll <- function(par, x, y) {
  
  b0 <- par[1]
  b1 <- par[2]
  sigma <- par[3]
  
  if(sigma <= 0) return(Inf)  # positive sigma
  
  y_pred <- b0 + b1 * x
  -sum(dnorm(y, mean = y_pred, sd = sigma, log = TRUE)) # negative log likelihood for normal errors
}

initial_pars_nll <- c(b0_grid_mle, b1_grid_mle, sd(y - (b0_hat + b1_hat * x))) # start with vals from grid search

mle_optim <- optim(par = initial_pars_nll, fn = Min_nll, x = x, y = y) # run optimization

mle_optim
# Convergence


# C # 
# Repeat automated optimization process from above
set.seed(333)

initial_matrix_mle <- matrix( # define number of random starting values to check
  c(
    runif(10, -5, 5), # random intercept guesses
    runif(10, -5, 5), # random slope guesses
    runif(10, 0.1, 5) # positive sigma
  ),
  ncol = 3
)
colnames(initial_matrix_mle) <- c("b0_init", "b1_init", "sigma_init")

initial_matrix_mle

results <- data.frame() # empty df to store results

for (i in 1:nrow(initial_matrix_mle)) { # loop over all randomized starting points
  
  # Apply structure of object function, parameters  
  res <- optim(
    par = initial_matrix_mle[i, ],
    fn = Min_nll,
    x = x,
    y = y
  )
  
  # Store
  results <- rbind(results,
                   data.frame(
                     init_b0 = initial_matrix_mle[i, 1],
                     init_b1 = initial_matrix_mle[i, 2],
                     init_sigma = initial_matrix_mle[i, 3],
                     est_b0 = res$par[1],
                     est_b1 = res$par[2],
                     est_sigma = res$par[3],
                     converged = res$convergence
                   ))
}

results
### All results achieved convergence, demonstrating the optimization process is robust to differences in starting values for MLE.



### --- OBJECTIVE 4: COMPARISON --- ###

compared_ests <- data.frame(
  Method = c(
    "Analytical OLS", 
    "Grid Search OLS", 
    "Optim OLS",
    "Grid Search MLE", 
    "Optim MLE"
  ),
    Intercept = c(
      b0_hat, # analytical
      b0_grid_ols, # grid search OLS
      ols_optim$par[1], # Ooptim OLS
      b0_grid_mle, # grid search MLE
      mle_optim$par[1] # optim MLE
    ),
    Slope = c(
      b1_hat,    
      b1_grid_ols, 
      ols_optim$par[2], 
      b1_grid_mle,  
      mle_optim$par[2]  
    ),
    Sigma = c(
      NA,  # no sigma for ols
      NA, 
      NA,
      sigma_grid_mle, # grid search MLE
      mle_optim$par[3] # optim MLE
    )
  )

compared_ests
# Estimates are nearly identical, illustrating that OLS and MLE are the same when errors are normally distributed.