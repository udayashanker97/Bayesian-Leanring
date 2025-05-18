# # Load necessary libraries
# library(BayesLogit) # For Polya-Gamma sampling
# library(coda)       # For MCMC diagnostics
# 
# # Gibbs sampler for logistic regression with Polya-Gamma augmentation
# gibbs_logit <- function(X, y, tau=3, n_iter=10000, burnin=2000) {
#   n <- nrow(X)
#   p <- ncol(X)
#   
#   # Initialize parameters
#   beta <- matrix(0, nrow=p, ncol=1)
#   omega <- rep(1, n)
#   
#   # Store samples
#   beta_samples <- matrix(NA, nrow=n_iter, ncol=p)
#   
#   # Precompute X'X for efficiency
#   XtX <- crossprod(X)
#   
#   # Gibbs sampling
#   for (iter in 1:n_iter) {
#     # Sample omega from Polya-Gamma distribution
#     psi <- X %*% beta
#     omega <- rpg(n, 1, psi)
#     
#     # Sample beta from multivariate normal
#     Omega <- diag(omega)
#     V <- solve(XtX + diag(1/tau^2, p))
#     m <- V %*% crossprod(X, y - 0.5)
#     beta <- m + t(chol(V)) %*% rnorm(p)
#     
#     # Store samples after burn-in
#     if (iter > burnin) {
#       beta_samples[iter - burnin, ] <- beta
#     }
#     
#   }
#   
#   return(beta_samples)
# }
# 
# # Function to calculate inefficiency factors
# calculate_IF <- function(samples) {
#   n_samples <- nrow(samples)
#   p <- ncol(samples)
#   IFs <- numeric(p)
#   
#   for (j in 1:p) {
#     acf_vals <- acf(samples[,j], plot=FALSE, lag.max=100)$acf
#     IFs[j] <- 1 + 2 * sum(acf_vals[-1])
#   }
#   
#   return(IFs)
# }
# 
# # Example usage:
# # Assuming X is your design matrix and y is your response vector (0/1)
# # beta_samples <- gibbs_logit(X, y, tau=3)
# # IFs <- calculate_IF(beta_samples)
# # plot(mcmc(beta_samples)) # For trace plots
# 
# 
# # Wrapper function to run Gibbs sampling on subsets
# run_subsets <- function(X, y, tau=3, subsets=c(10, 40, 50), n_iter=10000, burnin=2000) {
#   results <- list()
#   
#   for (m in subsets) {
#     if (m > nrow(X)) {
#       warning(paste("Subset size", m, "larger than dataset. Using full dataset."))
#       X_sub <- X
#       y_sub <- y
#     } else {
#       X_sub <- X[1:m, ]
#       y_sub <- y[1:m]
#     }
#     
#     # Run Gibbs sampler
#     samples <- gibbs_logit(X_sub, y_sub, tau, n_iter, burnin)
#     
#     # Calculate diagnostics
#     IFs <- calculate_IF(samples)
#     
#     # Store results
#     results[[paste0("m=", m)]] <- list(
#       samples = samples,
#       IFs = IFs,
#       trace_plot = function() plot(mcmc(samples))
#     )
#   }
#   
#   return(results)
# }
# 
# # Example usage:
# # subset_results <- run_subsets(X, y, tau=3, subsets=c(10, 40, 50))
# 
# # Imports
# library(mvtnorm)
# library(BayesLogit)
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# 
# # Variables
# X <- read.csv("Disease.csv")
# 
# # Sampling Functions
# beta_prior <- function(n, b, B){
#   rmvnorm(n, mean = b, sigma = B)
# } 
# 
# 
# # Ineffeciency factor relies on the assumption that all auto correlation values are positive
# # Geyer's inefficiency therefore terminates the summation when faced with the first
# # negative auto-correlation value
# compute_if_geyer <- function(x, lag.max = 100) {
#   acf_vals <- acf(x, lag.max = lag.max, plot = FALSE)$acf[-1] # compute auto correlations
#   m <- floor(length(acf_vals) / 2)
#   psum <- acf_vals[seq(1, 2 * m, by = 2)] + acf_vals[seq(2, 2 * m, by = 2)]
#   k <- which(psum <= 0)
#   if (length(k) > 0) {
#     m <- k[1] - 1
#   }
#   pos_rho <- acf_vals[1:(2 * m)]
#   IF <- 1 + 2 * sum(pos_rho)
#   return(max(1, IF))
# }
# 
# run_experiment <- function(X){
#   d <- ncol(X)
#   n <- nrow(X)
#   # Split X and y, add y-intercept to X, type conversion
#   y <- X$class_of_diagnosis; X$class_of_diagnosis <- NULL; X <- as.matrix(X); X <- cbind(1, X); y <- as.numeric(y)
#   tau <- 3
#   I <- diag(d)
#   kappa <- y - rep(-1/2, n)
#   b <- rep(0, d)
#   B <- tau^2*I
#   
#   # Execution context
#   beta_vec <- beta_prior(n = 1, b = b, B = B) # sampling from prior
#   
#   beta_matrix <- matrix(nrow = n, ncol = d)
#   # Gibbs sampling
#   for(i in 1:n){
#     omega_vec <- rpg(n, 1, z = X[i,] %*% beta_vec)
#     omega_matrix <- diag(omega_vec)
#     V_w <- solve(t(X) %*% omega_matrix %*% X + solve(B))
#     m_w <- V_w %*% (t(X) %*% kappa + solve(B) %*% b)
#     beta_vec <- rmvnorm(1, mean = m_w, sigma = V_w)
#     beta_matrix[i,] <- beta_vec
#   }
#   
#   beta_frame <- data.frame(beta_matrix)
#   colnames(beta_frame) <- paste0("Beta", seq_len(ncol(beta_frame)))
#   beta_frame$Iteration <- 1:nrow(beta_frame)  # Add iteration column
#   beta_long <- pivot_longer(beta_frame, 
#                             cols = starts_with("Beta"), 
#                             names_to = "Coefficient", 
#                             values_to = "Value")
#   
#   # Plot
#   p <- ggplot(beta_long, aes(x = Iteration, y = Value, color = Coefficient)) +
#     geom_line(alpha = 0.9) +
#     theme_minimal(base_size = 14) +
#     labs(
#       title = "Gibbs Sampling Trace Plots for β Coefficients",
#       x = "Iteration",
#       y = expression(beta),
#       color = "Coefficient"
#     ) +
#     theme(legend.position = "right")
#   print(p)
#   
#   print(sapply(beta_frame[1:7], compute_if_geyer))
# }
# 
# for(m in c(313, 10, 40, 80)){
#   X_current <- X[1:m, ]
#   print(m)
#   run_experiment(X_current)
# }


# Install Polya-Gamma sampler
if (!require(BayesLogit)) install.packages("BayesLogit")
library(BayesLogit)

# Prepare data (same as Lab 2)
data <- read.csv("Disease.csv")
data$Age_std <- scale(data$age)
data$Duration_std <- scale(data$duration_of_symptoms)
data$White_std <- scale(data$white_blood)

X <- cbind(1,
           data$Gender,
           data$Age_std,
           data$Duration_std,
           data$Dyspnoea,
           data$White_std)
y <- data$class_of_diagnosis
n <- nrow(X)
p <- ncol(X)

# Prior
tau2 <- 9
B0 <- diag(1 / tau2, p)

# Gibbs sampler settings
n_iter <- 5000
beta_draws <- matrix(0, nrow = n_iter, ncol = p)
beta <- rep(0, p)

for (i in 1:n_iter) {
  eta <- X %*% beta
  omega <- rpg(n, h = rep(1, n), z = eta)
  
  Omega <- diag(omega)
  V_inv <- t(X) %*% Omega %*% X + B0
  V <- solve(V_inv)
  m <- V %*% (t(X) %*% (y - 0.5))
  
  beta <- as.numeric(rmvnorm(1, mean = m, sigma = V))
  beta_draws[i, ] <- beta
}

# Plot trace plots
par(mfrow = c(3, 2))
for (j in 1:p) {
  plot(beta_draws[,j], type = "l", main = paste("Beta", j - 1))
}
par(mfrow = c(1, 1))

# Inefficiency factor
library(coda)
ineff <- apply(beta_draws, 2, function(x) {
  1 + 2 * sum(acf(x, plot = FALSE)$acf[-1])
})
print(round(ineff, 2))

# Required libraries
library(BayesLogit)
library(mvtnorm)
library(coda)

# Function to run Gibbs sampler for given m
run_gibbs_subset <- function(m, n_iter = 5000, tau2 = 9) {
  subset_data <- data[1:m, ]
  X_sub <- cbind(1,
                 subset_data$Gender,
                 scale(subset_data$Age),
                 scale(subset_data$Duration_of_symptoms),
                 subset_data$Dyspnoea,
                 scale(subset_data$White_blood))
  y_sub <- subset_data$Class_of_diagnosis
  p <- ncol(X_sub)
  
  B0 <- diag(1 / tau2, p)
  beta <- rep(0, p)
  beta_draws <- matrix(0, nrow = n_iter, ncol = p)
  
  for (i in 1:n_iter) {
    eta <- X_sub %*% beta
    omega <- rpg(m, h = rep(1, m), z = eta)
    
    Omega <- diag(omega)
    V_inv <- t(X_sub) %*% Omega %*% X_sub + B0
    V <- solve(V_inv)
    m_post <- V %*% (t(X_sub) %*% (y_sub - 0.5))
    
    beta <- as.numeric(rmvnorm(1, mean = m_post, sigma = V))
    beta_draws[i, ] <- beta
  }
  
  list(draws = beta_draws)
}

# Run for m = 10, 40, 80
set.seed(123)
results <- lapply(c(10, 40, 80), function(m) run_gibbs_subset(m))

# Trace plots
par(mfrow = c(3, 2))
for (j in 1:6) {
  plot(results[[1]]$draws[,j], type = "l", main = paste("Beta", j-1, " (m = 10)"))
}
par(mfrow = c(3, 2))
for (j in 1:6) {
  plot(results[[2]]$draws[,j], type = "l", main = paste("Beta", j-1, " (m = 40)"))
}
par(mfrow = c(3, 2))
for (j in 1:6) {
  plot(results[[3]]$draws[,j], type = "l", main = paste("Beta", j-1, " (m = 80)"))
}

# Inefficiency Factors for each
inefficiency_factors <- lapply(results, function(res) {
  apply(res$draws, 2, function(x) {
    1 + 2 * sum(acf(x, plot = FALSE)$acf[-1])
  })
})
names(inefficiency_factors) <- c("m=10", "m=40", "m=80")
print(inefficiency_factors)

