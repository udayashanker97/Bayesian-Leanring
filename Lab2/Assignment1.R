library(mvtnorm)

# Part A
data <- read.csv("temp_linkoping.csv")
data$time <- as.numeric(as.Date(data$date) - as.numeric(as.Date("2023-07-01"))) / 365
time_values <- seq(min(data$time), max(data$time), length=100)

conjugate_prior <- function(mu, omega, vu, sigma, time_values, n_curves=100) {
  sigma_square <- sigma * vu / rchisq(n_curves, df=vu)
  beta <- matrix(NA, nrow=n_curves, ncol=3)
  plot(data$time, data$temp, main="Prior Regression Curves",
       xlab="Time (years)", ylab="Temperature (°C)", pch=20, col="gray")
  for (value in 1:length(sigma_square)) {
    beta[value,] <- rmvnorm(1, mean=mu, sigma=sigma_square[value] * solve(omega))
    y_values <- beta[value,1] + beta[value,2]*time_values + beta[value,3]*time_values^2
    lines(time_values, y_values, col=rgb(0,0,1,0.2))
  }
  return(list(beta=beta, sigma=sigma_square))
}

mu <- c(0, -100, 100)
omega <- 0.01 * diag(3)
vu <- 1
sigma <- 1 
prior_sim <- conjugate_prior(mu, omega, vu, sigma, time_values)

# Adjusted prior
new_mu <- c(10, 0, -10)
new_omega <- 0.001 * diag(3)
new_vu <- 5
new_sigma <- 5
prior_sim_adj <- conjugate_prior(new_mu, new_omega, new_vu, new_sigma, time_values)

# Part B
simulate_posterior <- function(y, x, mu, omega, vu, sigma) {
  new_omega_value <- t(x) %*% x + omega
  new_mu_value <- solve(new_omega_value) %*% (t(x) %*% y + omega %*% mu)
  new_vu_value <- vu + length(y)
  new_sigma_value <- (vu * sigma + sum(y^2) + 
                        t(mu) %*% omega %*% mu - 
                        t(new_mu_value) %*% new_omega_value %*% new_mu_value) / new_vu_value
  
  new_sigma_square <- as.numeric(new_sigma_value) * new_vu_value / rchisq(10000, df=new_vu_value)
  posterior_beta <- matrix(NA, nrow=10000, ncol=3)
  for (i in 1:10000) {
    posterior_beta[i, ] <- rmvnorm(1, mean=as.vector(new_mu_value),
                                   sigma=new_sigma_square[i] * solve(new_omega_value))
  }
  list(sim_beta = posterior_beta,sim_sigma = new_sigma_square)
}


x <- cbind(1, data$time, data$time^2)
y <- data$temp

posterior <- simulate_posterior(y,x,mu,omega,vu,sigma)

hist(posterior$sim_beta[,1], main="Posterior of beta_zero", xlab="Value")
hist(posterior$sim_beta[,2], main="Posterior of beta_one", xlab="Value")
hist(posterior$sim_beta[,3], main="Posterior of beta_two", xlab="Value")
hist(posterior$sim_sigma, main="Posterior of sigma square", xlab="Value")

posterior_predictions <- matrix(NA, nrow=10000, ncol=length(time_values))
for (i in 1:10000) {
  posterior_predictions[i,] <- posterior$sim_beta[i,1] + posterior$sim_beta[i,2]*time_values + 
    posterior$sim_beta[i,3]*time_values^2
  
}

prediction_median <- apply(posterior_predictions, 2, median)
prediction_lower <- apply(posterior_predictions, 2, quantile, 0.05)
prediction_upper <- apply(posterior_predictions, 2, quantile, 0.95)

plot(data$time, data$temp, main="Posterior Regression", 
     xlab="Time", ylab="Temperature", pch=20)
lines(time_values, prediction_median, col="blue", lwd=2)
lines(time_values, prediction_lower, col="red", lty=2)
lines(time_values, prediction_upper, col="red", lty=2)

# Part C
min_time <- -posterior$sim_beta[,2] / (2 * posterior$sim_beta[,3])
hist(min_time, main="Posterior of Time with Minimum Temperature",
     xlab="Time", breaks=50)

# Part D
order <- 10
mu_poly <- rep(0, order+1)
omega_poly <- diag(order+1) * 0.01

# Create design matrix
x_poly <- matrix(1, nrow=length(y), ncol=order+1)
for (j in 1:order) {
  x_poly[,j+1] <- data$time^j
}

# Posterior calculations
omega_poly_post <- t(x_poly) %*% x_poly + omega_poly
mu_poly_post <- solve(omega_poly_post) %*% (t(x_poly) %*% y + omega_poly %*% mu_poly)
vu_poly_post <- new_vu + length(y)
sigma_poly_post <- (new_vu * new_sigma + sum(y^2) + 
                      t(mu_poly) %*% omega_poly %*% mu_poly - 
                      t(mu_poly_post) %*% omega_poly_post %*% mu_poly_post) / vu_poly_post

# Simulate from posterior
n_curves <- 100
sigma_poly_sim <- as.numeric(sigma_poly_post) * vu_poly_post / rchisq(n_curves, df=vu_poly_post)
beta_poly <- matrix(NA, nrow=n_curves, ncol=order+1)
for (i in 1:n_curves) {
  beta_poly[i,] <- rmvnorm(1, mean=as.vector(mu_poly_post),
                           sigma=sigma_poly_sim[i] * solve(omega_poly_post))
}

# Plot prior curves
x_poly_new <- matrix(1, nrow=length(time_values), ncol=order+1)
for (j in 1:order) {
  x_poly_new[,j+1] <- time_values^j
}

plot(data$time, data$temp, main="10th Order Polynomial Fit",
     xlab="Time", ylab="Temperature", pch=20, col="gray")
for (i in 1:n_curves) {
  y_poly <- x_poly_new %*% beta_poly[i,]
  lines(time_values, y_poly, col=rgb(0,0,1,0.2))
}