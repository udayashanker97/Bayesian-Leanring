library(rstan)

#Part A
ar1_simulation <- function(mu, phi, sigma_square, T) {
  results <- numeric(T)
  results[1] <- mu
  for (value in 2:T) {
    results[value] <- mu + phi * (results[value - 1] - mu) +
      rnorm(1, mean = 0, sd = sqrt(sigma_square))
  }
  return(results)
}

set.seed(123)
mu <- 5
sigma_square <- 9
T <- 300
phi_values <- c(-0.5, 0.1, 0.5)

par(mfrow = c(3,1))
for (phi_value in phi_values) {
  result <- ar1_simulation(mu, phi_value, sigma_square, T)
  plot(result, type = 'l', main = paste("AR(1) for phi =", phi_value),
       xlab = "Time (t)", ylab = "Value")
}
par(mfrow = c(1,1))

# for phi = -0.5, it shows an alternative behaviour this is because of negative 
# autocorrelation.Tendency of the values to move in opposite direction from the 
# previous value.

# Now for phi = 0.1, when you compare nearby values, in most of the case current
# value doesn't much depend on the previous value. So therefore it is having weak
# autocorrelation.

# Now for phi = 0.5,values are showing some tendency to change smoothly as the 
# time move forward. so in a way it can be considered as a moderate positive 
# autocorrelation.

#Part B
x_results_phi4 <- ar1_simulation(mu, 0.4, sigma_square, T)
x_results_phi98 <- ar1_simulation(mu, 0.98, sigma_square, T)

stan_code <- "
  data {
    int<lower=1> T;
    vector[T] results;
  }
  
  parameters {
    real mu;
    real<lower=-1, upper=1> phi;
    real<lower=0> sigma;
  }
  
  model {
    results[1] ~ normal(mu, sigma);
    for (t in 2:T) {
      results[t] ~ normal(mu + phi * (results[t - 1] - mu), sigma);
    }
  }
"

stan_model <- stan_model(model_code = stan_code)

stan_code_data4 <- list(T = T, results = x_results_phi4)
stan_code_data98 <- list(T = T, results = x_results_phi98)

fit_phi4 <- sampling(stan_model, data = stan_code_data4, chains = 5, iter = 5000)

fit_phi98 <- sampling(stan_model, data = stan_code_data98, chains = 5, iter = 5000)

# (i)

cat("\nResults when phi is 0.4:\n")
summary(fit_phi4, pars = c("mu", "phi", "sigma"), probs = c(0.025, 0.975))

cat("\nResults when phi is 0.98:\n")
summary(fit_phi98, pars = c("mu", "phi", "sigma"), probs = c(0.025, 0.975))

# Results shows that the posterior summaries from stan were able to successfully
# estimate the true values of ar1 stimulation with some minor changes of aroung
# 0.25 difference for each paramter value.

#(ii)

#Convergence of the plots

traceplot(fit_phi4, pars = c("mu", "phi")) +
  ggtitle("Convergence Plot when phi is 0.4")


traceplot(fit_phi98, pars = c("mu", "phi"))+
  ggtitle("Convergence Plot when phi is 0.98")

posterior_phi4 <- as.matrix(fit_phi4)

plot(
  posterior_phi4[, "mu"], posterior_phi4[, "phi"],
  main = "Joint Posterior of mu and phi (phi = 0.4)",
  xlab = expression(mu),
  ylab = expression(phi),
  pch = 20, col = rgb(0, 0, 1, 0.2)
)

posterior_phi98 <- as.matrix(fit_phi98)

plot(
  posterior_phi98[, "mu"], posterior_phi98[, "phi"],
  main = "Joint Posterior of mu and phi (phi = 0.98)",
  xlab = expression(mu),
  ylab = expression(phi),
  pch = 20, col = rgb(1, 0, 0, 0.2)
)
