##ASSIGNMENT 1

library(mvtnorm)
library(BayesLogit)

disease_data <- read.csv("Disease.csv")
x <- cbind(
  Intercept = 1,
  age_new = (disease_data$age - mean(disease_data$age))/sd(disease_data$age),
  gender = disease_data$gender,
  duration_new = (disease_data$duration_of_symptoms - mean(disease_data$duration_of_symptoms))/sd(disease_data$duration_of_symptoms),
  dyspnoea = disease_data$dyspnoea,
  white_Blood_new = (disease_data$white_blood - mean(disease_data$white_blood))/sd(disease_data$white_blood)
)

y <- disease_data$class_of_diagnosis
size <- nrow(x)
col_size <- ncol(x)

tau <- 3


gibbs_implementation <- function(x, y, n_iter=1000) {
  beta_samples <- matrix(0, n_iter, ncol(x))
  
  beta_value <- rep(0, ncol(x))
  size <- nrow(x)
  col_size <- ncol(x)
  
  prior_precision <- diag(1/tau^2, col_size)
  
  for (iter in 1:n_iter) {
    z <- x %*% beta_value
    poly_gamma <- rpg(size, 1, z)
    v <- solve(t(x) %*% diag(poly_gamma) %*% x + prior_precision)
    m <- v %*% t(x) %*% (y - 0.5)
    beta_value <- as.numeric(rmvnorm(1, m, v))
    beta_samples[iter,] <- beta_value
  }
  return(beta_samples)
}

set.seed(123)
results <- gibbs_implementation(x, y)

if_evaluation <- function(beta_values) {
  acf_value <- acf(beta_values, plot=FALSE, lag.max=50)$acf[,,1]
  if_value <- 1 + 2 * sum(acf_value[-1])
  return(if_value)
}

res_if <- apply(results, 2, if_evaluation)
cat("Inefficiency Factors for whole dataset:\n")
print(res_if)

plot(1, type = "n", 
     xlim = c(1, nrow(results)), 
     ylim = range(results),
     xlab = "Iteration", 
     ylab = "Value",
     main = "Gibbs Sampler for whole dataset")

cols <- rainbow(p)

for (i in 1:p) {
  lines(results[, i], col = cols[i], lwd = 1.5)
}
legend("topright", 
       legend = paste0("beta[", 0:(p-1), "]"), 
       col = cols, 
       lty = 1, 
       lwd = 1.5,
       cex = 0.8)  

results_for_10_observations <- gibbs_implementation(x[1:10,], y[1:10])
results_for_40_observations <- gibbs_implementation(x[1:40,], y[1:40])
results_for_80_observations <- gibbs_implementation(x[1:80,], y[1:80])

res_if_10_observations <- apply(results_for_10_observations, 2, if_evaluation)
cat("Inefficiency Factors for 10 Observations:\n")
print(res_if_10_observations)

res_if_40_observations <- apply(results_for_40_observations, 2, if_evaluation)
cat("Inefficiency Factors for 40 Observations:\n")
print(res_if_40_observations)

res_if_80_observations <- apply(results_for_80_observations, 2, if_evaluation)
cat("Inefficiency Factors for 80 Observations:\n")
print(res_if_80_observations)

par(mfrow = c(1,3))
#10 Observations
plot(1, type = "n", 
     xlim = c(1, nrow(results_for_10_observations)), 
     ylim = range(results_for_10_observations),
     xlab = "Iteration", 
     ylab = "Value",
     main = "Gibbs Sampler Trajectories - All Parameters")

cols <- rainbow(p)

for (i in 1:p) {
  lines(results_for_10_observations[, i], col = cols[i], lwd = 1.5)
}
legend("topright", 
       legend = paste0("beta[", 0:(p-1), "]"), 
       col = cols, 
       lty = 1, 
       lwd = 1.5,
       cex = 0.8) 

#40 Observations

plot(1, type = "n", 
     xlim = c(1, nrow(results_for_40_observations)), 
     ylim = range(results_for_40_observations),
     xlab = "Iteration", 
     ylab = "Value",
     main = "Gibbs Sampler Trajectories - All Parameters")

cols <- rainbow(p)

for (i in 1:p) {
  lines(results_for_40_observations[, i], col = cols[i], lwd = 1.5)
}
legend("topright", 
       legend = paste0("beta[", 0:(p-1), "]"), 
       col = cols, 
       lty = 1, 
       lwd = 1.5,
       cex = 0.8)  

#80 Observations
plot(1, type = "n", 
     xlim = c(1, nrow(results_for_80_observations)), 
     ylim = range(results_for_80_observations),
     xlab = "Iteration", 
     ylab = "Value",
     main = "Gibbs Sampler Trajectories - All Parameters")

cols <- rainbow(p)

for (i in 1:p) {
  lines(results_for_80_observations[, i], col = cols[i], lwd = 1.5)
}
legend("topright", 
       legend = paste0("beta[", 0:(p-1), "]"), 
       col = cols, 
       lty = 1, 
       lwd = 1.5,
       cex = 0.8)  

par(mfrow = c(1,1))
