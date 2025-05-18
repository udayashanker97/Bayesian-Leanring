library(mvtnorm)
library(dplyr)

# Part A
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

posterior_function <- function(beta, x, y) {
  linear_predictor <- x %*% beta
  log_likelihood <- sum(y * linear_predictor - log(1 + exp(linear_predictor)))
  log_prior <- sum(dnorm(beta, mean = 0, sd = 2, log = TRUE))
  return(log_likelihood + log_prior)
}

optim_res <- optim(par = rep(0, ncol(x)),
                   fn = posterior_function,
                   x = x,
                   y = y,
                   control = list(fnscale = -1),
                   method = "BFGS",
                   hessian = TRUE)

posterior_mode <- optim_res$par
negative_hessian <- -optim_res$hessian
covariance_matrix <- solve(negative_hessian)

coefficient_age <- posterior_mode[2]
standard_error_age <- sqrt(covariance_matrix[2, 2])
ci_age <- coefficient_age + c(-1.96, 1.96) * standard_error_age

cat("\n95% Posterior interval for coefficient of age:", ci_age)

disease_data <- disease_data %>%
  mutate(age_new = (age - mean(age))/sd(age),
         duration_new = (duration_of_symptoms - mean(duration_of_symptoms))/sd(duration_of_symptoms),
         white_blood_new = (white_blood - mean(white_blood))/sd(white_blood))

glm_model <- glm(class_of_diagnosis ~ age_new + gender + duration_new + dyspnoea + white_blood_new,
                 data = disease_data, family = binomial)

print("MLE from glm:")
print(coef(glm_model))
print("Posterior mode Values:")
print(posterior_mode)
#The values are nearly identical, confirming the Bayesian approximation is reasonable


## Part B

# ... [Previous code remains the same until Part B] ...

## Part B

posterior_predictive_distribution <- function(x_values) {
  beta_values <- rmvnorm(10000, mean = posterior_mode, sigma = covariance_matrix)
  probabilities <- as.vector(plogis(x_values %*% t(beta_values)))
  return(probabilities)
}

# Calculate standardized values using the same standardization as before
new_age <- (38 - mean(disease_data$age)) / sd(disease_data$age)
new_duration <- (10 - mean(disease_data$duration_of_symptoms)) / sd(disease_data$duration_of_symptoms)
new_white <- (11000 - mean(disease_data$white_blood)) / sd(disease_data$white_blood)

# Create the predictor vector (matches the order in your design matrix)
x_values <- c(1, new_age, 1, new_duration, 0, new_white)

posterior_predictive_values <- posterior_predictive_distribution(x_values)

# Plotting
hist(posterior_predictive_values, breaks = 30, col = "lightblue",
     main = "Posterior Predictive Distribution",
     xlab = "Probability of Disease",
     xlim = c(0, 1))
abline(v = mean(posterior_predictive_values), col = "red", lwd = 2)
legend("topright", legend = paste("Mean:", round(mean(posterior_predictive_values), 3)),
       col = "red", lwd = 2, bty = "n")

# Summary statistics
cat("\nSummary of posterior predictive distribution:\n")
print(summary(posterior_predictive_values))

cat("\nMean:", mean(posterior_predictive_values), "\n")

cat("\n95% interval:\n")
print(quantile(posterior_predictive_values, c(0.025, 0.975)))