y <- c(0, 2, 5, 5, 7, 1, 4)
sum <- sum(y)
lambda_values <- seq(0, 10, length.out = 1000)

prior_values <- log(sqrt(2) / (5 * sqrt(pi))) - (lambda_values^2) / (2 * 5^2)
likelihood_values <- -n * lambda_values + sum_y * log(lambda_values)
posterior_values <- prior_values + likelihood_values

posterior_unnormalized <- exp((posterior_values-max(posterior_values)))

posterior_normalized <- posterior_unnormalized / 
                          sum(posterior_unnormalized * diff(lambda_values)[1])

plot(lambda_values, posterior_normalized, type = "l", col = "blue", lwd = 2,
     xlab = "Lambda", ylab = "Posterior Value", main = "Posterior of lambda")

posterior_mode <- lambda_values[which.max(posterior_normalized)]

cat("The posterior mode is:", posterior_mode)
