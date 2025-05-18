library(bayestestR)
#Part A
incomes <- c(22, 33, 31, 49, 65, 78, 17, 24)
mu <- 3.65
size <- length(incomes)

income_log_values <- log(incomes)
chi_square<- sum((income_log_values - mu)^2) / size

posterior_samples <- (size * chi_square)/rchisq(10000,size) 

plot(density(posterior_samples), 
     main = "Posterior Distribution of sigma square",
     xlab = "sigma square",ylab = "Density")

#Part B
sigma_samples <- sqrt(posterior_samples)

gini_samples <- 2 * pnorm(sigma_samples / sqrt(2)) - 1

plot(density(gini_samples), 
     main = "Posterior Distribution of Gini Coefficient",
     xlab = "Gini Coefficient",ylab = "Density")

#Part C
ci <- quantile(gini_samples, probs = c(0.025, 0.975))
cat("95% Equal Tail Credible Interval:", ci, "\n")


#Part D
hdi_ci <- hdi(gini_samples, ci = 0.95)
cat("95% HPDI:", hdi_ci$CI_low, hdi_ci$CI_high, "\n")

# #Part D
# #-------------------------------------------------------------------------------
# # Install and load bayestestR if not already installed
# if (!require(bayestestR)) install.packages("bayestestR")
# library(bayestestR)
# 
# # Calculate HPDI
# hpdi <- hdi(gini_samples, ci = 0.95)
# cat("95% Highest Posterior Density Interval for G:", hpdi$CI_low, hpdi$CI_high, "\n")
# 
# # Comparison
# cat("\nComparison of intervals:\n")
# cat("Equal tail CI: (", equal_tail_ci[1], ",", equal_tail_ci[2], ")\n")
# cat("HPDI: (", hpdi$CI_low, ",", hpdi$CI_high, ")\n")
# 
# # The HPDI will typically be narrower than the equal tail CI for skewed distributions
