#-------------------------------------------------------------------------------
#partA

f <- 35
n <- 78
alpha_zero <- 7
beta_zero <- 7

s <- n - f

posterior_alpha <- alpha_zero + s 
posterior_beta <- beta_zero + f  

true_mean <- posterior_alpha / (posterior_alpha + posterior_beta)
true_var <- (posterior_alpha * posterior_beta) /((posterior_alpha + posterior_beta)^2 
                                                 * (posterior_alpha + posterior_beta + 1))
true_sd <- sqrt(true_var)

print(true_sd)
set.seed(123)
nDraw <- 10000
samples <- rbeta(nDraw, posterior_alpha, posterior_beta)

means <- cumsum(samples) / (1:nDraw)
vars <- cumsum((samples - means)^2) / (1:nDraw)
sds <- sqrt(vars)

# Plot convergence of mean
plot(1:nDraw, means, type = "l", 
     xlab = "Draws", ylab = "Mean",
     main = "Convergence of Posterior Mean")
abline(h = true_mean, col = "red", lty = 2)
legend("topright", legend = c("Sample mean", "True mean"), 
       col = c("black", "red"), lty = 1:2)

# Plot convergence of standard deviation
plot(1:nDraw, sds, type = "l", 
     xlab = "draws", ylab = "SD",
     main = "Convergence of Posterior Standard Deviation")
abline(h = true_sd, col = "red", lty = 2)
legend("topright", legend = c("Sample Stand. Deviation", "True Stand. Deviation"), 
       col = c("black", "red"), lty = 1:2)


##------------------------------------------------------------------------------
#part B

# Empirical probability (from samples)
actual_probability <- mean(samples > 0.5)

# Exact probability (using pbeta)
exact_probability <- 1 - pbeta(0.5, posterior_alpha, posterior_beta)

# Print results
cat("Actual probability ", est_probability, "\n")
cat("Exact probability ", pbeta_probability, "\n")

##------------------------------------------------------------------------------
#part C

# Calculate odds from theta draws
phi_values <- samples / (1 - samples)

hist(phi_values, breaks = 50, probability = TRUE, 
     main = "Posterior of Odds (φ = θ / (1 - θ))", xlab = "φ")
lines(density(phi_values), col = "blue", lwd = 2)