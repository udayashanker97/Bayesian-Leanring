# Part A
library(mvtnorm)

data <- read.csv("Disease.csv")

standardize <- function(x) (x - mean(x)) / sd(x)
data$age <- standardize(data$age)
data$duration_of_symptoms <- standardize(data$duration_of_symptoms)
data$white_blood <- standardize(data$white_blood)

X <- model.matrix(~ gender + age + duration_of_symptoms + dyspnoea + white_blood, data = data)
y <- data$class_of_diagnosis

neg_log_post <- function(beta, X, y, tau=2) {
  eta <- X %*% beta
  loglik <- sum(y * eta - log(1 + exp(eta)))
  logprior <- -sum(beta^2) / (2 * tau^2)
  return(loglik + logprior)
}

#Optimization
beta_init <- rep(0, ncol(X))
optim_fit <- optim(beta_init, neg_log_post, X = X, y = y,
             control=list(fnscale=-1), method = "BFGS", hessian = TRUE)


beta_mode <- optim_fit$par
Hessian <- optim_fit$hessian
cov_matrix <- solve(Hessian)

names(beta_mode) <- colnames(X)
colnames(cov_matrix) <- rownames(cov_matrix) <- colnames(X)

beta_se <- sqrt(diag(-solve(Hessian)))
names(beta_se) <- colnames(X)

print(beta_mode)
print(cov_matrix)

# 95% CI for age
beta_age <- beta_mode["age"]
se_age <- beta_se["age"]
ci_age <- beta_age + c(-1.96, 1.96) * se_age

cat("Posterior mode for beta (Age):", beta_age, "\n")
cat("95% CI:", ci_age, "\n")

glmModel <- glm(class_of_diagnosis ~ 0 + ., data = data, family = binomial)
glmModel$coefficients
confint(glmModel)

# Part B
posterior_predictive <- function(new_patient, beta_mode, cov_matrix, train_data) {
  mean_age <- mean(train_data$age)
  sd_age <- sd(train_data$age)
  
  mean_duration <- mean(train_data$duration_of_symptoms)
  sd_duration <- sd(train_data$duration_of_symptoms)
  
  mean_wbc <- mean(train_data$white_blood)
  sd_wbc <- sd(train_data$white_blood)
  
  age_std <- (new_patient$age - mean_age) / sd_age
  duration_std <- (new_patient$duration - mean_duration) / sd_duration
  wbc_std <- (new_patient$white_blood - mean_wbc) / sd_wbc
  
  x_new <- c(1, new_patient$gender, age_std, duration_std, new_patient$dyspnoea, wbc_std)
  
  n_sim <- 10000
  beta_samples <- rmvnorm(n_sim, mean = beta_mode, sigma = cov_matrix)
 
  eta_new <- beta_samples %*% x_new
  pred_probs <- 1 / (1 + exp(-eta_new))
  hist(pred_probs, breaks = 50, col = "red",
       main = "Posterior Predictive for New Patient",
       xlab = "Pr(Disease = 1)")
  return(pred_probs)
}

new_patient <- list(
  age = 38,
  gender = 1, #female
  duration = 10,
  dyspnoea = 0,
  white_blood = 11000
)

# Call the function
data2 <- read.csv("Lab_2/Disease.csv")
pred_probs <- posterior_predictive(new_patient, beta_mode, cov_matrix, data2)
# hist(pred_probs, breaks = 40, col = "skyblue",
#      main = "Posterior Predictive for New Patient",
#      xlab = "Pr(Disease = 1)")

cat("Posterior mean:", mean(pred_probs), "\n")
cat("95% CI:", quantile(pred_probs, c(0.025, 0.975)), "\n")
