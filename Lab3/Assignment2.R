library(mvtnorm)
library(MASS)
data<-read.table("eBayNumberOfBidderData_2025.dat",header=TRUE)
X <- model.matrix(~ PowerSeller + VerifyID + Sealed + 
                    Minblem + MajBlem + LargNeg + LogBook + 
                    MinBidShare, data = data)
y <- data$nBids
################################### Part A ############################
poisson_glm<-glm(y ~ . -1,data=as.data.frame(X),family=poisson)
cat("Summary\n")
print(summary(poisson_glm))
cat("\nAs observed from summary, Intercept along with VerifyID, Sealed, \n 
MajBlem, LogBook and MinBidShare are significant covariates.")
################################### Part B ############################
# Log-posterior function (with Zellner's g-prior)
log_posterior <- function(beta, X, y) {
  lambda <- exp(X %*% beta)       # exp(Xβ) = predicted Poisson rates
  log_lik <- sum(dpois(y,lambda,log=TRUE))    # Sum of log-likelihoods
  log_prior <- dmvnorm(beta,mean=rep(0,dim(X)[2]), 
                       sigma=100*solve(t(X)%*%X),log=TRUE)
  return(log_lik + log_prior)
}
# Find posterior mode
optim_result<-optim(rep(0,dim(X)[2]),log_posterior,X=X,y=y, 
                      control=list(fnscale=-1),hessian=TRUE)
beta_tilde<-optim_result$par
J_inv<-solve(-optim_result$hessian)  # Posterior covariance
names(beta_tilde)<-colnames(X)
colnames(J_inv)<-rownames(J_inv)<-colnames(X)
cat("\nBeta_tilde:\n")
print(beta_tilde)
cat("\nCovariance Matrix (Jy(-1) beta_tilde)\n")
print(J_inv)
################################### Part C ############################
RWMSampler <- function(logPostFunc,theta_init,cov_prop,n_iter,burn_in,c,...) {
  # logPostFunc: Function to compute the log-posterior density. First argument 
  # must be `theta`.
  # theta_init: Initial parameter vector.
  # cov_prop Proposal covariance matrix.
  # n_iter Total number of iterations.
  # burn_in Burn-in samples.
  # c Step size scaling factor.  # Tune for 25-30% acceptance
  p<-length(theta_init)
  theta<-matrix(NA, n_iter, p)
  theta[1, ]<-theta_init
  log_post_current<-logPostFunc(theta[1, ], ...)
  n_accept<-0
  for (i in 2:n_iter) {
    # Propose new theta
    theta_prop<-MASS::mvrnorm(1,mu=theta[i-1, ],Sigma=c*cov_prop)
    # Compute log-posterior at proposal
    # Log acceptance probability (avoid numerical overflow)
    log_alpha<-logPostFunc(theta_prop,...)-logPostFunc(theta[i-1, ], ...)
    if(log(runif(1)) < log_alpha){
      theta[i, ]<-theta_prop
      n_accept<-n_accept+1
    } else{
      theta[i, ]<-theta[i-1, ]
    }
  }
  theta<-theta[(burn_in+1):n_iter, ]
  acceptance_rate<-n_accept/n_iter
  cat("Acceptance rate:", acceptance_rate,"\n")
  return(theta)
}

LogPostPoisson <- function(theta, X, y) {
  lambda<-exp(X%*%theta)
  log_lik<-sum(dpois(y,lambda,log=TRUE))
  log_prior<-mvtnorm::dmvnorm(theta,mean=rep(0,ncol(X)),
                              sigma=100*solve(t(X)%*%X),log=TRUE)
  return(log_lik+log_prior)
}
samplesMH<-RWMSampler(logPostFunc=LogPostPoisson,theta_init=beta_tilde,
                     cov_prop=J_inv,n_iter=10000,burn_in=500,c=0.7,X=X,y=y)
colnames(samplesMH)<-colnames(X)
plot.ts(samplesMH, main = "Trace Plot for MH estimated Beta(s)")
cat("\nSummary Results\n")
print(data.frame(Covariate=colnames(samplesMH),
           Mean=round(colMeans(samplesMH),4),
           Std=round(apply(samplesMH,2,sd),4)),row.names=FALSE)
cat("\nComparing the values with those in part B, we find while values for
    Intercept, Sealed , LogBook and MinBidShare are similar, other covariates
    show significant differences in values\n")

################################### Part D ############################
x_new<-c(1,1,0,1,0,1,0,1.3,0.7)
lambda_new<-exp(samplesMH %*% x_new)
y_new<-rpois(nrow(samplesMH),lambda_new)
hist(y_new,breaks=30,main="Predictive Distribution for Number of Bidders")
cat("\nProbability of no bidders (y_new = 0):",mean(y_new==0),"\n")
