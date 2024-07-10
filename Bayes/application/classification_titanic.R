################################################################################
## Classification of titanic dataset using bayes ############################### 
################################################################################
library(rmcmc)

# importing the data
# remember to set working directory
df <- read.csv("data/train_data.csv")
df_test <- read.csv("data/test_data.csv")


## data preprocessing ##########################################################
df <- df[,-2:-1]

y <- df[,1]
X <- as.matrix(df[,-1])


## application of bayesian logistic regression #################################
prior <- function(beta){prod(dnorm(beta))}

beta_samples <- bayes_logit_reg_slice(1000, y, X,
                                      prior = prior,
                                      intveral_width = rep(10, length(X[1,])))


distribution_report <- function(var, burn) {
  sample <- beta_samples[burn:1000,var]
  plot(sample, type="l")
  hist(sample)
  acf(sample)
}

distribution_report(3, 300)

beta_mean <- apply(beta_samples, 2, mean)


compute_accuracy <- function(beta, X, y) {
  # predicting labels
  eta <- X %*% beta
  y_hat <- expit(eta)
  y_hat_pred <- as.numeric(y_hat >= 0.5)
  
  # computing accuracy
  accuracy = sum(y_hat_pred == y) / length(y_hat_pred)
  return(accuracy)
}

compute_accuracy(beta_mean, X, y)

## Assessing test accuracy #####################################################
y_test <- df_test$Survived

X_test <- as.matrix(df_test[,-3:-1])

compute_accuracy(beta_mean, X_test, y_test)

