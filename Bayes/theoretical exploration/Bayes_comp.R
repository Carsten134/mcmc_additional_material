# dependencies
library(tidyverse)
library(rmcmc)


## simulation of sample code #############################
# config
n <- 10
# simulation
sample_y <- function(p) {as.integer(runif(1) <= p)}

# reimplement expit to controll for different thresholds
expit <- function(x, t =1) {1/(1+exp(-(x-t)))}

sim_sample <- function(n, t) {
  x <- seq(t-6,t+6, length.out = n)
  ps <- expit(x, t)
  y <- sapply(ps, sample_y)
  return (data.frame(x = x,
                     y = y,
                     p = ps))
}

df <- sim_sample(n, 2)

## plotting sample #########################################
theme_set(theme_minimal())
df %>% ggplot(aes(x = x,
                  y = y)) +
  geom_point() +
  geom_line(aes(y = p), size = 1, color="red") +
  ggtitle("logistic regression: Assumptions and stochastic process")

## comparison between MLE and Bayes, plotting lines ########################

# training glm
lm <- glm(y ~ x,
          data = df,
          family = "binomial")

plot(lm$model$x, lm$fitted.values, type="l")


# conversion to a matrix
X <- matrix(data = c(rep(1, n), df$x), ncol = 2)

# using prior knowledge
prior <- function(beta) {
  first <- dnorm(beta[1], 2)
  second <- dnorm(beta[2])
  return(first*second)
}

# estimation

# estimate the mode given sample x
mode <- function(x) {
  den <- density(x)
  return(x[which.max(den$y)])
}

# combine mode function with implemented bayes_logit_reg_slice
est_beta_mode <- function(n,y, X, prior) {
  lm_b <- bayes_logit_reg_slice(1000, y, X,
                                prior = prior)
  return(apply(lm_b, 2, mode))
}

beta_mode <- est_beta_mode(1000, df$y, X, prior)

ps_bayes <- expit(X%*% beta_mode)

# assessing results
df$ps_bayes <- ps_bayes
df$ps_ols <- lm$fitted.values

df %>% ggplot(aes(x = x,
                  y = y)) +
  geom_point() +
  geom_line(aes(y = p), size = 1, color="red", alpha = 0.5) +
  geom_line(aes(y = ps_bayes), size = 1, color="green") +
  geom_line(aes(y = ps_ols), size = 1, color = "blue") +
  ggtitle("logistic regression: OLS vs Bayes with right priors")


## comparison between OLS and Bayes, accuracy comparison ######################

# config
ns <- c(5, 10, 20, 30, 50, 100)

# benchmarking given a sample size
ols_bayes_acc <- function(n_train, n_test, t) {
  # sim samples
  train <- sim_sample(n_train, t)
  test <- sim_sample(n_test, t)

  # OLS estimation
  lm <- glm(y ~ x,
            data = train,
            family = binomial)

  # regressor matrix
  X <- matrix(data = c(rep(1, n_train), train$x), ncol = 2)

  # Bayes estimation
  prior <- function(beta) {
    return (dnorm(beta[1], t)*dnorm(beta[2], 1))
  }
  bayes_mode <- est_beta_mode(1000, train$y, X, prior)

  # predict labels
  labs_ols <- as.integer(predict(lm, newdata = test) >= 0.5)
  X_test <- matrix(data = c(rep(1, n_test), test$x), ncol = 2)

  labs_bayes <- predict_labels_bayes_logit(beta_mode, X_test)

  # compute accuracy
  acc_ols <- sum(labs_ols == test$y)/n_test
  acc_bayes <- sum(labs_bayes == test$y) / n_test

  return(c(acc_ols, acc_bayes))
}

# computing results

acc_ols <- numeric(length(ns))
acc_bayes <- numeric(length(ns))

for (i in 1:length(ns)) {
  data <- replicate(20, ols_bayes_acc(ns[i], 100, 1))

  means <- apply(data, 1, mean)
  acc_ols[i] <- means[1]
  acc_bayes[i] <- means[2]
}


acc <- data.frame(acc_ols = acc_ols,
                  acc_bayes = acc_bayes,
                 n = ns)
# plotting results
acc %>% ggplot(aes(x = n, y = acc_ols)) +
  geom_line(color = "red", size = 1) +
  geom_line(aes(y = acc_bayes), color = "green", size = 1)+
  xlab("Sample size") +
  ylab("test accuracy") +
  ggtitle("Comparing performance MLE and Bayes")

