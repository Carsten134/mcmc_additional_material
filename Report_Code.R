## Monte Carlo Simulations: Theory and Practise
## Markov Chain Monte Carlo

## Working directory:

setwd("~/Studium/Master Statistik/Monte Carlo Simulations/Programming")

## Packages:

library("mvtnorm")
library("MASS")
library("microbenchmark")
library("condMVNorm")

################################################################################
## Slice Sampling:
################################################################################

## Simplest Slicing Algorithm for one dimension:

## First try using inverse density for determining the slice:
## it: Number of iterations
## x0: Initial value
## f: f(x), density function (or proportional function)
## inv.l, inv.u functions for lower and upper bounds of slice depending on u.
simpleslice <- function(it, x0, f, inv.l, inv.u) {
  x <- numeric(it)   # vector of samples
  x[1] <- x0   
  for(i in 2:it) {
    ## Sampling of auxiliary u
    u <- runif(1, min = 0, max = f(x[i - 1]))
    
    ## Sampling from slice
    x[i] <- runif(1, min = inv.l(u), max = inv.u(u))
  }
  
  return(x)
}


## Try it out:
## Standard normal distribution with bad initial value:
set.seed(74765)
x.firsttest <- simpleslice(10000, x0 = 15, f = function(x) {exp(-x^2/2)/sqrt(2 * pi)},
                     inv.l = function(x) {- sqrt(-2 * log(x * sqrt(2 * pi)))},
                     inv.u = function(x) {sqrt(-2 * log(x * sqrt(2 * pi)))})

hist(x.firsttest[-1], breaks = 100, main = "First try sampler", xlab = "x")
which.max(x.firsttest[-1])    
x.firsttest[1:10]
## really fast burn-in. only first 3 samples are unlikely for target 
## distribution.
## Better choice of x0 results in faster burn-in


################################################################################
## One-dimensional Slice Sampler with stepping-out and shrinkage procedure
################################################################################

## Input:
## it: number of iterations
## x0: starting value
## w: tuning parameter: Initial slice size
## f: function of x, density of desired distribution or proportional to that
##    distribution.
## L.lim, R.lim: lower and upper limit of the slice. [-Inf, Inf] is default,
##               but has to be specified if distribution is not defined on for
##               all real numbers (e.g. exp-distribution needs L.lim = 0)

OneDSlice <- function(it, x0, w, f, L.lim = -Inf, R.lim = Inf) {

  ## Vector of simulated data:
  x <- numeric(it)
  x[1] <- x0
  
  for(i in 2:it) {
    ## uniformly generate auxiliary variable
    u <- runif(1, min = 0, max = f(x[i - 1])) 
    
    ## Determine sampling interval:
    L <- x[i - 1] - w * runif(1)
    ## Set L to L.lim if Limit is crossed
    if(L < L.lim) {
      L <- L.lim
    }
    
    ## Same for R.lim
    R <- L + w
    if(R > R.lim) {
      R <- R.lim
    }
    
    ## Expand interval until lower and upper bound are outside slice:
    while(f(L) > u) {
      L <- L - w
      if(L < L.lim) {
        L <- L.lim
        break
      }
    }
    while(f(R) > u) {
      R <- R + w
      if(R > R.lim) {
        R <- R.lim
        break
      }
    }
    
    ## sample from [L, R]. Accept, if xstar is inside slice, else shrink interval
    ## to [xstar, R] or [L, xstar] and repeat
    repeat {
      xstar <- runif(1, L, R)
      if(f(xstar) > u) {
        x[i] <- xstar
        break
      }
      else {
        if(xstar < x[i - 1]) {L <- xstar}
        if(xstar > x[i - 1]) {R <- xstar}
      }
    }
  
  }
  return(x)
}


################################################################################
## Running some examples:

## Standard Normal distribution:
## We choose w = sd/8, so it scales to spread of desired distribution
set.seed(645287)
xnormal <- OneDSlice(10000, x0 = 1, w = 1/8,  
                   f = function(x) {exp(-x^2/2)/sqrt(2 * pi)})

hist(xnormal, breaks = 100, main = "Sampled N(0, 1)-distribution", 
     xlab = "x")
## distribution looks correct

par(mfrow = c(2, 1), mar = c(4.5, 4, 1, 1))
plot(xnormal[1:500], type = "l", ylab = "x", 
     xlab = "Trace plot and ACF for sampled N(0, 1)")
acf(xnormal, main = "")
## Low autocorrelation. Algorithm works well in this simple case


################################################################################
## Standard Normal distribution with bad initial value:
xnormal2 <- OneDSlice(100000, x0 = 20, w = 1,  
                     f = function(x) {exp(-x^2/2)/sqrt(2 * pi)})

par(mfrow = c(1, 1))
hist(xnormal2, breaks = 100, xlab = "x", 
     main = "Initial value x0 = 20" )

par(mfrow = c(2, 1)) 
plot(xnormal2[1:500], type = "l",  ylab = "x", 
     xlab = "Trace plot and ACF for bad initial value")

acf(xnormal2, main = "")
## Fast burn-in, might miss high-density region a couple of times, but once its
## hit, algorothm likely stays there.


################################################################################
## Multimodal case:
## Use mixture of two normal distributions.
xmm <- OneDSlice(10000, x0 = 1, w = 1/8,
                 f = function(x) {exp(-x^2/2)/sqrt(2 * pi) +
                                  exp(-(x - 4)^2/2)/sqrt(2 * pi)})

par(mfrow = c(1, 1))
hist(xmm, freq = FALSE, ylim = c(0, 0.25), breaks = 100, xlab = "x",
     main = "two normal distributions")
curve((exp(-x^2/2)/sqrt(2 * pi) +
         exp(-(x - 4)^2/2)/sqrt(2 * pi))/2, 
      from = -4, to = 8, add = TRUE)
## works as well.

## This plot is interesting, as it shows us, how switch between the two peaks
## of the distribution at 0 and 4.
## It makes the trace plot harder to interpret though. We might read it as 
## high autocorrelation, even if thats not necessarily the case.
plot(xmm[1:200], type = "l", main = "Trace plot of 2 normal distributions", 
     ylab = "x")

## autocorrelation:
acf(xmm)
## Looks decent. not too high after lag 7




################################################################################
## What if distribution is bounded, like exponential distribution?
## In this case we need to set the lower limit for the slice L.lim = 0.

xexp <- OneDSlice(10000, x0 = 1, w = 1/8, 
                  f = function(x) {0.5 * exp(-x * 0.5)},
                  L.lim = 0)

hist(xexp, freq = FALSE, breaks = 100, xlab = "x", 
     main = "exp-distribution")
curve(0.5 * exp(-x * 0.5), from = 0, to = 10, add = TRUE)

plot(xexp[1:500], type = "l",, ylab = "x", 
     main = "Trace plot of exp-dsitribution")

## autocorrelation:
acf(xexp)
## Bit higher than for the standard normal case.




################################################################################
## Multidimensional case:
################################################################################


## Sizes w_i again should be chosen with regards to spread of the distribution,
## so scaled to sd_i

## Input is the same, only this time x0, w, L.lim and R.lim are vectors 
## Ouput is matrix with rows being the individual samples and columns being the
## components.

MVSlice <- function(it, x0, w, f, L.lim, R.lim) {
  n <- length(x0)  # number of dimensions
  x <- matrix(numeric(it * n) , nrow = it) # matrix of samples
  x[1 ,] <- x0 # initial value
  
  for(i in 2:it) {
    ## uniformly generate auxiliary variable
    u <- runif(1, min = 0, max = f(x[i - 1,])) 
    
    ## Determine sampling region:
    L <- x[i - 1,] - w * runif(n)
    L[L < L.lim] <- L.lim[L < L.lim]   ## Set border to L.lim if necessary 
    
    R <- L + w
    R[R > R.lim] <- R.lim[R > R.lim]   ## Set border to R.lim if necessary 
    
    ## sample from [L, R]. Accept, if x_star is in [L, R], else shrink region
    ## to coordinates of sample
    repeat {
      xstar <- sapply(1:n, FUN = function(y) {runif(1, min = L[y], max = R[y])})
      
      if(f(xstar) > u) {
        x[i ,] <- xstar
        break
      }
      
      else {
        L[xstar < x[i - 1 ,]] <- xstar[xstar < x[i - 1 ,]]
        R[xstar > x[i - 1 ,]] <- xstar[xstar > x[i - 1 ,]]
      }
    }
  }
  
  return(x)
}


################################################################################
## Running some examples:

## Define f:
## 2D normal distributions with correlation 0.9.
fnormal <- function(x) {
  mu = c(0, 0)
  Sigma = matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  n <- length(x)
  fx <- 1/sqrt((2 * pi)^n * det(Sigma)) *
         exp(-0.5 * t(x - mu) %*% solve(Sigma) %*% (x - mu))
  return(fx)
}

mvnormaltest <- MVSlice(10000, x0 = c(0, 0), w = c(10, 10), f = fnormal, 
                        L.lim = c(-Inf, -Inf), R.lim = c(Inf, Inf))

par(mfrow = c(1, 1))
plot(mvnormaltest, ylim = c(-5, 5), xlim = c(-5, 5),
     col = "black", ylab = "x2", xlab = "x1", 
     main = "Slice Sampling")
## Looks good.

## Trace plot: 
par(mfrow = c(2, 1))
plot(1:1000, mvnormaltest[1:1000, 1], type = "l")
plot(1:1000, mvnormaltest[1:1000, 2], type = "l")
## Looks good

################################################################################
## Examine convergence with bad initial value:

set.seed(739264)
mvnormalbadinit <- MVSlice(1000, x0 = c(-8, 8), w = c(50, 50), f = fnormal, 
                           L.lim = c(-Inf, -Inf), R.lim = c(Inf, Inf))

## Plot whole data:
par(mfrow = c(1, 1))
plot(mvnormalbadinit, ylab = "x2", xlab = "x1",
     col = rep(c("red", "black"), c(1, 999)), 
     pch = rep(c(19, 1), c(1, 999)))
text(mvnormalbadinit[1 , 1], mvnormalbadinit[1 , 2],
     labels = "inital value", pos = 4)   
## Just a couple "bad" values before hitting the right distribution

## trace plots:
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
plot(1:300, mvnormalbadinit[1:300, 1], type = "l", 
     xlab = "t", ylab = "x1")
plot(1:300, mvnormalbadinit[1:300, 2], type = "l",
     xlab = "t", ylab = "x2")

## autocorrelations:
acf(mvnormalbadinit[, 1], ylab = "ACF for x1", 
    xlab = "ACF plots for w = c(50, 50)")
acf(mvnormalbadinit[, 2], ylab = "ACF for x2", xlab = "")

par(mfrow = c(1, 1))

################################################################################
## Try out different choices of w:
## Choice of w is really important, w too small leads to bad results
## Larger w results in longer runtime but gives better results.

## c(1, 1) is still kind of ok but c(0.2, 0.2) is a desaster.
mvnormalsmall <- MVSlice(1000, x0 = c(0, 0), w = c(1, 1), f = fnormal, 
                        L.lim = c(-Inf, -Inf), R.lim = c(Inf, Inf))

mvnormaltiny <- MVSlice(1000, x0 = c(0, 0), w = c(0.2, 0.2), f = fnormal, 
                         L.lim = c(-Inf, -Inf), R.lim = c(Inf, Inf))

par(mfrow = c(1, 2))
plot(mvnormalsmall, ylim = c(-2.5, 2.5), xlim = c(-2.5, 2.5),
     main = "w = (1, 1)", ylab = "x2", xlab = "x1",
     col = rep(rainbow(10), each = 100))
plot(mvnormaltiny, ylim = c(-2.5, 2.5), xlim = c(-2.5, 2.5), 
     main = "w = (0.2, 0.2)", ylab = "x2", xlab = "x1",
     col = rep(rainbow(10), each = 100))

## Plot on the left looks decent. limits to stepsize a little visible due to
## coloring of points (100 consecutive iterations have same color)
## Right side looks really bad in comparison. Step size too small, high 
## autocorrelation, and it does not reproduce the desired distribution. Combine
## that with a bad initial value and you have a recipe for desaster.


################################################################################
## Comparing runtimes:

## Compare runtimes for w = (0.2, 0.2), (1, 1), (10, 10) and (50, 50) 
time50 <- system.time(replicate(100, MVSlice(100, x0 = c(0, 0), w = c(50, 50), 
                      f = fnormal, 
                      L.lim = c(-Inf, -Inf), R.lim = c(Inf, Inf))))[3]

time10 <- system.time(replicate(100, MVSlice(100, x0 = c(0, 0), w = c(10, 10), 
                     f = fnormal,
                     L.lim = c(-Inf, -Inf), R.lim = c(Inf, Inf))))[3]

time1 <- system.time(replicate(100, MVSlice(100, x0 = c(0, 0), w = c(1, 1), 
                     f = fnormal, 
                     L.lim = c(-Inf, -Inf), R.lim = c(Inf, Inf))))[3]

time0.2 <- system.time(replicate(100, MVSlice(100, x0 = c(0, 0), 
                       w = c(0.2, 0.2), f = fnormal, 
                       L.lim = c(-Inf, -Inf), R.lim = c(Inf, Inf))))[3]


data.frame("w_i" = c(0.2, 1, 10, 50), 
           "runtime" = c(time0.2, time1, time10, time50))
## runtime for w = 50 only about 4 time larger than for w = 1. Not dramatic
## But this could be higher in higher dimensions.


## Conclusion: Better to choose a w that is a bit too big, and risk longer 
## runtime than choosing small w for shorter time and bad sampling.



################################################################################
## Application to the real world:
################################################################################

################################################################################
## The Titanic data set:

## Required packages:
library(rmcmc)
library(tidyverse)

## Import data:
titanic <- read.csv("titanic_train.csv")
titanic_test  <- read.csv("titanic_test.csv")


## Preprocessing of training data:

## NA treatment: set NAs to mean:
sapply(titanic, function(x) any(is.na(x))) ## only Age contains NA

titanic$Age[is.na(titanic$Age)] <- mean(titanic$Age, na.rm = TRUE)


## Convert categorical variables into dummies:
n <- nrow(titanic)
titanic_train <- data.frame(Pclass_1 = integer(n),
                       Pclass_2 = integer(n),
                       Pclass_3 = integer(n),
                       Embarked_Q = integer(n),
                       Embarked_S = integer(n),
                       Embarked_C = integer(n),
                       Sex = integer(n),
                       SibSp = titanic$SibSp,
                       Parch = titanic$Parch,
                       Fare = titanic$Fare,
                       Age = titanic$Age,
                       intercept = rep(1, n))

vals_emb <- c("Q", "S", "C")

for (i in 1:3) {
  titanic_train[,i] <- as.integer(titanic$Pclass == i)
  titanic_train[,(i+3)] <- as.integer(titanic$Embarked == vals_emb[i])
}

titanic_train$Sex = as.integer(titanic$Sex == "female")

## Normalize Parch, SibSp andFare so it is less skewed:

## logscaling skewed distributions
## Adding + 1 such that log(x) != -infinity
skewed_plus_one <- c("Parch", "SibSp", "Fare")

for (key in skewed_plus_one) {
  titanic_train[[key]] <- log(titanic_train[[key]] +1)
}


## Normalizing:
normalize_var <- function(df, key) {
  return((df[[key]] - mean(df[[key]], na.rm=T))/ sqrt(var(df[[key]], na.rm=T)))
}

to_normalize <- c("Age", "Fare", "Parch", "SibSp")
for (key in to_normalize) {
  titanic_train[[key]] <- normalize_var(titanic_train, key)
}


plot_distributions <- function(df, keys) {
  for (key in keys) {
    hist(titanic_train[[key]], main = paste("Histogram of: ", key))  
  }
}
plot_distributions(titanic_train, to_normalize)

## Log scaling didn't quite work for Parch and SibSp but will go with it for now

?MVSlice
?bayes_logit_reg_slice

## Define X and y of Bayesian logistic regression:
X <- as.matrix(titanic_train)
y <- titanic$Survived

beta_sample <- bayes_logit_reg_slice(10000, y, X)

head(beta_sample)

## playing around with properties of the samples
plot_props <- function(df, key) {
  plot(df[,key], type="l", main=paste("history of vals for: ", key ), 
       xlab = "t", ylab=key)
  plot(density(df[,key]), main=paste("distribution for:", key))
  acf(df[,key])
}



## Trace plots for all variables:
par(mfrow = c(4, 3), mar = c(1, 2, 1, 1))
sapply(1:11, FUN = function(key) {
  plot(beta_sample[,key], type="l", main = NA, 
       xlab = "t", ylab=key)
})

?sapply


## Assessing Accuracy:

?predict_labels_bayes_logit
accuracy <- function(beta){sum(predict_labels_bayes_logit(beta, X) == y)/
                           length(y)}

## Using mean:
beta_mean <- apply(beta_sample, 2, mean)
accuracy(beta_mean)

## Using median:
beta_median <- apply(beta_sample, 2, median)
accuracy(beta_median)

## Using mode:
mode <- function(x) {
  den <- density(x)
  return(x[which.max(den$y)])
}

beta_mode <- apply(beta_sample, 2, mode)
accuracy(beta_mode)


## Now test performance on titanic_test:

# doing random search to find best threshold:
n <- 1000
ps <- runif(n)
best_accuracy <- -1
best_p <- 0
for (p in ps) {
  acc <- sum(predict_labels_bayes_logit(beta_mean, X, p) == y) /length(y)
  
  if (acc > best_accuracy) {
    best_p <- p
    best_accuracy <- acc
  }
}

print(paste("best p: ", best_p))
print(paste("best accuracy: ", best_accuracy))

## Preprocess test data:

# like above exchange with mean
titanic_test$Age[is.na(titanic_test$Age)] <- mean(titanic_test$Age, na.rm=T)
titanic_test$Fare[is.na(titanic_test$Fare)] <- mean(titanic_test$Fare, na.rm=T)

preprocess_data <- function(df) {
  # prior information
  n <- length(df$Pclass)
  
  # instatiating result
  result <- data.frame(Pclass_1 = integer(n),
                       Pclass_2 = integer(n),
                       Pclass_3 = integer(n),
                       Embarked_Q = integer(n),
                       Embarked_S = integer(n),
                       Embarked_C = integer(n),
                       Sex = integer(n),
                       SibSp = df$SibSp,
                       Parch = df$Parch,
                       Fare = df$Fare,
                       Age = df$Age,
                       intercept = rep(1, n))
  
  vals_emb <- c("Q", "S", "C")
  
  for (i in 1:3) {
    result[,i] <- as.integer(df$Pclass == i)
    result[,(i+3)] <- as.integer(df$Embarked == vals_emb[i])
  }
  
  result$Sex = as.integer(df$Sex == "female")
  
  # log scaling
  # Adding + 1 such that log(x) != -infinity
  skewed_plus_one <- c("Parch", "SibSp", "Fare")
  
  
  for (key in skewed_plus_one) {
    result[[key]] <- log(result[[key]] +1)
  }
  
  # normalizing
  to_normalize <- c("Age", "Fare", "Parch", "SibSp")
  for (key in to_normalize) {
    result[[key]] <- normalize_var(result, key)
  }
  
  return(as.matrix(result))
}

## Process test data:
X_test <- preprocess_data(titanic_test)

## predict:
y_hat <- predict_labels_bayes_logit(beta_mode, X_test, best_p)

## Now how to assess the quality?
y_hat



################################################################################

## Try out Slicing with synthetic data of 8 variables, since that seems to be 
## the maximum dmcnorm can handle

## define f using the function dmvnorm from the package mvtnorm:
## use random covariances:

cov8 <- matrix(sample(seq(0.1, 0.9, by = 0.1), 64, replace = TRUE), 
                nrow = 8)
cov8[lower.tri(cov8)] <- t(cov8)[lower.tri(cov8)]
diag(cov8) <- rep(1, 8)

make_positive_semidefinite <- function(mat) {
  eig <- eigen(mat)
  eig$values[eig$values < 0] <- 0
  mat_psd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  return(mat_psd)
}

cov8 <- make_positive_semidefinite(cov8)
cov8


f8 <- function(x, mean = rep(0, 8), cov = cov8, regularization = 1e-6) {
  cov <- cov + diag(regularization, nrow(cov))
  return(dmvnorm(x, mean = mean, sigma = cov))
}

f8(rep(0, 8))


sapply(seq(from = -0.2, to = 0.2, by = 0.01), 
       FUN = function(x) {f8(rep(x, 8))} )





## Try out Bayesian Logistic Regression with 100 Variables:

## Use different distributions:


musigma <- cbind(sample(-2:2, 50, replace = TRUE),
                 rep(1, 50))

df <- do.call(cbind,lapply(1:50, FUN = function(x) {
  rnorm(300, mean = musigma[x, 1], sd = musigma[x, 2])
}))

beta <- rnorm(50)

lin_comb <- df %*% beta

probs <- 1 / (1 + exp(-lin_comb))

y <- rbinom(300, size = 1, prob = probs)


testsample <- bayes_logit_reg_slice(10000, y, df, 
                                    intveral_width = rep(50, 50),
                                    starting_values = rep(0, 50))

head(testsample)



par(mfrow = c(3, 2), mar = c(1, 2, 1, 1))
plot(1:10000, testsample[, 1], type = "l")
plot(1:10000, testsample[, 12], type = "l")
plot(1:10000, testsample[, 16], type = "l")
plot(1:10000, testsample[, 26], type = "l")
plot(1:10000, testsample[, 34], type = "l")
plot(1:10000, testsample[, 46], type = "l")


colMeans(df[, c(1, 12, 16, 26, 34, 46)])





## Less variables:

testsample2 <- bayes_logit_reg_slice(10000, y, df[, 1:10], 
                                    intveral_width = rep(50, 10),
                                    starting_values = rep(0, 10))

par(mfrow = c(3, 2), mar = c(1, 2, 1, 1))
plot(1:10000, testsample2[, 1], type = "l")
plot(1:10000, testsample2[, 2], type = "l")
plot(1:10000, testsample2[, 3], type = "l")
plot(1:10000, testsample2[, 4], type = "l")
plot(1:10000, testsample2[, 5], type = "l")
plot(1:10000, testsample2[, 6], type = "l")

colMeans(df[, 1:6])


## Covariance martix Utilities #################################################

## make positive semidefinite by spectral decomposition
make_positive_semidefinite <- function(mat) {
  eig <- eigen(mat)
  eig$values[eig$values < 0] <- 0
  mat_psd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  return(mat_psd)
}

## utility for different
GenSigma <- function(d, diag_val) {
  # generate entries randomly
  sigma <- matrix(sample(seq(0.1, 0.9, by = 0.1), d^2, replace = TRUE),
         nrow = d)
  # ensure symmetry
  sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]

  # set diagonal elements
  diag(sigma) <- rep(diag_val, d)

  sigma <- make_positive_semidefinite(sigma)
  return(sigma)
}

## Gibbs Implementation ########################################################
MVNormGibbs <- function(n, mu, simga, init) {
  d <- length(mu)

  sample_raw <- numeric(d*n)
  sample_raw[1:d] <- init

  for(i in 2:n) { # for each observation...
    for (j in 1:d) { # for each dimension...
      current <- d*(i-1)+j
      x.given <- sample_raw[(current-(d-1)):(current-1)]
      # condition the density for the dimension on the last k-1 values
      x <- condMVNorm::rcmvnorm(1, mu, sigma,
                                dependent.ind = j,
                                given.ind = (1:d)[-j],
                                X.given = x.given)

      # ...and then append the sample
      sample_raw[current] <- x
    }
  }
  sample <- list()
  for(j in 1:d) {
    sample[[as.character(j)]] <- sample_raw[seq(j, d*n, d)]
  }

  return(data.frame(sample))
}

## Simulation ##################################################################
sigma10 <- GenSigma(10, 6)
mu10 <- rep(0, 10)
sigma100 <- GenSigma(100, 6)
mu100 <- rep(0, 100)


Gibbs10 <- MVNormGibbs(1000, mu10, sigma10, rep(0, 10))
Gibbs100 <- MVNormGibbs(1000, mu100, sigma100, rep(0, 100))

## Benchmarking against the default normal #####################################
bnormal100 <- rmvnorm(1000, mu100, sigma100)
bnormal10 <- rmvnorm(1000, mu10, sigma10)

## plotting the results ########################################################

plot_comparison <- function(sim10, sim100, true10, true100) {
  par(mfrow = c(3, 2), mar = c(4, 4.2, 1.5, 1),
      cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  # histograms against true density
  hist_data <- hist(sim10[,1], plot = FALSE)
  hist(sim10[,1], xlab = "x", main = substitute(paste(bold("true 100 dimensional normal distribution"))), freq = FALSE)
  x_vals <- seq(min(sim10[,1]), max(sim10[,1]), length = 100)
  y_vals <- dnorm(x_vals, mean = mean(true10[,1]), sd = sd(true10[,1]))
  lines(x_vals, y_vals, col = "red", lwd = 2)

  hist_data <- hist(sim100[,1], plot = FALSE)
  hist(sim100[,1], xlab = "x",ylim = c(0, max(y_vals)), main = substitute(paste(bold("true 100 dimensional normal distribution"))), freq = FALSE)
  x_vals <- seq(min(sim100[,1]), max(sim100[,1]), length = 100)
  y_vals <- dnorm(x_vals, mean = mean(true100[,1]), sd = sd(true100[,1]))
  lines(x_vals, y_vals, col = "red", lwd = 2)

  # trace plots
  plot(sim10[1:1000, 1], type = "l", ylab = "x")
  plot(sim100[1:1000, 1], type = "l", ylab = "x")

  # assessing autocorrelation
  acf(sim10[, 1], main = "")
  acf(sim100[, 1], main = "")
}

plot_comparison(Gibbs10, Gibbs100, bnormal10, bnormal100)

## Comparison for easier covariance structures #################################

GenSigmaSimple <- function(d, sub_variance, variance) {
  sigma <- diag(variance, d)
  # Fill the superdiagonal
  for (i in 1:(d-1)) {
    sigma[i, i+1] <- sub_variance
  }

  # Fill the subdiagonal
  for (i in 2:d) {
    sigma[i, i-1] <- sub_variance
  }
  return(sigma)
}

sigma10 <- GenSigmaSimple(10, 0.5, 2)
sigma100 <- GenSigmaSimple(100, 0.5, 2)

Gibbs10 <- MVNormGibbs(1000, mu10, sigma10, rep(0, 10))
Gibbs100 <- MVNormGibbs(1000, mu100, sigma100, rep(0, 100))

## Benchmarking against the default normal #####################################
bnormal100 <- rmvnorm(1000, mu100, sigma100)
bnormal10 <- rmvnorm(1000, mu10, sigma10)

plot_comparison(Gibbs10, Gibbs100, bnormal10, bnormal100)

