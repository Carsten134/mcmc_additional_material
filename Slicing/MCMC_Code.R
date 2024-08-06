
## Working directory:

setwd("~/Studium/Master Statistik/Monte Carlo Simulations/Programming")

## Packages:

library("mvtnorm")



################################################################################
## Slice Sampling:
################################################################################


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
    
    ## Determine sampling region of size w:
    L <- x[i - 1,] - w * runif(n)
    R <- L + w
    
    L[L < L.lim] <- L.lim[L < L.lim]   ## Set border to L.lim if necessary 
    R[R > R.lim] <- R.lim[R > R.lim]   ## Set border to R.lim if necessary 
    
    ## sample from [L, R]. Accept, if x_star is in [L, R], else shrink region
    ## to coordinates of sample
    repeat {
      xstar <- sapply(1:n, FUN = function(y) {runif(1, min = L[y], max = R[y])})
      
      if(f(xstar) > u) {  ## acceptance condition
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
## Application to higher-dimensional synthetic data
################################################################################

## Try out Slicing with synthetic data of 10 variables:

## define f using the function dmvnorm from the package mvtnorm:
## use random covariances:

mu10 <- rep(0, 10)
cov10 <- matrix(sample(seq(0.1, 0.9, by = 0.1), 100, replace = TRUE), 
                nrow = 10)

## make it symmetric:
cov10[lower.tri(cov10)] <- t(cov10)[lower.tri(cov10)]
diag(cov10) <- rep(6, 10)

## make it positive semidefinite, if necessary
make_positive_semidefinite <- function(mat) {
  eig <- eigen(mat)
  eig$values[eig$values < 0] <- 0
  mat_psd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  return(mat_psd)
}

cov10 <- make_positive_semidefinite(cov10)
cov10
det(cov10) ## high determinant

## define density:
f10 <- function(x, mean = mu10, cov = cov10) {
  return(dmvnorm(x, mean = mean, sigma = cov10))
}


## Try out at different points:
f10(rep(0, 10))
sapply(seq(from = -2, to = 2, by = 0.1), 
       FUN = function(x) {f10(rep(x, 10))} )
## Seems to work fine.


## Apply slicing algorithm:
Slice10 <- MVSlice(10000, rep(0, 10), w = rep(50, 10), f10, 
                   L.lim = rep(-Inf, 10), 
                   R.lim = rep(Inf, 10))

## Examine results:

## Plot histogram of distribution of one variable, trace plot of first 1000 
## observations and ACF-Plot.
par(mfrow = c(3, 1), mar = c(4, 4, 1.5, 1))
hist(Slice10[, 1], breaks = 50, xlab = "x", 
     main = substitute(paste(bold("Histogram, Trace plot and ACF from 10-dim. normal distr."))))
plot(Slice10[1:1000, 1], type = "l", ylab = "x")
acf(Slice10[, 1], main = "")

## Try another variable
hist(Slice10[, 5], breaks = 50, xlab = "x",
     main = substitute(paste(bold("Histogram, Trace plot and ACF from 10-dim. normal distr."))))
plot(Slice10[1:1000, 5], type = "l", ylab = "x")
acf(Slice10[, 5], main = "")

## Higher Autocorrelation than in the univariate case, but not worse than 
## two-dimensional case.



## Lets try if algorithm breaks down for higher dimensions:

## 100 Dimensional Covariance Matrix:
mu100 <- rep(0, 100)
cov100 <- matrix(sample(seq(0.1, 0.9, by = 0.1), 10000, replace = TRUE), 
                 nrow = 100)

## make it symmetric:
cov100[lower.tri(cov100)] <- t(cov100)[lower.tri(cov100)]
diag(cov100) <- rep(6, 100)

## make it positive semidefinite, if necessary
make_positive_semidefinite <- function(mat) {
  eig <- eigen(mat)
  eig$values[eig$values < 0] <- 0
  mat_psd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  return(mat_psd)
}

cov100 <- make_positive_semidefinite(cov100)
cov100  
det(cov100) 

## Define density:
f100 <- function(x, mean = mu100, cov = cov100) {
  return(dmvnorm(x, mean = mean, sigma = cov))
}



## Apply slicing algorithm, use 10000 iterations again, it takes a while longer,
## but is manageable
Slice100 <- MVSlice(10000, rep(0, 100), w = rep(50, 100), f100, 
                    L.lim = rep(-Inf, 100), 
                    R.lim = rep(Inf, 100))

## Examine results again:

## Plot histogram of distribution of one variable, trace plot of first 1000 
## observations and ACF-Plot.
hist(Slice100[, 1], breaks = 50, xlab = "x", 
     main = substitute(paste(bold("Histogram, Trace plot and ACF from 100-dim. normal distr."))))
plot(Slice100[1:1000, 1], type = "l", ylab = "x")
acf(Slice100[, 1], main = "")

## Try another variable
hist(Slice100[, 56], breaks = 50, xlab = "x",
     main = substitute(paste(bold("Histogram, Trace plot and ACF from 100-dim. normal distr."))))
plot(Slice100[1:1000, 56], type = "l", ylab = "x")
acf(Slice100[, 56], main = "")

## Much higher autocorrelation, Trace Plot makes much smaller jumps, and 
## distribution in histogram looks slightly worse. The distribution spreads
## larger and is not as smooth as it could be



## Compare 10-dim and 100-dim case in one single figure:
par(mfrow = c(3, 2), mar = c(4, 4.2, 1.5, 1),
    cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(Slice10[, 1], breaks = 50, xlab = "x", xlim = c(-6, 6),
     main = substitute(paste(bold("10 dimensional normal distribution"))))
hist(Slice100[, 1], breaks = 50, xlab = "x", xlim = c(-6, 6),
     main = substitute(paste(bold("100 dimensional normal distribution"))))

plot(Slice10[1:1000, 1], type = "l", ylab = "x")
plot(Slice100[1:1000, 1], type = "l", ylab = "x")

acf(Slice10[, 1], main = "")
acf(Slice100[, 1], main = "")



## Can the high dimensional case handle bad starting values?

Slice10Bad <- MVSlice(10000, rep(10, 10), w = rep(50, 10), f10, 
                   L.lim = rep(-Inf, 10), 
                   R.lim = rep(Inf, 10))

par(mfrow = c(3, 1), mar = c(4, 4, 1.5, 1))
hist(Slice10Bad[, 1], breaks = 50, xlab = "x", 
     main = substitute(paste(bold("10 dimensional normal distr. with bad start"))))
plot(Slice10Bad[1:1000, 1], type = "l", ylab = "x")
acf(Slice10Bad[, 1], main = "")
## does not look too bad, fast burn-in, autocorrelation similar


Slice100Bad <- MVSlice(10000, rep(10, 100), w = rep(50, 100), f100, 
                    L.lim = rep(-Inf, 100), 
                    R.lim = rep(Inf, 100))

par(mfrow = c(2, 1), mar = c(4, 4, 1.5, 1))
hist(Slice100Bad[, 1], breaks = 50, xlab = "x", 
     main = substitute(paste(bold("100 dimensional normal distribution with bad starting value"))))
plot(Slice100Bad[1:1000, 1], type = "l", ylab = "x", 
     ylim = c(-3, 10))
abline(a = 0, b = 0, col = "red", lty = 2)
acf(Slice100Bad[, 1], main = "")
## Notice how we stay above 0 for 1000 iterations! 
## Really long burn-in, we need a lot of samples to get a decent distribution

par(mfrow = c(1, 1))
plot(Slice100Bad[1:1000, 1], type = "l", ylab = "x", 
     ylim = c(-5, 10))
abline(a = 0, b = 0, col = "red", lty = 2)



## Calculate Correlation matrix:
round(cov2cor(cov10), 2)
round(cov2cor(cov100), 2)



################################################################################
## Try out simpler Covariance Matrix:

## 10 dimensions:
sigma10 <- diag(1, 10, 10)

for (i in 1:(10 - 1)) {     # Fill the superdiagonal
  sigma10[i, i+1] <- 0.5
}
for (i in 2:10) {           # Fill the subdiagonal
  sigma10[i, i-1] <- 0.5
}

## Define density:
f10simple <- function(x, mean = mu10, cov = sigma10) {
  return(dmvnorm(x, mean = mean, sigma = cov))
}

## Apply slicing algorithm:
Slice10simple <- MVSlice(10000, rep(0, 10), w = rep(50, 10), f10simple, 
                         L.lim = rep(-Inf, 10), 
                         R.lim = rep(Inf, 10))



## 100 dimensions:
sigma100 <- diag(1, 100, 100)

for (i in 1:(100 - 1)) {     # Fill the superdiagonal
  sigma100[i, i+1] <- 0.5
}
for (i in 2:100) {           # Fill the subdiagonal
  sigma100[i, i-1] <- 0.5
}

## Define density:
f100simple <- function(x, mean = mu100, cov = sigma100) {
  return(dmvnorm(x, mean = mean, sigma = cov))
}

## Apply slicing algorithm:
Slice100simple <- MVSlice(10000, rep(0, 100), w = rep(50, 100), f100simple, 
                         L.lim = rep(-Inf, 100), 
                         R.lim = rep(Inf, 100))

## Plot:
par(mfrow = c(3, 2), mar = c(4, 4.2, 1.5, 1),
    cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(Slice10simple[, 1], breaks = 50, xlab = "x", xlim = c(-3, 3),
     main = substitute(paste(bold("10 dimensional normal distribution"))))
hist(Slice100simple[, 1], breaks = 50, xlab = "x", xlim = c(-3, 3),
     main = substitute(paste(bold("100 dimensional normal distribution"))))

plot(Slice10simple[1:1000, 1], type = "l", ylab = "x")
plot(Slice100simple[1:1000, 1], type = "l", ylab = "x")

acf(Slice10simple[, 1], main = "")
acf(Slice100simple[, 1], main = "")



