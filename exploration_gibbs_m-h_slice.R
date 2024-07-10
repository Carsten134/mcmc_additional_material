## Gibbs vs m.h. exploration

## dependencies
library(rmcmc)
library(rbenchmark)
library(tidyverse)
library(stats)


## M-H. application
p <- function(x) {dnorm(x)}
q <- function(x, x_prime) {dnorm(x_prime, x)}
q_sample <- function(x) {rnorm(1,x)}

data_mh  = data.frame(x = metropolis_hastings(10, p, q, q_sample, 1000))

# diagnostics
plot(data_mh$x, type="l")
acf(data_mh$x)
hist(data_mh$x)


## Gibbs application
cond <- list(x1 = function(x) {rnorm(1, x)},
             x2 = function(x) {rnorm(1,x)})

init <- list(x1 = 10,
             x2 = 10)

data_g <- gibbs(cond, init, 1000)

# diagnostics
plot(data_g$x2, type="l")
acf(data_g$x1)
hist(data_g$x2)

## Bad initial values
# runtime explodes for M-H. no output prduceable
data_mh$x2 <- metropolis_hastings(100, p, q, q_sample, 1000)

init <- list(x1 = 100000,
             x2 = 100000)

data_g_2 <- gibbs(cond, init, 1000)

# diagnostics
plot(data_g_2$x2, type="l")

# very good burn in phase and very good autocorrelation
plot(data_g_2$x1[3:1000], type="l")
acf(data_g_2$x2[3:1000])


## Assessing runtime
benchmark_dif_p <- function(vals, func, rep) {
  # instanciate results
  results <- data.frame(elapsed=rep(-1, length(vals)), n=vals)
  for (i in 1:length(results$elapsed)) {
    # benchmark for different sample sizes
    temp <- benchmark(func(vals[i]), replications = rep, columns = c("elapsed"))
    results$elapsed[i] <- temp$elapsed[1]
  }
  return(results)
}

# M-H.
func <- function (n){metropolis_hastings(1,p, q, q_sample, n)}
vals <- c(10, 100, 200, 500, 1000)
result <- benchmark_dif_p(vals, func, 5)

plot(result$n, result$elapsed, type="l")

# gibbs
func <- function(n) {gibbs(cond, init, n)}

results_g <- benchmark_dif_p(vals, func, 10)

# Slice
func <- function(n) {OneDSlice(n, 10, 2, p)}

results_s <- benchmark_dif_p(vals, func, 10)

# comparing results
plot(result$n,result$elapsed, type = "n", xlim = c(1, max(vals)), ylim = c(0, max(result$elapsed)), 
     xlab = "Samplesize n", ylab = "elapsed in Sec.", main = "Gibbs vs M.H. vs. One-D-Slice runtime")
# lines(result$n, result$elapsed, col="red")
lines(results_g$n, results_g$elapsed, col="green")
lines(results_s$n, results_s$elapsed, col="blue")
lines(result$n, result$elapsed, col="red")
legend("bottomright", legend = c("M.H.", "Gibbs", "Slicing"), 
       col = c("red", "green", "blue"), lty = 1)
