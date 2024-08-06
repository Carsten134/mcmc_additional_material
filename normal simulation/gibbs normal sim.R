## Dependencies ################################################################
library(condMVNorm)

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
