## Multidimensional case:

## Hyperrectangle with shrinkage analogous to one-dimensional case

## Choose Hyperrectangle with appropiate sizes w_i in every dimension i and 
## position it randomly around x_0. Dont expand, as its difficult. Just 
## shrink until point within S is found.

## Maybe choose size w_i as  2 * sd_i? So it scales to width of distribution

## (alternative: Doubling Procedure with random sampling, that stops, when point
## outside Slice is sampled)

MVSLice <- function(it, x0, w, f, L.lim, R.lim) {
  n <- length(x0)
  x <- matrix(numeric(it * n) , nrow = it)
  x[1 ,] <- x0
  
  # browser()
  for(i in 2:it) {
    ## uniformly generate auxiliary variable
    u <- runif(1, min = 0, max = f(x[i - 1,])) 
    
    ## Determine sampling region:
    L <- x[i - 1,] - w * runif(n)
    L[L < L.lim] <- L.lim[L < L.lim]   ## Replace all components outside borders
    
    R <- L + w
    R[R > R.lim] <- R.lim[R > R.lim]   ## Replace all components outside borders
    
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



## Try out multivariate normal distribution:

## define f:
fnormal <- function(x) {
  mu = c(0, 0)
  Sigma = matrix(c(9, 8, 8, 9), nrow = 2)
  n <- length(x)
  fx <- 1/sqrt((2 * pi)^n * det(Sigma)) *
    exp(-0.5 * t(x - mu) %*% solve(Sigma) %*% (x - mu))
  return(fx)
}


mvnormaltest <- MVSLice(1000, x0 = c(0, 0), w = c(100, 100), f = fnormal, 
                        L.lim = c(-Inf, -Inf), R.lim = c(Inf, Inf))

plot(mvnormaltest, 
     col = rep(rainbow(10), each = 100))

plot(mvrnorm(n = 10000, mu = c(0, 0),
             Sigma  =  matrix(c(9, 8, 8, 9), nrow = 2)))


## First impressions: Choice of w is really important, w too small leads
## to very bad results. Larger w seems to be a lot better.
## c(1, 1) is ok, c(0.5, 0.5) already shows irregularities and c(0.1, 0.1)
## is a desaster (autocorrelation to high?)
## Maybe choose something like w = diag(Sigma^2)


## Multimodal distribution:
plot(rbind(mvrnorm(n = 1000, mu = c(0, 0),
                   Sigma  =  matrix(c(1, 0.7, 0.7, 1), nrow = 2)) ,
           mvrnorm(n = 1000, mu = c(0, 3),
                   Sigma  =  matrix(c(1, 0.7, 0.7, 1), nrow = 2))) )

fmult <- function(x) {
  mu1 = c(0, 0)
  mu2 <- c(0, 5)
  Sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2)
  n <- length(x)
  fx <- 1/sqrt((2 * pi)^n * det(Sigma)) *
    exp(-0.5 * t(x - mu1) %*% solve(Sigma) %*% (x - mu1)) + 
    1/sqrt((2 * pi)^n * det(Sigma)) *
    exp(-0.5 * t(x - mu2) %*% solve(Sigma) %*% (x - mu2))
  return(fx)
}

mvmulttest <- MVSLice(1000, x0 = c(0, 0), w = c(9, 9), f = fmult, 
                      L.lim = c(-Inf, -Inf), R.lim = c(Inf, Inf))

plot(mvmulttest, 
     col = rep(rainbow(10), each = 100))

## Here choice of w even more important, if w too small, regions of the
## distribution might never be reached! Factor in distance between peaks!
## In general larger w seems to always work better in terms of autocorrelation
## and accuracy of distribution. But might come at the expense of efficiency at
## higher dimensions.