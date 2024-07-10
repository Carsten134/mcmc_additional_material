## Monte Carlo Simulations: Theory and Practise
## Markov Chain Monte Carlo


## Slice Sampling:
## Implementation of the Slicing Algorithm:

## First try using inverse density for determining the slice:
simpleslice <- function(it, x0, f, inv.l, inv.u) {
  x <- numeric(it)
  x[1] <- x0
  for(i in 2:it) {
    u <- runif(1, min = 0, max = f(x[i - 1]))

    ## Inverse has to be more flexible
    x[i] <- runif(1, min = inv.l(u), max = inv.u(u))
  }

  return(x)
}


## Try it out:
## Standard normal distribution with bad initial value:
set.seed(74765)
x.firsttest <- simple_slice(1000, x0 = 15, f = function(x) {exp(-x^2/2)/sqrt(2 * pi)},
                     inv.l = function(x) {- sqrt(-2 * log(x * sqrt(2 * pi)))},
                     inv.u = function(x) {sqrt(-2 * log(x * sqrt(2 * pi)))})

hist(x.firsttest[-1])
which.max(x.firsttest[-1])
x.firsttest[1:10]
## really fast burn-in. only first 3 samples are unlikely for target
## distribution.
## Better choice of x0 results in faster burn-in


## Try out another example from Roberts and Casella: exp-distribution:
x.exp <- simpleslice(10000, x0 = 1, f = function(x) {0.5 * exp(-sqrt(x))},
                     inv.l = function(x) {0},
                     inv.u = function(x) {log(2 * x)^2})

hist(x.exp[x.exp < 70], breaks = 0:70)
## looks good.


## What about more complex f? What if Slice is not continuous?

## expansion or doubling techniques. Important: Initial interval needs to
## be random: L = x_0 - w * runif(1), R = L + w.  (w interval length)
## Interval constructed around x[i - 1], as this is guaranteed to lie in
## interval.

## How to choose appropriate w? need a way to approximate size of slice
## Doubling may help here as it is faster, when w is chosen too small
## -> biggest problem so far.
## Also should w be adaptive to f(x)? Maybe just in dependence of sd?
## Maybe like 1/8 * sd?

## Input:
## it: number of iterations
## x0: starting value
## w: tuning parameter: Initial slice size
## f: function of x, density of desired distribution
## L.lim, R.lim: lower and upper limit of the slice. [-Inf, Inf] is default,
##               but has to be specified for something like exponential dist.
##               (Then L.lim = 0 for example)


OneDSlice <- function(it, x0, w, f, L.lim = -Inf, R.lim = Inf) {
  #browser()
  ## Vector of simulated data:
  x <- numeric(it)
  x[1] <- x0

  for(i in 2:it) {
    ## uniformly generate auxiliary variable
    u <- runif(1, min = 0, max = f(x[i - 1]))

    ## Determine sampling interval:
    L <- x[i - 1] - w * runif(1)
    if(L < L.lim) {
      L <- L.lim
    }

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

    ## sample from [L, R]. Accept, if x_star is in [L, R], else shrink interval
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


## Running examples:

## Normal distribution:
xnormal <- OneDSlice(10000, x0 = 1, w = 1/8,
                   f = function(x) {exp(-x^2/2)/sqrt(2 * pi)})

hist(xnormal)
## Works.

## Multimodal case:
xmm <- OneDSlice(10000, x0 = 1, w = 1/8,
                 f = function(x) {exp(-x^2/2)/sqrt(2 * pi) +
                                  exp(-(x - 4)^2/2)/sqrt(2 * pi)})

hist(xmm, freq = FALSE, ylim = c(0, 0.25), breaks = 100)
curve((exp(-x^2/2)/sqrt(2 * pi) +
         exp(-(x - 4)^2/2)/sqrt(2 * pi))/2,
      from = -4, to = 8, add = TRUE)
## works as well.

xmm2 <- OneDSlice(10000, x0 = 1, w = 1/8,
                  f = function(x) {exp(-x^2/2)/sqrt(2 * pi) +
                      exp(-(x - 4)^2/2)/sqrt(2 * pi) +
                      exp(-((x + 5)/3)^2/2)/sqrt(2 * pi * 9)})

hist(xmm2, freq = FALSE, ylim = c(0, 0.15), breaks = 100)
curve((exp(-x^2/2)/sqrt(2 * pi) +
        exp(-(x - 4)^2/2)/sqrt(2 * pi)+
        exp(-((x + 5)/3)^2/2)/sqrt(2 * pi * 9))/3,
      from = -15, to = 8, add = TRUE)
## Looks perfect.

## Try out case with bad initial value:
##...

## What if distribution is bounded, like exponential distribution?
## In this case we need to set the lower limit for the slice L.lim = 0.

xexp <- OneDSlice(10000, x0 = 1, w = 1/8,
                  f = function(x) {0.5 * exp(-x * 0.5)},
                  L.lim = 0)

hist(xexp, freq = FALSE, breaks = 100)
curve(0.5 * exp(-x * 0.5), from = 0, to = 10, add = TRUE)



## Diagnostics: (Autocorrelation etc)
## ...


## Multidimensional case:

## Hyperrectangle with shrinkage analogous to one-dimensional case

## Choose Hyperrectangle with appropiate sizes w_i in every dimension i and
## position it randomly around x_0. Dont expand, as its difficult. Just
## shrink until point within S is found.

## Maybe choose size w_i as  2 * sd_i? So it scales to width of distribution

## (alternative: Doubling Procedure with random sampling, that stops, when point
## outside Slice is sampled)

MVSLice <- function(it, x0, w, f) {
  x <- numeric(it)
  return(x)
}


## Advantages and disadvantages:

## + Only needs sampling from uniform distribution
## + Only needs knowledge of density f(x)
## + Only one tuning parameter (w)

## - Efficiency dependent on choice of w. supervision and knowledge necessary
## -

