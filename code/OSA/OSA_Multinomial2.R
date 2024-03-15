#######################################################
## Basic multinomial understanding check for OSA:
## Paul van Dam-Bates
## 3-15-2024
##
## SA: revised to match TMB. (1) use the TMB log calculations, (2) include the
## last category to match (needed? I think just so seed matches),
# (3) importantly: need to accumulate nll from previous observations
## from previously entered observations in OSA.
## Code is a bit of a mess now... but it matches :)
#######################################################

## Install Trijoulet package code:
# TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
# remotes::install_github("fishfollower/compResidual/compResidual", force=TRUE)

## Basic Example from R package vignette:
library(compResidual)

X <- matrix(c(
  9,    1,    5,    9,    9,    3,    6,    3,    3,     4,
  1,    5,    7,    1,    3,    6,    8,    9,    6,     5,
  3,    7,    4,    2,    1,    8,    4,    5,    1,     2,
  4,    2,    3,    5,    3,    1,    0,    1,    0,     2,
  1,    3,    3,    0,    5,    5,    4,    3,    1,     0,
  7,    7,    3,    8,    4,    2,    3,    4,   14,    12
), nrow = 6, ncol = 10, byrow = TRUE)

P <- matrix(c(
  0.32, 0.08, 0.16, 0.24, 0.32, 0.20, 0.20, 0.16, 0.16, 0.16,
  0.16, 0.16, 0.24, 0.20, 0.12, 0.16, 0.32, 0.28, 0.20, 0.20,
  0.12, 0.24, 0.20, 0.12, 0.04, 0.24, 0.16, 0.12, 0.04, 0.08,
  0.04, 0.16, 0.16, 0.12, 0.04, 0.08, 0.00, 0.08, 0.04, 0.20,
  0.12, 0.08, 0.12, 0.00, 0.12, 0.12, 0.12, 0.08, 0.04, 0.04,
  0.24, 0.28, 0.12, 0.32, 0.36, 0.20, 0.20, 0.28, 0.52, 0.32
), nrow = 6, ncol = 10, byrow = TRUE)

set.seed(123)
## Residuals from package
res <- resMulti(X, P)

## Each column is one sample from a multinomial distribution:
## This is my understanding of OSA for multinomials, based on a
## series of binomials and RQR transformation to z.

# 123 matches TMB::oneStepPredict seed, *but* it's hardcoded to NULL
# in resMulti() function!
rqrMulti <- function(X, P, seed = 123) {
  nll <- NULL
  row <- NULL
  col <- NULL
  px <- NULL
  Fx <- NULL
  K <- nrow(X)
  N <- ncol(X)
  nll_out <- NULL
  nlcdf.lower_out <- NULL
  nlcdf.upper_out <- NULL
  for (i in 1:N) {
    ## Equation for Multinomial to Binomials in Trijoulet paper:
    for (j in 1:(K)) {
      if (j == 1) {
        p <- P[1, i]
        n <- sum(X[, i])
      } else {
        p <- P[j, i] / (1 - sum(P[1:(j - 1), i]))
        n <- sum(X[j:K, i])
      }
      suppressWarnings({ # Infs/NaNs, dealt with below
        nll <- -dbinom(X[j, i], size = n, prob = p, log = TRUE)
        nll_out <- c(nll_out, nll)
        nll_up_to_obs <- sum(nll_out) - nll
        nlcdf.lower <- -pbinom(X[j, i], size = n, prob = p, log.p = TRUE) + nll_up_to_obs
        nlcdf.upper <- -pbinom(X[j, i], size = n, prob = p, log.p = TRUE, lower.tail = FALSE) +
          nll_up_to_obs
      })
      Fx <- c(Fx, 1 / (1 + exp(nlcdf.lower - nlcdf.upper))) # CDF of obs conditional on...
      px <- c(px, 1 / (exp(-nlcdf.lower + sum(nll_out)) +
          exp(-nlcdf.upper + sum(nll_out)))) # density of obs conditional on...
      if (is.nan(nll)) { # generally these are the last category
        nll <- 0
        nll_out[length(nll_out)] <- 0
        Fx[length(Fx)] <- 1
        px[length(px)] <- 1
      }
      nlcdf.upper_out <- c(nlcdf.upper_out, nlcdf.upper)
      nlcdf.lower_out <- c(nlcdf.lower_out, nlcdf.lower)
      row <- c(row, j)
      col <- c(col, i)
    }
  }
  set.seed(seed)
  U <- runif(length(Fx))
  z <- qnorm(Fx - U * px)
  data.frame(row = row, col = col, residual = z, Fx = Fx, px = px,
    nll = cumsum(nll_out), nlcdf.upper = nlcdf.upper_out, nlcdf.lower = nlcdf.lower_out)
}

## Our R function version:
res2 <- rqrMulti(X, P)
res2 <- res2[-seq(6, 60, by = 6), ] # last category

## Compare them:
qqnorm(res, pch = 4, col = 'blue')
qqline(res)
qqnorms <- qqnorm(res2$residual, plot.it = FALSE)
points(qqnorms$x, qqnorms$y, col = 'red', pch = 16)

## Same:
plot(as.numeric(res), res2$residual)
abline(0, 1)

max(as.numeric(res) - res2$residual)
