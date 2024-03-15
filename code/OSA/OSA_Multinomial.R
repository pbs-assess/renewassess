#######################################################
## Basic multinomial understanding check for OSA:
## Paul van Dam-Bates
## 3-15-2024
#######################################################

## Install Trijoulet package code:
# TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
# remotes::install_github("fishfollower/compResidual/compResidual", force=TRUE)

## Basic Example from R package vignette:
library(compResidual)

X<-matrix(c(
  9,    1,    5,    9,    9,    3,    6,    3,    3,     4,
  1,    5,    7,    1,    3,    6,    8,    9,    6,     5,
  3,    7,    4,    2,    1,    8,    4,    5,    1,     2,
  4,    2,    3,    5,    3,    1,    0,    1,    0,     2,
  1,    3,    3,    0,    5,    5,    4,    3,    1,     0,
  7,    7,    3,    8,    4,    2,    3,    4,   14,    12
), nrow=6, ncol=10, byrow=TRUE) 

P<-matrix(c(
  0.32, 0.08, 0.16, 0.24, 0.32, 0.20, 0.20, 0.16, 0.16, 0.16,
  0.16, 0.16, 0.24, 0.20, 0.12, 0.16, 0.32, 0.28, 0.20, 0.20,
  0.12, 0.24, 0.20, 0.12, 0.04, 0.24, 0.16, 0.12, 0.04, 0.08,
  0.04, 0.16, 0.16, 0.12, 0.04, 0.08, 0.00, 0.08, 0.04, 0.20,
  0.12, 0.08, 0.12, 0.00, 0.12, 0.12, 0.12, 0.08, 0.04, 0.04,
  0.24, 0.28, 0.12, 0.32, 0.36, 0.20, 0.20, 0.28, 0.52, 0.32
), nrow=6, ncol=10, byrow=TRUE) 

## Residuals from package
res <- resMulti(X,P)


## Each column is one sample from a multinomial distribution:
## This is my understanding of OSA for multinomials, based on a
## series of binomials and RQR transformation to z.
rqrMulti <- function(X,P){
  
  z <- NULL
  row <- NULL
  col <- NULL
  K <- nrow(X)
  N <- ncol(X)
  
  for( i in 1:N ){
    ## Equation for Multinomial to Binomials in Trijoulet paper:
    for( j in 1:(K-1) ){
      if( j == 1 ){
        p <- P[1,i]
        n <- sum(X[,i])
      }else{
        p <- P[j,i]/( 1 - sum(P[1:(j-1),i]) )
        n <- sum( X[j:K,i] )
      }
      ## Now the RQR classic part for adjusting discrete CDF (Dunn and Smyth 1996):
      b <- pbinom(X[j,i], size = n, prob = p)  # CDF at x
      a <- pbinom(X[j,i]-1, size = n, prob = p) ## CDF at x-1
      z <- c( z, qnorm(runif(1, a, b)) )  ## inverse CDF at unif(a,b) is the quantile.
      row <- c(row, j)
      col <- c(col, i)
    }
  }
  data.frame(row = row, col = col, residual = z)
}


## Paul's function
res2 <- rqrMulti(X,P)

## Compare them:
qqnorm(res, pch = 4, col = 'blue')
qqline(res)
qqnorms <- qqnorm(res2$residual, plot.it = FALSE)
points(qqnorms$x, qqnorms$y, col = 'red', pch = 16)

## Nearly the same...
plot(as.numeric(res), res2$residual)  ## A couple are very different.