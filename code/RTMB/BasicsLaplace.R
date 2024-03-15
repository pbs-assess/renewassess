## Understanding Laplace:

library(RTMB)

## Fake data:
x <- rnorm(100, 0, 10)
y <- rnorm(100, 0.5 + 1.5*x, 0.5)
grps <- sample(7, 100, replace = TRUE)
for(i in 1:7) y[grps == i] <- y[grps == i] + rnorm(1, 0, 0.25)

dat <- list()
dat$x = x
dat$y = y
dat$grps = grps

pars <- list()
pars$beta0 <- 0
pars$beta1 <- 0
pars$re <- rep(0, 7)
pars$logsd <- 0
pars$logsd.re <- 0

f <- function(pars){
  getAll(pars, dat)
  sd.re <- exp(logsd.re)
  sd <- exp(logsd)
  negLL <- -sum(dnorm(re, 0, sd.re, log = TRUE))
  mu <- beta0 + beta1*x + re[grps]
  negLL <- negLL - sum(dnorm(y, mu, sd, log = TRUE))
  return(negLL)
}

## Now let's try and do this ourselves:
## Need to make only pars random-effects
fmarg <- function(data, pars){
  finternal <- function(pars){
    getAll(pars, dat)
    sd.re <- exp(logsd.re)
    sd <- exp(logsd)
    negLL <- -sum(dnorm(re, 0, sd.re, log = TRUE))
    mu <- beta0 + beta1*x + re[grps]
    negLL <- negLL - sum(dnorm(y, mu, sd, log = TRUE))
    return(negLL)
  }

  objin <- MakeADFun(finternal, pars, silent = TRUE, 
    map = list(logsd = factor(NA), logsd.re = factor(NA), beta0 = factor(NA), beta1 = factor(NA)))
  opt <- nlminb(objin$par, objin$fn, objin$gr)
  Q <- objin$he()
  ## Laplace:  
  ## log like
  ans <- -objin$fn(opt$par) - 0.5 * log(det(Q)) + 0.5 * length(pars$re) * log(2*pi)
  return(-ans)
}

## Use Laplace to marginalize over the random effect
obj <- MakeADFun(f, pars, random = "re", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
summary(sdreport(obj))
opt$objective

## Do it manually with a nested RTMB function to get their AD still.
pl <- obj$env$parList()
pl$re <- rep(0, length(pl$re))
opt$objective - fmarg(dat, pl)  ## Pretty cool eh?
