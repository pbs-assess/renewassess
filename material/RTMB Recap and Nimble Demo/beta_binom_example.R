## Choose sample size:
set.seed(123)
n <- 50
N <- 20
p <- rbeta(n, 10, 2)
y <- rbinom(n, size = N, prob = p)

fn <- function(pars){
  getAll(pars, dat)
  a <- exp(loga)
  b <- exp(logb)
  n <- length(y)
  negll <- 0
  
  p <- 1/(1+exp(-re))
  for( i in 1:n )  {
    negll <- negll - dbeta(p[i], a, b, log = TRUE) - log(exp(-p[i])/(1+exp(-p[i]))^2)
    negll <- negll - dbinom(y[i], prob = p[i], size = N, log = TRUE)
  }
  return(negll)
}

pars <- list(loga = 0, logb = 0, re = rep(0, n))
dat <- list(y = y)

obj <- MakeADFun(fn, pars, random = "re", silent = TRUE)

ll.betabin <- function(pars){
  a <- exp(pars[1])
  b <- exp(pars[2])
  ll <- 0
  for( i in seq_along(dat$y))
  {
    ll <- ll - lbeta(a,b) + lchoose(N, dat$y[i]) + lbeta(a + dat$y[i], b + N-dat$y[i])
  }
  ll
}

dibeta <- function(a,b)
{
  da <-  digamma(a) - digamma(a+b) 
  db <-  digamma(b) - digamma(a+b)
  c(da, db)
}

gr.betabin <- function(pars){
  a <- exp(pars[1])
  b <- exp(pars[2])
  dll <- 0
  for( i in seq_along(dat$y))
  {
      dll <- dll - dibeta(a,b) + dibeta(a + dat$y[i], b + N-dat$y[i])
  }
  return(dll*c(a,b))
}

fit <- nlminb(obj$par, obj$fn, obj$gr)
fit.true <- nlminb(fit$par, ll.betabin, gr.betabin)
