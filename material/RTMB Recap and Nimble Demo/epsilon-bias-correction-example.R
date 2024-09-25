library(RTMB)

set.seed(21453)
dat <- NULL
re.true <- NULL
sd.re <- 6
lambda <- 1
loglambda <- log(lambda)
for (i in 1:100) {
  ngrpi <- rpois(1, 10)
  re <- rnorm(1, 0, sd.re)
  dat <- rbind(dat, data.frame(y = rpois(ngrpi, exp(re + loglambda)), id = i))
  re.true <- c(re, re.true)
}

nll <- function(pars) {
  getAll(pars, dat)
  sdre <- exp(logsdre)
  lambda <- exp(loglambda)
  ADREPORT(lambda)
  ADREPORT(sdre)

  negll <- -sum(dnorm(re, 0, sdre, log = TRUE))

  # this demonstrates the 'epsilon' bias correction:
  # https://doi.org/10.1016/j.fishres.2015.11.016
  mu <- exp(loglambda + re[id])
  ADREPORT(mu[1])
  ADREPORT(mu[20])

  negll <- negll - sum(dpois(y, mu, log = TRUE))
  negll
}
obj <- MakeADFun(nll, parameters = list(
  logsdre = 0, loglambda = 0,
  re = rep(0, max(dat$id))
), silent = TRUE, random = "re")
fit1 <- nlminb(obj$par, obj$fn, obj$gr)
sdrep1 <- sdreport(obj)
summary(sdrep1, "report")

sdrep2 <- sdreport(obj, bias.correct = TRUE)
summary(sdrep2, "report")

sdrep2 <- sdreport(obj, bias.correct = TRUE, bias.correct.control = list(sd = TRUE))
summary(sdrep2, "report")
