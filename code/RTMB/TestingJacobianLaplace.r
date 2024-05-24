library(RTMB)

vals <- vals0 <- NULL
for( ii in 1:1000 ){
  n <- 100
  U <- rgamma(100, shape = 1, rate = 0.1)
  J <- 200
  y <- matrix(0, nrow = n, ncol = J)
  for( i in 1:n){
    y[i,] <- rpois(J, U[i])
  }
  # hist(y)

  pars <- list()
  pars$z <- log(U)
  pars$logshape <- 0
  pars$lograte <- log(0.1)

  data <-  list()
  data$y <- y
  data$n <- n

  fn <- function(pars){
    getAll(pars, data)
    negll <- 0
    shape <- exp(logshape)
    rate <- exp(lograte)
    u <- exp(z)
    negll <- negll - sum(dgamma(u, shape = shape, scale = 1/rate, log = TRUE)) - sum(z)
    for( i in 1:n ) negll <- negll - sum(dpois(y[i,], u[i], log = TRUE))
    negll
  }

  obj <- MakeADFun(fn, pars, random = "z", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  vals <- rbind(vals, exp(opt$par))

  fnnojac <- function(pars){
    getAll(pars, data)
    negll <- 0
    shape <- exp(logshape)
    rate <- exp(lograte)
    u <- exp(z)
    negll <- negll - sum(dgamma(u, shape = shape, scale = 1/rate, log = TRUE))
    for( i in 1:n ) negll <- negll - sum(dpois(y[i,], u[i], log = TRUE))
    negll
  }
  
  ## Clearly the jacobian is working better given the need for a trycatch here.
  objnojac <- MakeADFun(fnnojac, pars, random = "z", silent = TRUE)
  optnojac <- tryCatch({ nlminb(objnojac$par, objnojac$fn, objnojac$gr)},  warning=function(w){ list(par = c(-Inf, -Inf))},
    error = function(e) {list(par = c(-Inf, -Inf))})
  optnojac$par
  vals0 <- rbind(vals0, exp(optnojac$par))
}
par(mfrow = c(2,1))
boxplot(vals[,1], vals0[vals0[,1] != 0,1])  ## Get rid of those that didn't converge.
abline(h = 1, col = 'red')
boxplot(vals[,2], vals0[vals0[,1] != 0,2])
abline(h = 0.1, col = 'red')

mean(vals0[,1] != 0)  ## Sometimes it was bad.

mean(vals[,1] - 1)
mean(vals0[vals0[,1] != 0,1] - 1) # Slightly more biased.
mean(vals[,2] - 0.1)
mean(vals0[vals0[,2] != 0,2] - 0.1) # Slightly more biased.

## Can actually write what the marginal distribution in this case because the gamma is a conjugate prior to the poisson.
## I'm lazy but you can check how well the Laplace approx is even doing by writing that out.

