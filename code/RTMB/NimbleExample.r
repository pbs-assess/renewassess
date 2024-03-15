###########################
## GLMM in Nimble with RTMB
###########################

library(nimble)
library(RTMB)

## Standard GLMM
code <- nimbleCode({
    beta0 ~ dnorm(0, sd = 10000)
    beta1 ~ dnorm(0, sd = 10000)
    sigma_RE ~ dunif(0, 1000)
    for (i in 1:N) {
        beta2[i] ~ dnorm(0, sd = sigma_RE)
        logit(p[i]) <- beta0 + beta1 * x[i] + beta2[i]
        r[i] ~ dbin(p[i], n[i])
    }
})

## constants, data, and initial values
constants <- list(N = 10)

data <- list(
    r = c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46),
    n = c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79),
    x = c(0,  0,  0,  0,  0,  1, 1,  1,  1,  1)
)

inits <- list(beta0 = 0, beta1 = 0, sigma_RE = 1)

## create the model object
glmmModel <- nimbleModel(code = code, constants = constants, data = data, 
                         inits = inits, check = FALSE)
glmmMCMC <- buildMCMC(glmmModel)
CglmmModel <- compileNimble(glmmModel)
CglmmMCMC <- compileNimble(glmmMCMC, project = glmmModel)
samples <- runMCMC(CglmmMCMC, niter = 10000, samplesAsCodaMCMC = TRUE)

## RTMB marginal likelihood:
dbinom_rtmb <- function(x){
  r <- x[1]
  n <- x[2]
  p <- plogis(x[3] + x[5])
  ans <- dnorm(x[5], 0, x[4], log=TRUE)
  ans <- ans + dbinom(r, prob = p, size = n, log = TRUE)
  -ans
}
fn <- MakeTape(dbinom_rtmb, c(0,1,0,0.1,0))
fnLaplace <- fn$laplace(random = 5)

fnLaplaceR <- nimbleRcall(function(x = double(1)){}, Rfun = 'fnLaplace',
                  returnType = double(0))
dbinMarg <- nimbleFunction(
  run = function(x=double(), size=double(), xbeta = double(), 
      sigma = double(), log = integer(0, default=0)){
    vec <- numeric(value = 0, length = 4)
    vec[1] <- x
    vec[2] <- size
    vec[3] <- xbeta
    vec[4] <- sigma
    
    ans <- fnLaplaceR(vec)
    if(is.nan(ans)) ans <- -Inf
    returnType(double())
    if(log) return(-ans) else return(exp(-ans))
  })

## Check compiles:
cdbin <- compileNimble(dbinMarg)
cdbin(1,2,1,1,1)

## Standard GLMM
codeRTMB <- nimbleCode({
    beta0 ~ dnorm(0, sd = 10000)
    beta1 ~ dnorm(0, sd = 10000)
    sigma_RE ~ dunif(0, 1000)
    for (i in 1:N) {
        xbeta[i] <- beta0 + beta1 * x[i]
        r[i] ~ dbinMarg(size=n[i], xbeta=xbeta[i], sigma = sigma_RE)
    }
})

constants <- list(N = 10)
data <- list(
    r = c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46),
    n = c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79),
    x = c(0,  0,  0,  0,  0,  1, 1,  1,  1,  1)
)

inits <- list(beta0 = 0, beta1 = 0, sigma_RE = 0.3)

## create the model object
glmmModelMarg <- nimbleModel(code = codeRTMB, constants = constants, data = data, 
                         inits = inits, check = FALSE)
glmmMargMCMC <- buildMCMC(glmmModelMarg)
CglmmModelMarg <- compileNimble(glmmModelMarg)
CglmmMargMCMC <- compileNimble(glmmMargMCMC, project = glmmModelMarg)
samplesMarg <- runMCMC(CglmmMargMCMC, niter = 10000,  samplesAsCodaMCMC = TRUE)


coda::traceplot(samples[,'sigma_RE'])
lines(samplesMarg[,'sigma_RE'], col='red')

coda::traceplot(samples[,'beta0'])
lines(samplesMarg[,'beta0'], col='red')

coda::traceplot(samples[,'beta1'])
lines(samplesMarg[,'beta1'], col='red')