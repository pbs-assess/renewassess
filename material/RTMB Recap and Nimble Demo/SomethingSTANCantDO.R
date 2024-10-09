## Choose sample size:
library(nimble)
load("harck")

modelcode <- nimbleCode({
  ## Make some priors. These give bounds for parameter transformations.
  beta ~ dgamma(1,1)
  sigma ~ dunif(0, 10)
  
  ## Bayesian non-parametric prior on log alpha.
  ## Creates an alpha that is drawn from a mixture of log alphas
  ## that are discretely chosen.
  mean_logalpha ~ dnorm(1.5, 1)
  mean_sd ~ dnorm(2.5, 0.5)
  for( i in 1:n ) {
    logalpha[i] ~ dnorm(mean=mean_logalpha, sd=mean_sd)
  }
  ## Sample it using the Chinese Restaurant Process
  id[1:n] ~ dCRP(conc = a, size = n)
  a ~ dgamma(shape = 1, rate = 1)

  ## Build the Likelihood
  for (i in 1:n){
    loga[i] <- logalpha[id[i]]
    E_logRS[i] <- loga[i] - beta*S[i]
    logRS[i] ~ dnorm(mean = E_logRS[i], sd = sigma)
  }
})
n <- nrow(harck)
constants <- list(n = nrow(harck), S = harck$S)
dat <- list(logRS = harck$logRS)
inits <- function(){list(id = 1:n, logalpha = 1+rnorm(n), a = 1, sigma = 3)}
## Create some initial values:
model <- nimbleModel(code = modelcode, constants = constants, data = dat, inits = inits())
cmodel <- compileNimble(model)

## Quick compare with MCMC version:
conf <- configureMCMC(model, monitors = c('id','loga', 'logalpha', 'beta', 'sigma', 'a', "mean_logalpha"), print = TRUE)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)
cmcmc$run(50000)
mvSamples <- cmcmc$mvSamples
samples <- as.matrix(mvSamples)[-(1:5000),] # for burn in.

## How many log alpha groups are there in the data?
nGroups <- apply(samples, 1, function(x)  length(unique(x)))
hist(nGroups)

## Look at some of the other posteriors
out <- coda::mcmc(samples)
mcmc.sum <- do.call("cbind", summary(out[, c("a", "sigma", "beta")]))
mcmc.sum

## Time series plot of logalpha to see what the actual impact is:
out.la <- coda::mcmc(samples[,grep("loga\\[", colnames(samples))])
mcmc.sum.la <- do.call("cbind", summary(out.la))
plot(harck$by, mcmc.sum.la[,"Mean"], type = "l", ylab = "Posterior mean log alpha", xlab = "Year")


## Let's see what years are most similar with a pairwise probability plot.
## We now know the probability from the model that two observations share the
## same productivity. Let's identify which years are the same.
idsamples <- samples[, grep('id', colnames(samples))]
pairMatrix <- apply(idsamples, 2, function(focal) {
                                   colSums(focal == idsamples)})
pairMatrix <- pairMatrix / nrow(samples)
ord <- ggcorrplot:::.hc_cormat_order(pairMatrix, hc.method ="complete") ## use hclust from gg to order them.
pairMatrix <- pairMatrix[ord, ord]
year <- harck$by[ord]
library(fields)
collist <- colorRampPalette(c('white', 'grey', 'black'))
image.plot(1:n, 1:n, pairMatrix , col = collist(6),
           xlab = 'Year', ylab = 'Year', axes = FALSE)
axis(1, at = 1:n, labels = year, tck = -.02, las = 2)
axis(2, at = 1:n, labels = year, tck = -.02, las = 2)
axis(3, at = 1:n, tck = 0, labels = FALSE)
axis(4, at = 1:n, tck = 0, labels = FALSE)
