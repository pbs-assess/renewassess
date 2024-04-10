## Nested models and AIC fun:

## Check AIC model choosing:
## True model y = b0 + b1*x
set.seed(1234)
n <- 100    ## Play around with sample size.
sigma <- 1.5  ## And standard deviation. 
x1 <- runif(n, -3, 3)
x2 <- runif(n, -3, 3)
x3 <- runif(n, -3, 3)
x4 <- runif(n, -3, 3)
x5 <- runif(n, -3, 3)
b0 <- 30
b1 <- -3  ## Smaller effect... what happens then?
options(na.action = "na.fail")
delta <- npar <- nested <- NULL

for( i in 1:1000 ){
  y <- b0 + b1*x1 + rnorm(n, 0, sd = sigma)
  fit <- lm(y ~ x1 + x2 + x3 + x4 + x5)
  ## Fit all combination of parameters with the very dubious dredge function.
  ## DO NOT RECOMMEND THIS FOR REAL DATA ANALYSIS.
  all.models <-  MuMIn::dredge(global.model = fit, fixed = "(Intercept)", rank = "AIC")  
  ## Which one is b0 + b1*x1?
  tru <- which(all.models$x1 != 0 & all.models$df == 3)
  ## delta AIC from top model.
  delta <- c(delta, all.models$AIC[1] - all.models$AIC[tru])
  
  npar <- c(npar, all.models$df[1]-1)
  ## Is this a nested model comparison?
  nested <- c(nested, !is.na(all.models$x1[1]))
}

boxplot(delta)

## How often was the correct model the best model?
mean(delta == 0)*100  ## Not great.

## How often was the correct model within 2?
mean(delta > -2)*100  ## Pretty good.

## Were the nested models close?
## Proportion of within 2 wrong answers that were close and not nested (e.g. top model didn't have x1 at all).
mean(delta > -2 & delta != 0 & !nested)*100 ## This split is very related to the choice of sigma.

## Proportion of within 2 wrong answers that were nested (included x2).
mean(delta > -2 & delta != 0 & nested)*100