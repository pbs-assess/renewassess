``` r
# special branch:
# devtools::install_github("pbs-assess/sdmTMB", ref = "lognormal")

library(sdmTMB)
d <- subset(pcod_2011, density > 0)

# with mean_adjust = TRUE (only option on main branch)
fit <- sdmTMB(
  density ~ 1,
  data = d, spatial = "off",
  family = lognormal(link = "log")
)
fit
#> Model fit by ML ['sdmTMB']
#> Formula: density ~ 1
#> Mesh: NULL (isotropic covariance)
#> Data: d
#> Family: lognormal(link = 'log')
#>  
#>             coef.est coef.se
#> (Intercept)      4.5    0.09
#> 
#> Dispersion parameter: 1.41
#> ML criterion at convergence: 2350.450
#> 
#> See ?tidy.sdmTMB to extract these values as a data frame.

p <- exp(predict(fit)$est[1])
p
#> [1] 90.38351
mean(d$density)
#> [1] 88.52506

# now with mean_adjust = FALSE

fit <- sdmTMB(
  density ~ 1,
  data = d, spatial = "off",
  family = lognormal(link = "log", mean_adjust = FALSE)
)
fit
#> Model fit by ML ['sdmTMB']
#> Formula: density ~ 1
#> Mesh: NULL (isotropic covariance)
#> Data: d
#> Family: lognormal(link = 'log')
#>  
#>             coef.est coef.se
#> (Intercept)     3.51    0.07
#> 
#> Dispersion parameter: 1.41
#> ML criterion at convergence: 2350.450
#> 
#> See ?tidy.sdmTMB to extract these values as a data frame.

p <- exp(predict(fit)$est[1])
p
#> [1] 33.28309

# what we have is the geometric mean:
exp(mean(log(d$density)))
#> [1] 33.28309

# but this is the mean of the data:
mean(d$density)
#> [1] 88.52506

# to get back to the mean_adjust = TRUE prediction, you'd have to do this:
phi <- exp(get_pars(fit)$ln_phi)
exp(log(p) + phi^2/2)
#> [1] 90.38351


# Now try with a linear regression on log(y):
fit <- sdmTMB(
  log(density) ~ 1,
  data = d, spatial = "off",
  family = gaussian()
)
fit
#> Model fit by ML ['sdmTMB']
#> Formula: log(density) ~ 1
#> Mesh: NULL (isotropic covariance)
#> Data: d
#> Family: gaussian(link = 'identity')
#>  
#>             coef.est coef.se
#> (Intercept)     3.51    0.07
#> 
#> Dispersion parameter: 1.41
#> ML criterion at convergence: 787.198
#> 
#> See ?tidy.sdmTMB to extract these values as a data frame.
p <- exp(predict(fit)$est[1])
p
#> [1] 33.28309

# back to the geometric mean

# and so the lognormal behaves similarly to standard GLM families by
# 'baking' in the adjustment

fit <- sdmTMB(
  density ~ 1,
  data = d, spatial = "off",
  family = Gamma(link = "log")
)
fit
#> Model fit by ML ['sdmTMB']
#> Formula: density ~ 1
#> Mesh: NULL (isotropic covariance)
#> Data: d
#> Family: Gamma(link = 'log')
#>  
#>             coef.est coef.se
#> (Intercept)     4.48    0.06
#> 
#> Dispersion parameter: 0.63
#> ML criterion at convergence: 2406.626
#> 
#> See ?tidy.sdmTMB to extract these values as a data frame.

p <- exp(predict(fit)$est[1])
p
#> [1] 88.52506
mean(d$density)
#> [1] 88.52506

# So with the adjustment, you're modelling the mean of the response (only exact
# if the data are exactly lognormal and the adjustment is perfect). This is in
# the same spirit as how we expect other GLM families to behave (and therefore
# what makes sense to me).
# Without the adjustment, you're modelling the geometric mean - same as log(y) ~
# 1, family = gaussian(), which requires adding on the sigma^2/2 if you want to
# go from geometric mean to mean.
```

<sup>Created on 2024-09-26 with [reprex v2.1.1](https://reprex.tidyverse.org)</sup>
  
