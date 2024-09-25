``` r
# devtools::install_github("pbs-assess/sdmTMB", ref = "lognormal")
library(sdmTMB)
d <- subset(pcod_2011, density > 0)
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

# but this is the mean of the data
mean(d$density)
#> [1] 88.52506

# to get there, you'd have to do this
phi <- exp(get_pars(fit)$ln_phi)
exp(log(p) + phi^2/2)
#> [1] 90.38351

# and so the lognormal behaves like other distributions by 'baking' in the 
# adjustment

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
```

<sup>Created on 2024-09-25 with [reprex v2.1.1](https://reprex.tidyverse.org)</sup>
