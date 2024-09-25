library(sdmTMB)
set.seed(1)
predictor_dat <- data.frame(
  X = runif(20000), Y = runif(20000)
)
mesh <- make_mesh(predictor_dat, c("X", "Y"), cutoff = 0.1)

# Poisson:
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = poisson(),
  range = 0.5,
  phi = 0.8,
  sigma_O = 0,
  seed = 42,
  B = 1
)
fit <- sdmTMB(observed ~ 1, family = poisson(), data = sim_dat, spatial = "off")

coef(fit)
mean(sim_dat$observed)
p <- predict(fit)
mean(p$est)
mean(exp(p$est))
mean(rpois(1e6, 1))

# lognormal:
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = lognormal(),
  range = 0.5,
  phi = 0.8,
  sigma_O = 0,
  seed = 42,
  B = 1
)
fit <- sdmTMB(observed ~ 1, family = lognormal(), data = sim_dat, spatial = "off")

coef(fit)
exp(get_pars(fit)$ln_phi)
mean(sim_dat$observed)
p <- predict(fit)
mean(p$est)
mean(exp(p$est))
mean(rlnorm(1e6, meanlog = 1, sdlog = 0.8))
mean(rlnorm(1e6, meanlog = 1 - 0.8^2/2, sdlog = 0.8))

# so the simulation/estimation is set up such that
# the mean of the linear predictor (in log space here)
# is 1 (the intercept value)
mean(log(rlnorm(1e6, meanlog = 1, sdlog = 0.8)))

# Generalized gamma:
# 'gengamma Q' is set at 0.5 here by default
# (a bit heavier tailed than the Gamma (phi = Q)
# but not as much as the lognormal (Q = 0))
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = gengamma(),
  range = 0.5,
  phi = 0.8,
  sigma_O = 0,
  seed = 42,
  B = 1
)
fit <- sdmTMB(observed ~ 1, family = gengamma(), data = sim_dat, spatial = "off")

coef(fit)
exp(get_pars(fit)$ln_phi)
get_pars(fit)$gengamma_Q
mean(sim_dat$observed)
p <- predict(fit)
mean(p$est)
mean(exp(p$est))
