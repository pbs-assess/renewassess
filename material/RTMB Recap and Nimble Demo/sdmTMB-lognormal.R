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
exp(coef(fit))
mean(sim_dat$observed)
p <- predict(fit)
mean(p$est)
mean(exp(p$est))
mean(rpois(1e6, 1))

# lognormal with mean_adjustment = FALSE
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = lognormal(mean_adjust = FALSE),
  range = 0.5,
  phi = 0.8,
  sigma_O = 0,
  seed = 42,
  B = 1
)
fit <- sdmTMB(observed ~ 1, family = lognormal(mean_adjust = FALSE), data = sim_dat, spatial = "off")
logLik(fit)

coef(fit)
exp(get_pars(fit)$ln_phi)
mean(sim_dat$observed)
p <- predict(fit)
mean(p$est)

# these match:
mean(exp(p$est))
mean(rlnorm(1e6, meanlog = 1 - 0.8^2/2, sdlog = 0.8))

# these match:
mean(exp(p$est + 0.8^2/2))
mean(rlnorm(1e6, meanlog = 1, sdlog = 0.8))

# redo the fit with mean_adjustment = TRUE
fit2 <- sdmTMB(observed ~ 1, family = lognormal(mean_adjust = TRUE),
  data = sim_dat, spatial = "off")
logLik(fit)
logLik(fit2)

coef(fit2)
exp(get_pars(fit2)$ln_phi)
mean(sim_dat$observed)
p <- predict(fit2)
mean(p$est)

# now these match:
mean(exp(p$est))
mean(rlnorm(1e6, meanlog = 1, sdlog = 0.8))

# now these match:
mean(exp(p$est - 0.8^2/2))
mean(rlnorm(1e6, meanlog = 1 - 0.8^2/2, sdlog = 0.8))

# mean_adjust = TRUE in sim and fit:

# lognormal with mean_adjustment = FALSE
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = lognormal(mean_adjust = TRUE),
  range = 0.5,
  phi = 0.8,
  sigma_O = 0,
  seed = 42,
  B = 1
)
fit3 <- sdmTMB(observed ~ 1, family = lognormal(mean_adjust = TRUE), data = sim_dat, spatial = "off")
logLik(fit3)

coef(fit3)
exp(get_pars(fit3)$ln_phi)
mean(sim_dat$observed)
p <- predict(fit3)
mean(p$est)

# these match:
mean(exp(p$est))
mean(rlnorm(1e6, meanlog = 1 - 0.8^2/2, sdlog = 0.8))

# these match:
mean(exp(p$est + 0.8^2/2))
mean(rlnorm(1e6, meanlog = 1, sdlog = 0.8))

