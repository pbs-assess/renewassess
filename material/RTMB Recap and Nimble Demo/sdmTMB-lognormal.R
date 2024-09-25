library(sdmTMB)
set.seed(1)
predictor_dat <- data.frame(
  X = runif(20000), Y = runif(20000)
)
mesh <- make_mesh(predictor_dat, c("X", "Y"), cutoff = 0.1)
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
head(sim_dat)

fit <- sdmTMB(observed ~ 1, family = lognormal(), data = sim_dat, spatial = "off")

coef(fit)
exp(get_pars(fit)$ln_phi)
