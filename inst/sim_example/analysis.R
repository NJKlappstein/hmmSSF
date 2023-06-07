
library(hmmSSF)

# load data
load("inst/sim_example/sim_data.RData")

# add random locations
nc <- 50
data <- random_locs(obs = track, n_controls = nc, dist = "gamma")
data$cov1 <- extract(cov1, as.matrix(data[, c("x", "y")]))
data$tod <- rep(lubridate::hour(track$time)[-(1:2)], each = nc + 1)

# set model formulation
ssf_formula <- ~ step + cos(angle) + cov1
tpm_formula <- ~ cos(2 * pi * tod / 24) + sin(2 * pi * tod / 24)
n_states <- 2
ssf_par0 <- matrix(c(-8, -0.5,
                     0.2, 5,
                     0, 0),
                   ncol = 2, byrow = TRUE)

# fit model
mod <- fitHMMSSF(ssf_formula = ssf_formula,
                 tpm_formula = tpm_formula,
                 n_states = n_states,
                 data = data,
                 ssf_par0 = ssf_par0,
                 optim_opts = list(trace = 1, maxit = 1e4))
