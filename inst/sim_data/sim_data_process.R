library(hmmSSF)
library(terra)
library(ggplot2)

# load data
track <- readRDS("inst/sim_data/sim_data.RData")
cov_data <- readRDS("inst/sim_data/cov_raster.RData")
par <- readRDS("inst/sim_data/sim_par.RData")
#track <- track[c(1:200),]

# add random locations
data <- random_locs(obs = track,
                    n_controls = 20,
                    dist = c("exp", "vm"))
data$cov1 <- extract(cov_data, as.matrix(data[, c("x", "y")]))
data$tod <- rep(lubridate::hour(track$time)[-(1:2)], each = 21)

# set model formulation
ssf_formula <- ~ step + log(step) + cos(angle) + cov1
tpm_formula <- ~ cos(2 * pi * tod / 24) + sin(2 * pi * tod / 24)
n_states <- 2
ssf_par0 <- matrix(c(-2, -2,
                     -1, 2,
                     0.2, 5,
                     3, -1),
                   ncol = 2, byrow = TRUE)
tpm_par0 <- matrix(c(-2, -2,
                     1, -2,
                     0.1, 1),
                   ncol = 2, byrow = TRUE)

# fit model
mod <- fitHMMSSF(ssf_formula = ssf_formula,
                 tpm_formula = tpm_formula,
                 n_states = n_states,
                 data = data,
                 ssf_par0 = ssf_par0,
                 tpm_par0 = tpm_par0,
                 optim_opts = list(trace = 1))


mod
confint(mod)
mod$par


states <- viterbi_decoding(mod)
sp <- local_decoding(mod)

new_data <- data.frame(tod = 0:24)
foo <- predict_tpm(mod, new_data = new_data)
bar <- predict_delta(mod, new_data = new_data)



