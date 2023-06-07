
source("inst/functions/simHMMSSF.R")
source("inst/functions/cov_df.R")
source("inst/functions/simRaster.R")
library(raster)
library(ggplot2)
set.seed(55)

# simulate covariate data
rl <- 250 # set raster limit

cov1 <- simRaster(rho = 10, lim = c(-rl, rl, -rl, rl))
cov_data <- list("cov1" = cov1)

# set model formulation
ssf_formula <- ~ step + cos(angle) + cov1
tpm_formula <- ~ cos(2 * pi * tod / 24) + sin(2 * pi * tod / 24)
n_states <- 2

# simulation parameters
ssf_par <- matrix(c(-10, -1,
                    0.5, 5,
                    1, 0),
                  ncol = n_states, byrow = TRUE)
tpm_par <- matrix(c(-2.2, -2.2,
                    0.75, -2,
                    0.1, 1),
                  ncol = n_states * (n_states - 1), byrow = TRUE)

# set simulation settings
n_zeros <- 1e4
n_obs <- 1000
rmax <- 5
n_iter <- 100
time <- data.frame(time = seq(as.POSIXct("2020-01-01 0:00", tz = "UTC"),
                              by = "hour", length.out = n_obs))
time$tod <- lubridate::hour(time$time)

# simulate track
track <- simHMMSSF(ssf_formula = ssf_formula,
                   tpm_formula = tpm_formula,
                   ssf_cov = cov_data,
                   tpm_cov = time,
                   par = list(betas = ssf_par, alphas = tpm_par),
                   n_states = n_states,
                   n_obs = n_obs,
                   n_zeros = n_zeros,
                   y1 = c(0, 0),
                   rmax = rmax)

# process data
track$ID <- 1
track$time <- time$time
names(track) <- c("x", "y", "state", "ID", "time")

track <- data.frame(ID = 1,
                    x = track$x,
                    y = track$y,
                    time = track$time,
                    state = track$state)

# plot data
plot(cov1, xlim = range(track$x) + c(-5, 5), ylim = range(track$y) + c(-5, 5))
points(track$x, track$y, type = "o", pch = 20, cex = 0.3)

save(cov1, track, file = "inst/sim_example/sim_data.RData")
