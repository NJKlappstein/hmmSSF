source("inst/functions/simHMMSSF.R")
source("inst/functions/cov_df.R")
source("inst/functions/simRaster.R")
library(raster)
library(ggplot2)

# simulate covariate data
rl <- 250 # set raster limit
set.seed(55)
cov1 <- simRaster(rho = 10, lim = c(-rl, rl, -rl, rl))
cov_data <- list("cov1" = cov1)

# set model formulation
ssf_formula <- ~ step + log(step) + cos(angle) + cov1
tpm_formula <- ~ cos(2 * pi * tod / 24) + sin(2 * pi * tod / 24)
n_states <- 2

# movement parameters
mu <- c(0.25, 2)
sd <- c(0.30, 1)
scale <- sd^2 / mu
shape <-  mu^2 / sd^2
beta1 <- -1 / scale
beta2 <- shape - 2
angle_beta <- c(0.25, 5)

# set habitat parameters
cov1_beta <- c(3, -1)

# set parameters
betas <- matrix(c(beta1,
                  beta2,
                  angle_beta,
                  cov1_beta),
                ncol = n_states,
                byrow = TRUE)

alphas <- matrix(c(-2.2, -2.2,
                   0.75, -2,
                   0.1, 1),
                 ncol = n_states^2 - n_states,
                 byrow = TRUE)

par <- list(betas = betas,
            alphas = alphas)


# set simulation settings
n_zeros <- 1e4
n_obs <- 2500
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
                   par = par,
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

# save(track, file = "data/track.RData")
# save(cov1, file = "data/cov1.RData")

# saveRDS(track, file = "inst/sim_data/sim_data.RData")
# saveRDS(cov1, file = "inst/sim_data/cov_raster.RData")
# saveRDS(par, file = "inst/sim_data/sim_par.RData")
#
# # plot data
# covmap <- data.frame(coordinates(cov1), val = values(cov1))
#
# ggplot(covmap, aes(x, y)) +
#   geom_raster(aes(fill = val)) +
#   coord_equal() +
#   xlab("Easting (km)") +
#   ylab("Northing (km)") +
#   scale_fill_viridis_c(option = "mako", guide = "none") +
#   xlim(c(-150, 150)) + ylim(c(-150, 150)) +
#   geom_path(aes(x, y), track, size = 0.35, color = "burlywood1", alpha = 0.6) +
#   geom_point(aes(x, y), track, size = 0.2, color = "burlywood1", alpha = 0.6) +
#   theme(legend.key.size = unit(1.4, "lines"))









