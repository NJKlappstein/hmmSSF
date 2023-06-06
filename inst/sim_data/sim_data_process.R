library(hmmSSF)
library(terra)
library(ggplot2)

# load data
track <- readRDS("inst/sim_data/sim_data.RData")
cov_data <- readRDS("inst/sim_data/cov_raster.RData")

track <- track[c(1:200),]


data <- random_locs(obs = track,
                    n_controls = 20,
                    dist = c("exp", "vm"))


hist(subset(data, obs ==0)$step)
hist(subset(data, obs ==0)$angle)

ggplot(subset(data, stratum %in% 3:4),
       aes(x,
           y,
           color = factor(stratum),
           group = factor(obs),
           shape = factor(obs),
           size = factor(obs))) +
  geom_point(alpha = 0.5, size = 0.5) +
  coord_equal()

data$cov1 <- extract(cov_data, as.matrix(data[, c("x", "y")]))
data$tod <- rep(lubridate::hour(track$time)[-(1:2)], each = 21)

# set model formulation
ssf_formula <- ~ step + log(step) + cos(angle) + cov1
tpm_formula <- ~ cos(2 * pi * tod / 24) + sin(2 * pi * tod / 24)
n_states <- 2
par0 <- list(betas = matrix(c(-2, -2,
                              -1, 2,
                              0.2, 5,
                              3, -1),
                            ncol = 2, byrow = TRUE),
             alphas = matrix(c(-2, -2,
                               1, -2,
                               0.1, 1),
                             ncol = 2, byrow = TRUE))

mod <- fitHMMSSF(ssf_formula = ssf_formula, tpm_formula = tpm_formula,
                 data = data, par0 = par0, n_states = n_states)


mod
names(mod)
mod$fit
mod$par_CI
mod$args
coefficients(mod)$ssf_par
lower(mod)






