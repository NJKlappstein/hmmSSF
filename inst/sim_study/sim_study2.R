library(raster)
library(tictoc)
library(hmmSSF)
library(pbmcapply)
source("inst/functions/simHMMSSF.R")
source("inst/functions/cov_df.R")
source("inst/functions/simRaster.R")

######################################
## implementation simulation study ##
#####################################

# simulate covariate data
rl <- 2000 # set raster limit
set.seed(435363)
cov1 <- simRaster(rho = 3, lim = c(-rl, rl, -rl, rl))
cov2 <- simRaster(rho = 3, lim = c(-rl, rl, -rl, rl))
cov_data <- list("cov1" = cov1, "cov2" = cov2)

# set model formulation
ssf_formula <- ~ step + cos(angle) + cov1
tpm_formula <- ~ cos(2 * pi * tod / 24) + sin(2 * pi * tod / 24)
n_states <- 2

# set parameters
betas <- matrix(c(-10, -1,
                  -1, 5,
                  3, -3),
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
n_obs <- 1000
rmax <- 5
n_iter <- 50
time <- data.frame(time = seq(as.POSIXct("2020-01-01 0:00", tz = "UTC"),
                              by = "hour", length.out = n_obs))
time$tod <- lubridate::hour(time$time)

out <- list()
out <- pbmclapply(1:n_iter, function(i) {
  # get initial location for track
  y1 <-  c(runif(1, c(-rl*0.5, rl*0.5)),
           runif(1, c(-rl*0.5, rl*0.5)))

  # simulate track
  track <- simHMMSSF(ssf_formula = ssf_formula,
                     tpm_formula = tpm_formula,
                     ssf_cov = cov_data,
                     tpm_cov = time,
                     par = par,
                     n_states = n_states,
                     n_obs = n_obs,
                     n_zeros = n_zeros,
                     y1 = y1,
                     rmax = rmax,
                     print = FALSE)

  # process data
  track$ID <- 1
  track$time <- time$time
  names(track) <- c("x", "y", "state", "ID", "time")

  # set starting values for ssf
  ssf_par0 <- matrix(c(-8, -2,
                       0, 2,
                       2, -2),
                     ncol = n_states,
                     byrow = TRUE)

  # set tpm starting values
  tpm_par0 <- matrix(c(-2, -2,
                       0, 0,
                       0, 0),
                     ncol = n_states^2 - n_states,
                     byrow = TRUE)


  # generate controls
  data <- random_locs(obs = track,
                      n_controls = 200,
                      distr = c("gamma", "vm"))

  # get covariates
  data$cov1 <- extract(cov_data$cov1, as.matrix(data[, c("x", "y")]))
  data$cov2 <- extract(cov_data$cov2, as.matrix(data[, c("x", "y")]))
  data$tod <- lubridate::hour(data$time)

  # fit HMM-SSF
  fit <- fitHMMSSF(ssf_formula = ssf_formula,
                   tpm_formula = tpm_formula,
                   n_states = n_states,
                   data = data,
                   ssf_par0 = ssf_par0,
                   tpm_par0 = tpm_par0,
                   optim_opts = list(maxit = 1e4))

  # put all data into data frame
  return(fit$par)
}, mc.cores = 4)

par(mfrow = c(3, 2))
for(k in 1:3) {
  for(l in 1:2) {
    hist(sapply(out, function(m) m$ssf[k,l]), main = "")
    abline(v = betas[k,l], col = "firebrick", lwd = 2)
  }
}

par(mfrow = c(2, 2))
for(k in 1:2) {
  for(l in 1:2) {
    hist(sapply(out, function(m) m$tpm[k,l]), main = "")
    abline(v = alphas[k,l], col = "firebrick", lwd = 2)
  }
}
