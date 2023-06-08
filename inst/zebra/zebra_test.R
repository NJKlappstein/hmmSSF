library(hmmSSF)
library(raster)

# load data
data <- readRDS("inst/zebra/zebra_processed.RData")
track <- subset(data, obs == 1)
hb <- raster("inst/zebra/vegetation2.grd")

# generate random locations using new function
data <- random_locs(obs = track, n_controls = 25, distr = "gamma")

# add needed covariates
notNA <-  which(!is.na(data$x))
data$veg[notNA] <- extract(hb, data[notNA, c("x", "y")], method = "simple")
data$veg <- factor(data$veg)
levels(data$veg) <- c("grassland", "bushed grassland", "bushland", "woodland")
data$tod <- as.numeric(format(data$time, "%H")) +
  as.numeric(format(data$time, "%M"))/60
data$step[which(data$step == 0)] <- data$step[which(data$step == 0)] +
  runif(length(data$step[which(data$step == 0)]), 0, min(data$step[which(data$step > 0)]))

# model settings
initial_par <- readRDS("inst/zebra/initial_par.RData")
initial_par <- initial_par[[3]]
n_states <- 2
ssf_formula <- ~ step + log(step) + cos(angle) + veg
tpm_formula <- ~ cos(2*pi*tod/24) + sin(2*pi*tod/24)

# fit model
fit <- fitHMMSSF(ssf_formula = ssf_formula,
                 tpm_formula = tpm_formula,
                 data = data,
                 ssf_par0 = initial_par$betas,
                 tpm_par0 = initial_par$alphas,
                 n_states = n_states,
                 optim_opts = list(trace = 1,
                                   maxit = 1e4))
fit
