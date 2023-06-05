library(hmmSSF)

# load data
data <- readRDS("inst/zebra/zebra_processed.RData")
initial_par <- readRDS("inst/zebra/initial_par.RData")
initial_par <- initial_par[[3]]

# model settings
n_states <- 2
ssf_formula <- ~ step + log(step) + cos(angle) + veg
tpm_formula <- ~ cos(2*pi*tod/24) + sin(2*pi*tod/24)

# fit model
fit <- fitHMMSSF(ssf_formula = ssf_formula,
                 tpm_formula = tpm_formula,
                 data = data,
                 par0 = initial_par,
                 n_states = n_states,
                 dist = "gamma",
                 optim_opts = list(trace = 1,
                                   maxit = 1e4))
fit$betas
