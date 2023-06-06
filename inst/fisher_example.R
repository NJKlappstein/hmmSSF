
library(ggplot2)
library(amt)
library(terra)
library(hmmSSF)
set.seed(1)

##################
## Prepare data ##
##################
tracks <- data.frame(ID = amt_fisher$id,
                     x = amt_fisher$x_,
                     y = amt_fisher$y_,
                     time = amt_fisher$t_)
tracks <- subset(tracks, ID == unique(ID)[1])
data <- random_locs(obs = tracks, n_controls = 20, distr = "gamma")

elev <- unwrap(amt_fisher_covar$elevation)
data$elev <- extract(elev, data[, c("x", "y")])$elevation

data$step <- data$step / 1000
data$x <- data$x / 1000
data$y <- data$y / 1000

######################
## HMM-SSF analysis ##
######################
f <- ~ step + cos(angle) + elev
par0 <- list(ssf_par = matrix(c(-10, -2,
                                0.5, 2,
                                0, 0),
                              ncol = 2, byrow = TRUE),
             tpm_par = matrix(c(-2, -2),
                              nrow = 1))
mod <- fitHMMSSF(ssf_formula = f, n_states = 2, data = data,
                 par0 = par0, optim_opts = list(trace = 1))

mod

################
## Plot model ##
################
data$state <- rep(factor(viterbi_decoding(mod)), each = 21)
ggplot(subset(data, obs == 1), aes(x, y, col = state, group = NA)) +
  geom_point() +
  geom_path() +
  coord_equal()
