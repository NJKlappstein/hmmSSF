
library(ggplot2)
theme_set(theme_bw())
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
tracks <- subset(tracks, ID %in% unique(ID)[c(1, 3)])
n_random <- 50
data <- random_locs(obs = tracks, n_controls = n_random, distr = "gamma")

elev <- unwrap(amt_fisher_covar$elevation)
data$elev <- extract(elev, data[, c("x", "y")])$elevation

data$step <- data$step / 1000
data$x <- data$x / 1000
data$y <- data$y / 1000

# ggplot(subset(data, obs == 1), aes(x, y, group = ID)) +
#   geom_point(size = 0.3) +
#   geom_path() +
#   coord_equal()

######################
## HMM-SSF analysis ##
######################
f <- ~ step + log(step) + cos(angle) + elev
ssf_par0 <- matrix(c(-10, -2,
                     0, 0,
                     0.5, 2,
                     0, 0),
                   ncol = 2, byrow = TRUE)
mod <- fitHMMSSF(ssf_formula = f, n_states = 2, data = data,
                 ssf_par0 = ssf_par0, optim_opts = list(trace = 1))

################
## Plot track ##
################
data$state <- rep(factor(viterbi_decoding(mod)), each = n_random + 1)
ggplot(subset(data, obs == 1), aes(x, y, col = state, group = ID)) +
  facet_wrap("ID") +
  geom_point() +
  geom_path() +
  coord_equal()

#################################
## Plot movement distributions ##
#################################
plot_ssf(mod, "step")
plot_ssf(mod, "angle")
plot_ssf(mod, "elev")
