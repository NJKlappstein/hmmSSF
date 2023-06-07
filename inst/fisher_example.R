
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
tracks <- subset(tracks, ID == unique(ID)[4])
n_random <- 50
data <- random_locs(obs = tracks, n_controls = n_random, distr = "gamma")

elev <- unwrap(amt_fisher_covar$elevation)
data$elev <- extract(elev, data[, c("x", "y")])$elevation

data$step <- data$step / 1000
data$x <- data$x / 1000
data$y <- data$y / 1000

ggplot(subset(data, obs == 1), aes(x, y, group = ID)) +
  geom_point(size = 0.3) +
  geom_path() +
  coord_equal()

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

mod

################
## Plot track ##
################
data$state <- rep(factor(viterbi_decoding(mod)), each = n_random + 1)
ggplot(subset(data, obs == 1), aes(x, y, col = state, group = NA)) +
  geom_point() +
  geom_path() +
  coord_equal()

#################################
## Plot movement distributions ##
#################################
step_grid <- seq(min(data$step, na.rm = TRUE),
                 max(data$step, na.rm = TRUE), length = 100)
new_data <- data.frame(step = step_grid, angle = pi/2, elev = 0)
ssf_MM <- model.matrix(f, new_data)[,-1]
par1 <- mod$par$ssf[which(mod$par$ssf$state == 1), "estimate"]
par2 <- mod$par$ssf[which(mod$par$ssf$state == 2), "estimate"]
new_data$pred1 <- exp(ssf_MM %*% par1)
new_data$pred2 <- exp(ssf_MM %*% par2)
new_data$pred1 <- new_data$pred1 / sum(new_data$pred1)
new_data$pred2 <- new_data$pred2 / sum(new_data$pred2)
ggplot(new_data, aes(step, pred1)) +
  geom_line(col = "firebrick") +
  geom_line(aes(step, pred2), col = "royalblue")

angle_grid <- seq(-pi, pi, length = 100)
new_data <- data.frame(step = 1, angle = angle_grid, elev = 0)
ssf_MM <- model.matrix(f, new_data)[,-1]
new_data$pred1 <- exp(ssf_MM %*% par1)
new_data$pred2 <- exp(ssf_MM %*% par2)
new_data$pred1 <- new_data$pred1 / sum(new_data$pred1)
new_data$pred2 <- new_data$pred2 / sum(new_data$pred2)
ggplot(new_data, aes(angle, pred1)) +
  geom_line(col = "firebrick") +
  geom_line(aes(angle, pred2), col = "royalblue")
