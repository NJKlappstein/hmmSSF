library(hmmSSF)
library(raster)
library(ggplot2)

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




###########################
###### plots to compare ###
###########################
CIs <- confint(fit)

# ssf plots
ssf_df <- data.frame(state = rep(c(1, 2), each = 3),
                     covariate = rep(c("bushed \n grassland",
                                       "bushland",
                                       "woodland"),
                                     times = 2),
                     estimate = c(fit$par$ssf[c(4:6, 10:12)]),
                     upper = c(CIs$ssf$upp[c(4:6, 10:12)]),
                     lower = c(CIs$ssf$low[c(4:6, 10:12)]))

pd <- position_dodge(0.25)
ggplot(ssf_df, aes(y = estimate,
                         x = covariate,
                         group = state,
                         color = as.factor(state))) +
  geom_point(position = pd, size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = pd,
                size = 0.5,
                width = 0.15) +
  xlab("habitat type") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_light(13)


#tpm plots
new_data <- data.frame(tod = seq(0, 24, length = 100))
delta <- predict_delta(fit, new_data)

df <- as.data.frame.table(delta$mle)
colnames(df) <- c("tod", "state", "value")
df$tod <- new_data$tod
df$low <- as.data.frame.table(delta$lower)[,3]
df$upp <- as.data.frame.table(delta$upper)[,3]
df$state <- ifelse(df$state == 1, "encamped", "exploratory")
pal <- c("#0072B2", "#D55E00")


ggplot(df, aes(tod, value, group = state, col = factor(state), fill = factor(state))) +
  geom_ribbon(aes(ymin = low, ymax = upp), col = NA, alpha = 0.3) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "time of day", y = "stationary probability") +
  scale_color_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), labels = c("0", "0.5", "1")) +
  # theme(legend.position = c(0.1, 0.8))
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10.5))




