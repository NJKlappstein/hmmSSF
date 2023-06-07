library(tictoc)
source("inst/functions/simHMMSSF.R")
source("inst/functions/cov_df.R")
source("inst/functions/simRaster.R")

######################################
## implementation simulation study ##
#####################################

# simulate covariate data
rl <- 2000 # set raster limit
set.seed(50)
cov1 <- simRaster(rho = 5, lim = c(-rl, rl, -rl, rl))
cov2 <- simRaster(rho = 25, lim = c(-rl, rl, -rl, rl))
cov_data <- list("cov1" = cov1, "cov2" = cov2)

# set model formulation
ssf_formula <- ~ step + log(step) + cos(angle) + cov1 + cov2
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
cov2_beta <- c(-1, 3)

# set parameters
betas <- matrix(c(beta1,
                  beta2,
                  angle_beta,
                  cov1_beta,
                  cov2_beta),
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
n_obs <- 750
rmax <- 5
n_iter <- 50
time <- data.frame(time = seq(as.POSIXct("2020-01-01 0:00", tz = "UTC"),
                              by = "hour", length.out = n_obs))
time$tod <- as.numeric(substr(time$time, 12, 13))

#initialise arrays for parameters
ssf  <- NULL
tpm <- NULL
vit <- NULL

tic("sim time")
for(i in 1:n_iter) {
  cat("\rIteration", i, "of", n_iter)

  # set seed so I can recreate each track if needed
  set.seed(i + 501)

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
  step_beta0    <- c(runif(1, -6, -2), runif(1, -3, -1))
  logstep_beta0 <- c(runif(1, -1, 1), runif(1, 0, 2))
  angle_beta0   <- c(runif(1, 0, 1), runif(1, 1, 5))
  cov1_beta0    <- c(runif(1, 1, 5), runif(1, 0, 4))
  cov2_beta0    <- c(runif(1, 0, 4), runif(1, 1, 5))

  ssf_par0 <- matrix(c(step_beta0,
                       logstep_beta0,
                       angle_beta0,
                       cov1_beta0,
                       cov2_beta0),
                     ncol = n_states,
                     byrow = TRUE)

  # set tpm starting values
  tpm_par0 <- matrix(c(-2, -2,
                       runif(1, 0, 1.5), runif(1, -3, -1),
                       runif(1, 0, 0.5), runif(1, 0, 2)),
                     ncol = n_states^2 - n_states,
                     byrow = TRUE)


  # generate controls
  data <- random_locs(obs = track,
                      n_controls = 35,
                      distr = "gamma")

  # get covariates
  data$cov1 <- extract(cov_data$cov1, as.matrix(data[, c("x", "y")]))
  data$cov2 <- extract(cov_data$cov2, as.matrix(data[, c("x", "y")]))
  data$tod <- as.numeric(substr(data$time, 12, 13))

  # fit HMM-SSF
  fit <- fitHMMSSF(ssf_formula = ssf_formula,
                   tpm_formula = tpm_formula,
                   n_states = n_states,
                   data = data,
                   ssf_par0 = ssf_par0,
                   tpm_par0 = tpm_par0)

  # put all data into data frame
  ssf_est <- data.frame(track = i,
                        state = rep(c(1, 2), each = 5),
                        cov = rownames(fit$par$ssf),
                        estimate = as.vector(fit$par$ssf))

  tpm_est <- data.frame(track = i,
                        transition = rep(c("1-2", "2-1"), each = 3),
                        cov = rownames(fit$par$tpm),
                        estimate = as.vector(fit$par$tpm))

  #bind with other estimates
  ssf <- rbind(ssf, ssf_est)
  tpm <- rbind(tpm, tpm_est)

  # decoding and convergence
  decode <- c(NA, NA, viterbi_decoding(fit))
  decode_perc <- length(which(decode == track$state)) / n_obs

  decode <- data.frame(track = i,
                       convergence = fit$fit$convergence,
                       decoded = decode_perc)

  vit <- rbind(vit, decode)
}
toc()

write.csv(ssf, "results/simulations/csv_out/ssf_est.csv", row.names = F)
write.csv(tpm, "results/simulations/csv_out/tpm_est.csv", row.names = F)
write.csv(vit, "results/simulations/csv_out/vit_est.csv", row.names = F)


betas_out <- read.csv("results/simulations/csv_out/TPM_sim_betas.csv")
alphas_out <- read.csv("results/simulations/csv_out/TPM_sim_alphas.csv")
decoding <- read.csv("results/simulations/csv_out/TPM_decode.csv")


betas_out$se <- (betas_out$upper - betas_out$estimate) / 1.96
betas_out$cov <- factor(betas_out$cov,
                        levels = c("step", "log(step)", "cos(angle)", "cov1", "cov2"))
b_summary <- betas_out %>%
  group_by(state, cov) %>%
  summarise(mean_se = mean(se),
            est_sd = sd(estimate),
            median = median(estimate),
            mean = mean(estimate))
b_summary$truth <- as.vector(par$betas)
b_summary$bias <- (b_summary$truth - b_summary$mean) / b_summary$truth
b_summary


mean(decoding$decoded)

alphas_out$se <- (alphas_out$upper - alphas_out$estimate) / 1.96
a_summary <- alphas_out %>%
  group_by(transition, cov) %>%
  summarise(mean_se = mean(se),
            est_sd = sd(estimate),
            median = median(estimate),
            mean = mean(estimate))
a_summary$truth <- as.vector(par$alphas)
a_summary$bias <- (a_summary$truth - a_summary$mean) / a_summary$truth
a_summary


alphas_out %>% group_by(transition, cov) %>%
  summarise(mean = mean(estimate), median = median(estimate))
sim_par$alphas


#############################
##### covariate plots ######
###########################
cov1_plot <- data.frame(coordinates(cov1),
                        val = values(cov1),
                        cov = "covariate 1")
cov2_plot <- data.frame(coordinates(cov2),
                        val = values(cov2),
                        cov = "covariate 2")
cov_plot <- rbind(cov1_plot, cov2_plot)
cov_plot <- subset(cov_plot, x <= 25 & x >= -25)
cov_plot <- subset(cov_plot, y <= 25 & y >= -25)

pdf(file = "writing/figures/sim_hab.pdf", width = 6, height = 3)
ggplot(cov_plot, aes(x, y)) +
  geom_raster(aes(fill = val)) +
  facet_wrap("cov") +
  coord_equal() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_scico(palette = "davos") +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))
dev.off()














