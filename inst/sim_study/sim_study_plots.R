library(ggplot2)
library(dplyr)
library(CircStats)
library(cowplot)

sim_par <- readRDS("inst/sim_study/sim_par.RData")
betas_out <- read.csv("inst/sim_study/ssf_est.csv")
alphas_out <- read.csv("inst/sim_study/tpm_est.csv")
vit <- read.csv("inst/sim_study/vit_est.csv")

# get true shape and scale parameters
sim_shape <- sim_par$betas[2,] + 2
sim_scale <- - 1 / sim_par$betas[1,]
sim_conc <- sim_par$betas[3,]

# set plot loops
step <- seq(0.05, 6, by = 0.05)
angle <- seq(-pi, pi, by = 0.01)
step_df <- NULL
angle_df <- NULL
states <- c(1:2)

for(i in 1:50) {
  for(k in 1:2) {
    par <- subset(betas_out,
                  track == i &
                    state == states[k])

    # get shape and scale parameters and step density
    shape <- par[which(par$cov == "log(step)"),]$estimate + 2
    scale <- - 1 / par[which(par$cov == "step"),]$estimate
    step_dens <- dgamma(step, shape = shape, scale = scale)

    # get angle density
    angle_dens <- dvm(angle,
                      mu = 0,
                      kappa = par[which(par$cov == "cos(angle)"),]$estimate)

    #create df
    step_foo <- data.frame(track = i,
                           state = k,
                           step = step,
                           density = step_dens)
    angl_foo <- data.frame(track = i,
                           state = k,
                           angle = angle,
                           density = angle_dens)

    # bind with main df
    step_df <- rbind(step_df, step_foo)
    angle_df <- rbind(angle_df, angl_foo)

  }
}



true_step <- data.frame(state = rep(c(1:2), each = length(step)),
                        step = rep(step, times = 2),
                        density = c(dgamma(step,
                                           shape = sim_shape[1],
                                           scale = sim_scale[1]),
                                    dgamma(step,
                                           shape = sim_shape[2],
                                           scale = sim_scale[2])))

pal <- c("#02808f", "#FDD262", "#f09537")
p_step <- ggplot(step_df, aes(x = step, y = density, group = interaction(track, state),
                    color  = as.factor(state))) +
  geom_line(alpha = 0.3) +
  geom_line(data = true_step,
            aes(x = step,
                y = density,
                group = as.factor(state)),
            size = 0.5, linetype = "dashed", color = "black") +
  scale_color_manual("state", values = pal[1:2]) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  ylab("density \n") +
  theme(legend.position = "none",
        legend.title = element_text(size=17),
        legend.text = element_text(size = 13),
        legend.background = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=16),
        panel.grid.major = element_line(color="grey90", size = 0.1),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        legend.key=element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha=1)))


true_angl <- data.frame(state = rep(c(1:2), each = length(angle)),
                        angle = rep(angle, times = 2),
                        density = c(dvm(angle, mu = 0, kappa = sim_conc[1]),
                                    dvm(angle, mu = 0, kappa = sim_conc[2])))

p_angle <- ggplot(angle_df, aes(x = angle, y = density, group = interaction(track, state),
                     color  = as.factor(state))) +
  geom_line(alpha = 0.3) +
  geom_line(data = true_angl,
            aes(x = angle,
                y = density,
                group = as.factor(state)),
            size = 0.5, linetype = "dashed", color = "black") +
  scale_color_manual("state", values = pal[1:2]) +
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size = 13),
        legend.background = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title=element_text(size=16),
        panel.grid.major = element_line(color="grey90", size = 0.1),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        legend.key=element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha=1)))

#################################
#### habitat parameters ########
################################
hline_b <- data.frame(state = rep(c(1, 2), each = 2),
                      cov = rep(c("covariate 1", "covariate 2"), times = 2),
                      beta = c(sim_par$betas[4,], sim_par$betas[5,]))

betas_cov <- subset(betas_out, cov == "cov1" | cov == "cov2")
betas_cov$cov <- ifelse(betas_cov$cov == "cov1", "covariate 1", "covariate 2")


p_cov <- ggplot(betas_cov, aes(x = as.factor(state), y = estimate, fill = as.factor(state))) +
  geom_boxplot() +
  geom_point(data = hline_b, aes(x = state, y = beta), size = 3, shape = 4) +
  scale_fill_manual("state", values = pal[1:2]) +
  facet_wrap(facets = "cov") +
  xlab("state") +
  theme(legend.position = "none",
        legend.background = element_blank(),
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        strip.text = element_text(size = 17),
        strip.background = element_blank(),
        axis.title = element_text(size=15),
        panel.grid.major = element_line(color="grey90", size = 0.1),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        plot.title = element_text(size = 18))

plot_grid(p_step, p_angle, p_cov)











