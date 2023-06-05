library(hmmSSF)
library(terra)
library(ggplot2)

# load data
track <- readRDS("inst/sim_data/sim_data.RData")
cov_data <- readRDS("inst/sim_data/cov_raster.RData")

track <- track[c(1:500),]


data <- sim_controls2(obs = track,
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
