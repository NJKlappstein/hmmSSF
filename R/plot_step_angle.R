#' Plot histograms with estimated state densities
#'
#' @param mod Fitted model object, as returned by \code{\link{fitHMMSSF}}
#' @param pal color palette to use (default = NULL)
#' @param as_list whether or not to return plots as a list
#' (default = FALSE to simply plot)
#'
#' @export

plot_step_angle <- function(mod,
                            pal = NULL,
                            as_list = FALSE) {

  states <- viterbi_decoding(mod)
  n_states <- mod$args$n_states
  obs <- subset(mod$args$data, obs == 1 & !is.na(step))

  # get weights for each state
  weights <- c(length(which(states == 1)) / nrow(obs),
               length(which(states == 2)) / nrow(obs))

  # set colour palette
  if(is.null(pal)) {
    pal <- wes_palette("Chevalier1")
    pal <- c("#3A8248", pal[2], pal[4])
  }

  #########################
  #### plot step dist ####
  #######################

  # get step parameters from ssf_par
  step_par <- beta_to_mean(fit, n_states)
  gamma_par <- get_shape_scale(fit, n_states)

  # get grid of step lengths
  max_step <- max(obs$step, na.rm = T)
  step_grid <-  seq(0, max_step, by = 0.01)
  step_lines <- NULL
  for(i in 1:n_states) {

    dens_line <- dgamma(step_grid,
                        shape = gamma_par$shape[i],
                        scale = gamma_par$scale[i]) * weights[i]
    dens_df <- data.frame(step = step_grid,
                          density = dens_line,
                          state = i)

    step_lines <-  rbind(step_lines, dens_df)
  }

  max_x <- quantile(obs$step, 0.99)
  step_obs <- subset(obs, step <= as.numeric(max_x))
  max_dens <- max(hist(step_obs$step,
                       breaks = seq(0, max_x, length.out = 21),
                       plot = FALSE)$density)
  max_y <- max_dens + max_dens * 0.25



  p_step <- ggplot(step_obs, aes(step)) +
    geom_histogram(aes(y = ..density..),
                   fill = "grey80",
                   color = 'grey15',
                   boundary = 0,
                   bins = 20) +
    geom_line(data = step_lines, aes(x = step, y = density, color = as.factor(state)),
              size = 1, alpha = 0.8) +
    scale_color_manual("state", values = pal[1:n_states]) +
    scale_y_continuous(limits = c(0, max_y)) +
    scale_x_continuous(limits = c(0, max_x)) +
    theme_light() +
    xlab("step length") +
    theme(legend.position = "none")

  #########################
  #### plot angl dist ####
  #######################

  # get grid of angles
  angle_par <- mod$par$ssf$estimate[which(mod$par$ssf$cov == "cos(angle)")]
  angle_grid <- seq(-pi, pi, by = 0.01)
  angle_lines <- NULL
  for(i in 1:n_states) {
    if(angle_par[i] < 0) {
      dens_line <- dvm(theta = angle_grid,
                       mu = pi,
                       kappa = -angle_par[i]) * weights[i] }
    else {
      dens_line <- dvm(theta = angle_grid,
                       mu = 0,
                       kappa = angle_par[i]) * weights[i]
    }
    dens_df <- data.frame(angle = angle_grid,
                          density = dens_line,
                          state = i)

    angle_lines <-  rbind(angle_lines, dens_df)
  }

  p_angle <- ggplot(obs, aes(angle)) +
    geom_histogram(aes(y = ..density..),
                   fill = "grey80",
                   color = 'grey15',
                   bins = 20) +
    geom_line(data = angle_lines, aes(x = angle, y = density, color = as.factor(state)),
              size = 1, alpha = 0.8) +
    scale_color_manual("state", values = pal[1:n_states]) +
    scale_x_continuous(limits = c(-pi, pi), breaks = c(-pi, -pi/2, 0, pi/2, pi),
                       labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
    theme_light() +
    xlab("turning angle") +
    theme(legend.position = "none")

  if(as_list) {
    p <- list(step = p_step, angle = p_angle)
  } else {
    p <- plot_grid(p_step, p_angle)
  }
  return(p)
}







