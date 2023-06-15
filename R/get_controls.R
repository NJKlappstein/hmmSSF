#' Simulate random locations for model fitting
#'
#' @param obs Data frame of observations (ID, x, y, time)
#' @param n_controls Number of control steps ("controls") per observation
#' @param distr Name(s) of distribution(s) used to generate control steps.
#' Spatially uniform ("uniform") points are the default. Otherwise, the distributions
#' should be entered as a vector of the step length and turning angle distributions
#' (in that order). If only one distribution is entered, this is used as the
#' step length distribution. Current options for step length: "gamma", "exponential".
#' Current options for turning angle: "vm", "uniform", "wrpcauchy".
#'
#' @export

get_controls <- function(obs,
                        n_controls,
                        distr = "uniform")
{
  # Setup sampling distributions
  setup <- setup_samp(distr = distr, obs = obs)

  # Loop over IDs
  data_all <- NULL
  for(id in unique(obs$ID)) {
    # Subset to this ID
    sub_obs <- subset(obs, obs$ID == id)

    # Matrix of locations
    xy <- as.matrix(sub_obs[, c("x", "y")])
    n_obs <- nrow(xy)

    # Get a few movement metrics
    steps <- sqrt(rowSums((xy[-1,] - xy[-n_obs,])^2))
    bears <- atan2(diff(xy[,2]), diff(xy[,1]))
    angles <- shift_angle(diff(bears))

    # Loop over observed steps
    data_list <- lapply(3:n_obs, function(i) {
      # Simulate step lengths and turning angles for random steps
      sim_step <- setup$r_step(n = n_controls, par = setup$par_step)
      sim_angle <- shift_angle(setup$r_angle(n = n_controls,
                                             par = setup$par_angle))
      # Importance weights
      weights <- setup$d_step(x = sim_step, par = setup$par_step) *
        setup$d_angle(x = sim_angle, par = setup$par_angle)
      weight_obs <- setup$d_step(x = steps[i-1], par = setup$par_step) *
        setup$d_angle(x = angles[i-2], par = setup$par_angle)
      # Derive coordinates of random points
      sim_bear <- bears[i-2] + sim_angle
      pts <- rep(xy[i-1,], each = n_controls) +
        sim_step * cbind(cos(sim_bear), sin(sim_bear))

      data.frame(ID = id,
                 stratum = i,
                 obs = c(1, rep(0, n_controls)),
                 x = c(xy[i,1], pts[,1]),
                 y = c(xy[i,2], pts[,2]),
                 step = c(steps[i-1], sim_step),
                 angle = c(angles[i-2], sim_angle),
                 time = as.POSIXct(sub_obs$time[i],
                                   tz = attr(sub_obs$time, "tzone")),
                 w = c(weight_obs, weights))
    })

    # Add data for this ID to data frame
    data_all <- rbind(data_all, do.call(rbind, data_list))
  }

  # Save importance weights as attributes and remove column
  attr(data_all, "weights") <- data_all$w
  data_all$w <- NULL

  return(data_all)
}
