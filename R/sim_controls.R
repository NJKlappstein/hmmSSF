##' Simulate control steps
##'
##' @param obs Data frame of observations (ID, x, y, time)
##' @param n_controls Number of control steps ("controls") per observation
##' @param step_dist distribution from which to generate control locations
##' (can be "uniform" or "gamma")
##'
##' @export


sim_controls <- function(obs, n_controls, step_dist)
{
  ## Number of observations
  n_obs <- nrow(obs)

  ## Derive availability radius from observed step lengths (requires moveHMM)
  step_data <- prepData(obs[,c("ID", "x", "y")], type = "UTM")
  steps <- step_data$step[which(!is.na(step_data$step))]
  rmax <- 1.1 * max(steps)
  mean <- mean(steps, na.rm = T)
  sd <- sd(steps, na.rm = T)

  ## Prepare data frame
  n_all <- n_obs * (n_controls + 1)
  data_all <- data.frame(obs = rep(0, n_all),
                         x = rep(NA, n_all),
                         y = rep(NA, n_all),
                         time = as.POSIXct(NA, tz = tz(obs$time)),
                         stratum = rep(1:n_obs, each = n_controls + 1),
                         ID = rep(obs$ID, each = n_controls + 1))

  ## Matrix of locations
  xy <- as.matrix(obs[, c("x", "y")])

  ## Loop over observed steps
  for(i in 1:n_obs) {
    ## Index in big data frame
    ind_obs <- (i-1) * (n_controls + 1) + 1
    data_all[ind_obs, "obs"] <- 1
    data_all[ind_obs, c("ID", "x", "y", "time")] <- obs[i, c("ID", "x", "y", "time")]

    ind_controls <- (i-1) * (n_controls + 1) + (2 : (n_controls + 1))
    data_all[ind_controls, "ID"] <- obs[i, "ID"]
    data_all[ind_controls, "time"] <- obs[i, "time"]

    ## Only simulate control steps within tracks
    if(i >= 2) {
      if(obs$ID[i] == obs$ID[i-1]) {

        if(step_dist == "uniform") {
        ## Simulate uniformly from disc of radius rmax
        sim_steps <- sqrt(runif(n_controls, 0, rmax^2))
        }

        if(step_dist == "gamma") {
          sim_steps <- rgamma(n_controls,
                              shape = mean^2 / sd^2,
                              scale = sd^2 / mean)
        }

        #get uniform turning angles
        sim_bearings <- runif(n_controls, -pi, pi)

        ## Derive control steps from step lengths and bearings
        controls <- matrix(rep(xy[i-1,], each = n_controls), ncol = 2) +
          sim_steps * cbind(cos(sim_bearings), sin(sim_bearings))

        data_all[ind_controls, c("x", "y")] <- controls
      }
    }
  }

  return(data_all)
}
