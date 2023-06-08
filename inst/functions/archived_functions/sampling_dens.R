##' Get the densities for each location given the sampling distribution
##'
##' @param data data with column for step length
##' @param dist distribution of the controls
##' @param dist_par a vector of distribution parameters (default = NULL)
##'
##' @export

sampling_dens <- function(data, dist, dist_par = NULL) {

  #get mu/sd of step data
  mean <- mean(subset(data, obs == 1)$step, na.rm = T)
  sd   <- sd(subset(data, obs ==1)$step, na.rm = T)

  # no correction needed for uniform - just use 1
  if(dist == "uniform") {
    densities <- rep(1, times = nrow(data))
  }

  # get pdf of gamma distribution used to generate steps
  if(dist == "gamma") {
    if(is.null(dist_par)) {
      densities <- dgamma(data$step,
                          shape = mean^2 / sd^2,
                          scale = sd^2 / mean)
    } else {
      densities <- dgamma(data$step,
                          shape = dist_par$shape,
                          scale = dist_par$scale)
    }
    densities <- densities / (2 * pi * data$step)
  }

  return(densities)
}
