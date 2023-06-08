
#' Log-sum-exp utility function
logsumexp <- function(x) {
  xmax <- max(x)
  val <- xmax + log(sum(exp(x - xmax)))
  return(val)
}

#' Shift angles to be between -pi and pi
shift_angle <- function(angles) {
  angles[which(angles < -pi)] <- angles[which(angles < -pi)] + 2 * pi
  angles[which(angles > pi)] <- angles[which(angles > pi)] - 2 * pi
  return(angles)
}

#' Simulate raster
#'
#' @export
sim_raster <- function(rho = 11, lim = c(-30, 30, -30, 30), res = 1) {
  # make rho odd
  rho <- ifelse(rho %% 2 == 0, rho + 1, rho)

  # initialise raster
  r <- rast(xmin = lim[1], xmax = lim[2], ymin = lim[3], ymax = lim[4],
            resolution = res)
  r[] <- rnorm(ncell(r), mean = 0, sd = 1)

  # run moving average over square
  w <- abs(row(matrix(1, rho, rho)) - (rho + 1) / 2) +
    abs(col(matrix(1, rho, rho)) - (rho + 1) / 2) <= (rho - 1) / 2
  r <- focal(r, w = 1 * w, fun = mean)
  r <- scale(r)
  return(r)
}
