
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
