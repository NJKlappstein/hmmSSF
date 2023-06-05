
#' Log-sum-exp utility function
logsumexp <- function(x) {
  xmax <- max(x)
  val <- xmax + log(sum(exp(x - xmax)))
  return(val)
}
