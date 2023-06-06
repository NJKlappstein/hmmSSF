
#' Local decoding using the forward-backward algorithm
#' (modified from hmmTMB function 'HMM$state_probs')
#'
#' @param mod Fitted model object, as returned by \code{\link{fitHMMSSF}}
#'
#' @return Matrix of state probabilities, with one row for each
#' observation time, and one column for each state. The (i, j)-th
#' element is the probability of being in state j at time i.
#'
#' @export
local_decoding <- function(mod) {
  # unpack model
  data <- mod$args$data
  n_states <- mod$args$n_states

  # get observed locations
  obs <- subset(data, obs == 1)
  n_obs <- nrow(obs)
  n_by_ID <- as.numeric(table(obs$ID))
  cumul_n <- cumsum(n_by_ID)

  # log forward and backward probabilities
  fb <- forward_backward(mod = mod)
  log_forward <- fb$log_forward
  log_backward <- fb$log_backward

  # initialise matrix of state probabilities
  state_probs <- matrix(0, nrow = n_obs, ncol = n_states)

  # loop over tracks
  k <- 0
  for (ind in 1:length(n_by_ID)) {
    llk <- logsumexp(log_forward[, cumul_n[ind]])
    # loop over time steps
    for (i in 1:n_by_ID[ind]) {
      state_probs[k + i,] <-
        exp(log_forward[, k + i] + log_backward[, k + i] - llk)
    }
    k <- k + n_by_ID[ind]
  }

  colnames(state_probs) <- paste0("state", 1:n_states)
  return(state_probs)
}
