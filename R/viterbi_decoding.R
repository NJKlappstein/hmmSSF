
#' Generate most probable state sequence using the viterbi algorithm
#' (modified from moveHMM function 'viterbi')
#'
#' @param mod Fitted model object, as returned by \code{\link{fitHMMSSF}}
#'
#' @export

viterbi_decoding <- function(mod) {
  # unpack model
  ssf_formula <- mod$args$ssf_formula
  tpm_formula <- mod$args$tpm_formula
  data <- mod$args$data
  n_states <- mod$args$n_states

  # get observed locations
  obs <- subset(data, obs == 1)
  n_obs <- nrow(obs)

  # separate parameters from list
  ssf_par <-  mod$par$ssf
  tpm_par <- mod$par$tpm

  # get model matrix (without intercept)
  options(na.action = 'na.pass')
  ssf_MM <- model.matrix(ssf_formula, data)
  ssf_MM <- ssf_MM[,!colnames(ssf_MM) == "(Intercept)"]

  # calculate linear predictors
  ssf_LP <- ssf_MM %*% ssf_par

  # get sampling densities
  sampling_densities <- attr(data, "weights")

  # get state-dependent densities
  densities <- state_dens_rcpp(linear_pred = ssf_LP,
                               stratum = data$stratum,
                               n_states = n_states,
                               sampling_densities = sampling_densities,
                               n_obs = n_obs)
  densities <- ifelse(is.na(densities), 1, densities)

  # get Gamma from TP model matrix
  options(na.action = 'na.pass')
  tpm_MM <- model.matrix(tpm_formula, obs)
  Gamma <-  moveHMM:::trMatrix_rcpp(nbStates = n_states,
                                        beta = tpm_par,
                                        covs = tpm_MM)

  # get delta from Gamma
  delta <- solve(t(diag(n_states) - Gamma[,,1] + 1), rep(1, n_states))

  # create matrix for (what is xi?)
  xi <- matrix(0, nrow = n_obs, ncol = n_states)
  v <- delta * densities[1,]
  xi[1,] <- v / sum(v)
  for (t in 2:n_obs) {
    v <- apply(xi[t-1,] * Gamma[,,t], 2, max) * densities[t,]
    xi[t,] <- v / sum(v) }

  # most probable state sequence
  states <- rep(NA, n_obs)
  states[n_obs] <- which.max(xi[n_obs,])

  for (t in (n_obs-1):1) {
    states[t] <- which.max(Gamma[, states[t+1], t+1] * xi[t,])
  }

  return(states)
}

