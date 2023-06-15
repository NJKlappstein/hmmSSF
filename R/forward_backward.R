
#' Forward-backward algorithm
#' (modified from hmmTMB function 'HMM$forward_backward')
#'
#' @param mod Fitted model object, as returned by \code{\link{hmmSSF}}
#'
#' @return List with two elements: \code{log_forward}
#' (log forward probabilities) and \code{log_backward}
#' (log backward probabilities)

forward_backward <- function(mod) {
  # unpack model
  ssf_formula <- mod$args$ssf_formula
  tpm_formula <- mod$args$tpm_formula
  data <- mod$args$data
  n_states <- mod$args$n_states

  # get observed locations
  obs <- subset(data, obs == 1)
  n_obs <- nrow(obs)
  n_by_ID <- as.numeric(table(obs$ID))

  # separate parameters from list
  ssf_par <- mod$par$ssf
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

  # initialise log-forward/backward probabilities
  log_forward <- matrix(0, nrow = n_states, ncol = n_obs)
  log_backward <- matrix(0, nrow = n_states, ncol = n_obs)

  # loop over ID (tracks)
  k <- 1
  for (ind in 1:length(n_by_ID)) {
    # forward algorithm
    p <- delta * densities[k,]
    psum <- sum(p)
    llk <- log(psum)
    p <- p / psum
    log_forward[, k] <- log(p) + llk
    for (i in 2:n_by_ID[ind]) {
      p <- p %*% Gamma[,, k + i - 2] * densities[k + i - 1,]
      psum <- sum(p)
      llk <- llk + log(psum)
      p <- p / psum
      log_forward[, k + i - 1] <- log(p) + llk
    }

    # backward algorithm
    log_backward[, k + n_by_ID[ind] - 1] <- rep(0, n_states)
    p <- rep(1 / n_states, n_states)
    llk <- log(n_states)
    for (i in (n_by_ID[ind] - 1):1) {
      p <- Gamma[, , k + i - 1] %*% (densities[k + i, ] * p)
      log_backward[, k + i - 1] <- log(p) + llk
      psum <- sum(p)
      p <- p / psum
      llk <- llk + log(psum)
    }

    k <- k + n_by_ID[ind]
  }

  return(list(log_forward = log_forward, log_backward = log_backward))
}
