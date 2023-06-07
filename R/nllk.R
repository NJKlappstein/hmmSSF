#' Negative log likelihood of the HMM-SSF (estimated with forward algorithm)
#'
#' @param par Vector of working paramters (ssf_par, tpm_par)
#' @param ssf_MM model matrix based on SSF formula
#' @param tpm_MM model matrix based on transition probability formula
#' @param ID track ID from data
#' @param statum stratum ID from data
#' @param n_states number of states in the HMM
#' @param sampling_densities dens of the controls for given sampling distribution
#' @param n_obs number of observed locations in the data
#'
#' @export

nllk <- function(par,
                 ssf_MM,
                 tpm_MM,
                 ID,
                 stratum,
                 n_states,
                 sampling_densities,
                 n_obs)
{
  # format parameters into matrices
  par <- format_par(par = par, n_states = n_states,
                    ssf_cov = colnames(ssf_MM),
                    tpm_cov = colnames(tpm_MM))

  # calculate linear predictors for the ssf
  ssf_LP <- ssf_MM %*% par$ssf

  # get state-dependent densities
  densities <- state_dens_rcpp(linear_pred = ssf_LP,
                               stratum = stratum,
                               n_states = n_states,
                               sampling_densities = sampling_densities,
                               n_obs = n_obs)
  densities <- ifelse(is.na(densities), 1, densities)

  # get array of all time-dependent transition probability matrices
  Gamma <- moveHMM:::trMatrix_rcpp(nbStates = n_states,
                                   beta = par$tpm,
                                   covs = tpm_MM)

  # get delta from stationary distribution
  delta <- solve(t(diag(n_states) - Gamma[,,1] + 1), rep(1, n_states))

  ## Run forward algorithm looped over track ID
  llk <- 0
  n_ID <- length(unique(ID)) # track names
  obs_ID <- ID[c(1, which(stratum[-1] != stratum[-length(stratum)])+1)]

  for(k in 1:n_ID) {
    # get initial density for this track
    v <- delta * densities[1,]

    # Loop over observations for this track
    for(i in which(obs_ID == unique(ID[k]))) {
      v <- v %*% Gamma[,,i] * densities[i,]
      llk <- llk + log(sum(v))
      v <- v/sum(v)
    }
  }

  return(-llk)

}
