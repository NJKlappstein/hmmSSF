
#' Predict stationary state probabilities (with uncertainty)
#'
#' @param mod Fitted model object, as returned by \code{\link{hmmSSF}}
#' @param new_data Data frame with covariate values to use for prediction
#' @param return_CI Logical (default: TRUE). Should confidence intervals be
#' computed?
#'
#' @return List with elements 'mle', 'upper', and 'lower'.
#' Each element is an array, where each layer is a transition
#' probability matrix corresponding to a row of 'new_data'.
#'
#' @export

predict_delta <- function(mod,
                          new_data,
                          return_CI = TRUE) {

  n_states <- mod$args$n_states
  tpm_formula <- mod$args$tpm_formula

  # separate parameters from list
  tpm_par <- mod$par$tpm

  # get Gamma from TP model matrix
  options(na.action = 'na.pass')
  tpm_MM <- model.matrix(tpm_formula, new_data)
  Gamma <-  moveHMM:::trMatrix_rcpp(nbStates = n_states,
                                    beta = tpm_par,
                                    covs = tpm_MM)
  delta <- t(apply(Gamma, 3, function(tpm) {
    solve(t(diag(n_states) - tpm + 1), rep(1, n_states))
  }))
  colnames(delta) <- 1:n_states
  out <- list(mle = delta)

  if(return_CI) {
    # get covariance matrix
    par_ind <- (length(mod$par$ssf) + 1) :
      (length(mod$par$ssf) + length(mod$par$tpm))
    Sigma <- solve(mod$fit$hessian)[par_ind, par_ind]

    # for differentiation to obtain confidence intervals (delta method)
    get_delta <- function(alpha, modmat, n_states, i) {
      tpm <- moveHMM:::trMatrix_rcpp(nbStates = n_states,
                                     beta = alpha,
                                     covs = modmat)[,,1]
      delta <- solve(t(diag(n_states) - tpm + 1), rep(1, n_states))
      return(delta[i])
    }

    # copy structure of Gamma for upper/lower CI bounds
    lower <- delta
    upper <- delta

    # loop over transition probabilities
    for(i in 1:n_states) {
      # derive confidence intervals using the delta method
      dN <- t(apply(tpm_MM, 1, function(x)
        numDeriv::grad(get_delta, tpm_par,
                       modmat = matrix(x, nrow = 1),
                       n_states = n_states, i = i)))

      se <- t(apply(dN, 1, function(x)
        suppressWarnings(sqrt(x %*% Sigma %*% x))))

      # transform estimates and standard errors to R, to derive CI
      # on working scale, then back-transform to [0,1]
      width <- 1.96 * se / (delta[,i] - delta[,i]^2)
      lower[,i] <- plogis(qlogis(delta[,i]) - width)
      upper[,i] <- plogis(qlogis(delta[,i]) + width)
    }

    out$lower <- lower
    out$upper <- upper
  }

  return(out)
}
