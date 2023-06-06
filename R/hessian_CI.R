##' Calculate CIs for the HMM-SSF parameters via the Hessian matrix
##'
##' @param fit fitted model output from optim
##' @param CI_range range of the CIs (default = 0.95)
##' @param n_states number of fitted states
##' @param ssf_MM model matrix for SSF
##' @param tpm_MM model matrix of transition probabilities
##'
##' @importFrom MASS ginv

hessian_CI <- function(fit, n_states, CI_range = 0.95, ssf_MM, tpm_MM) {

  # get number of covariates for each formula
  ssf_cov <- colnames(ssf_MM)
  tpm_cov <- colnames(tpm_MM)

  # get estimates
  estimates <- format_par(fit$par,
                          n_states,
                          length(ssf_cov),
                          length(tpm_cov))

  # inverse of Hessian
  Sigma <- ginv(fit$hessian)
  var <- diag(Sigma)

  # define appropriate quantile
  quant <- qnorm(1-(1-CI_range)/2)

  # compute lower and upper for working parameters
  lower_working <- fit$par - quant * sqrt(var)
  upper_working <- fit$par + quant * sqrt(var)

  # transform to natural parameters
  lower <- format_par(lower_working,
                      n_states,
                      length(ssf_cov),
                      length(tpm_cov))

  upper <- format_par(upper_working,
                     n_states,
                     length(ssf_cov),
                     length(tpm_cov))


  # format ssf_par
  ssf_par <- data.frame(state = rep(seq(1:n_states), each = length(ssf_cov)),
                      cov = rep(ssf_cov, times = n_states),
                      estimate = as.vector(estimates$ssf_par),
                      lower = as.vector(lower$ssf_par),
                      upper = as.vector(upper$ssf_par))

  #format tpm_par
  transitions <- NULL
  for(i in 1:n_states) {
    transitions_i <- paste(i, seq(1:n_states)[-i], sep = "-")
    transitions <- c(transitions, transitions_i)
  }
  tpm_par <- data.frame(transition = rep(transitions, each = length(tpm_cov)),
                       cov = rep(tpm_cov, times = n_states),
                       estimate = as.vector(estimates$tpm_par),
                       lower = as.vector(lower$tpm_par),
                       upper = as.vector(upper$tpm_par))


  CI <- list("ssf_par" = ssf_par, "tpm_par" = tpm_par)

  return(CI)
}
