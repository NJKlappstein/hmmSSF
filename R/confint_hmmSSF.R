#' Calculate CIs for the HMM-SSF parameters via the Hessian matrix
#' @method confint hmmSSF
#'
#' @param mod Fitted model output returned by \code{\link{fitHMMSSF}}
#' @param range range of the CIs (default = 0.95)
#'
#' @importFrom MASS ginv
#'
#' @export

confint.hmmSSF <- function(mod, range = 0.95) {
  # inverse of Hessian
  Sigma <- ginv(mod$fit$hessian)
  var <- diag(Sigma)

  # define appropriate quantile
  quant <- qnorm(1 - (1 - range) / 2)

  # compute lower and upper bounds
  mle <- mod$fit$par
  low <- mle - quant * sqrt(var)
  upp <- mle + quant * sqrt(var)

  # format as matrices
  low_mat <- format_par(par = low, n_states = mod$args$n_states,
                        ssf_cov = mod$args$ssf_cov,
                        tpm_cov = mod$args$tpm_cov)
  upp_mat <- format_par(par = upp, n_states = mod$args$n_states,
                        ssf_cov = mod$args$ssf_cov,
                        tpm_cov = mod$args$tpm_cov)

  return(list(low = low_mat, upp = upp_mat))
}
