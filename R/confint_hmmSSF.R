#' Calculate CIs for the HMM-SSF parameters via the Hessian matrix
#' @method confint hmmSSF
#'
#' @param mod Fitted model output returned by \code{\link{hmmSSF}}
#' @param range range of the CIs (default = 0.95)
#' @param pretty controls how the CIs are stored. If set to TRUE, the CIs will
#' be stored as they appear in the printed model. Otherwise, they are stored
#' as a list with elements for ssf and tpm.
#'
#' @importFrom MASS ginv
#' @importFrom stats qnorm
#'
#' @export

confint.hmmSSF <- function(mod, range = 0.95, pretty = FALSE) {
  # inverse of Hessian
  Sigma <- ginv(mod$fit$hessian)
  var <- diag(Sigma)

  # define appropriate quantile
  quant <- qnorm(1 - (1 - range) / 2)

  # compute lower and upper bounds
  mle <- mod$fit$estimate
  low <- mle - quant * sqrt(var)
  upp <- mle + quant * sqrt(var)

  # format as matrices
  low_mat <- format_par(par = low, n_states = mod$args$n_states,
                        ssf_cov = mod$args$ssf_cov,
                        tpm_cov = mod$args$tpm_cov)
  upp_mat <- format_par(par = upp, n_states = mod$args$n_states,
                        ssf_cov = mod$args$ssf_cov,
                        tpm_cov = mod$args$tpm_cov)

  out <- list(ssf = list(low = low_mat$ssf, upp = upp_mat$ssf),
              tpm = list(low = low_mat$tpm, upp = upp_mat$tpm))

  if(pretty) {
    rows_ssf <- paste0(rep(rownames(out$ssf$low), 2), ".",
                       rep(colnames(out$ssf$low), each = nrow(out$ssf$low)))
    mat_ssf <- matrix(c(mod$par$ssf, out$ssf$low, out$ssf$upp), ncol = 3,
                      dimnames = list(rows_ssf, c("mle", "low", "upp")))
    rows_tpm <- paste0(rep(rownames(out$tpm$low), 2), ".",
                       rep(colnames(out$tpm$low), each = nrow(out$tpm$low)))
    mat_tpm <- matrix(c(mod$par$tpm, out$tpm$low, out$tpm$upp), ncol = 3,
                      dimnames = list(rows_tpm, c("mle", "low", "upp")))
    out <- list(ssf = mat_ssf, tpm = mat_tpm)
  }

  return(out)
}
