##' Format parameters back into matrix form
##'
##' @param par vector of parameters (ssf_par, tpm_par)
##' @param n_states number of states in the HMM
##' @param n_ssf_cov number of SSF covariates
##' @param n_tpm_cov number of TPM covariates

format_par <- function(par,
                       n_states,
                       n_ssf_cov,
                       n_tpm_cov) {

  # unpack and format ssf_par
  last_beta <- n_states * n_ssf_cov
  ssf_par <- matrix(c(par[1 : last_beta]),
                  ncol = n_states)

  # unpack and format tpm_par
  tpm_par <- matrix(c(par[(last_beta + 1) : length(par)]),
                   nrow = n_tpm_cov)

  return(list("ssf_par" = ssf_par,
              "tpm_par" = tpm_par))
}
