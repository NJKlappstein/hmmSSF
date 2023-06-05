##' Format parameters back into matrix form
##' 
##' @param par vector of parameters (betas, alphas)
##' @param n_states number of states in the HMM
##' @param n_ssf_cov number of SSF covariates
##' @param n_tpm_cov number of TPM covariates

format_par <- function(par, 
                       n_states, 
                       n_ssf_cov,  
                       n_tpm_cov) {

  # unpack and format betas
  last_beta <- n_states * n_ssf_cov
  betas <- matrix(c(par[1 : last_beta]),
                  ncol = n_states)
  
  # unpack and format alphas
  alphas <- matrix(c(par[(last_beta + 1) : length(par)]), 
                   nrow = n_tpm_cov)
  
  return(list("betas" = betas, 
              "alphas" = alphas))
}
