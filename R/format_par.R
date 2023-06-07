#' Format parameters back into matrix form
#'
#' @param par Vector of parameters (ssf_par, tpm_par)
#' @param n_states Number of states in the HMM
#' @param ssf_cov Names of SSF covariates
#' @param tpm_cov Names of TPM covariates

format_par <- function(par,
                       n_states,
                       ssf_cov,
                       tpm_cov) {

  n_ssf_cov <- length(ssf_cov)
  n_tpm_cov <- length(tpm_cov)

  # unpack and format ssf_par
  last_beta <- n_states * n_ssf_cov
  ssf_par <- matrix(c(par[1 : last_beta]),
                    ncol = n_states)
  colnames(ssf_par) <- paste0("S", 1:n_states)
  rownames(ssf_par) <- ssf_cov

  # unpack and format tpm_par
  tpm_par <- matrix(c(par[(last_beta + 1) : length(par)]),
                    nrow = n_tpm_cov)
  tpm_names <- paste0(rep(1:n_states, each = n_states),
                      ">", rep(1:n_states, n_states))
  colnames(tpm_par) <- tpm_names[-which(diag(n_states) == 1)]
  rownames(tpm_par) <- tpm_cov

  return(list(ssf = ssf_par,
              tpm = tpm_par))
}
