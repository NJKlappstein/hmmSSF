##' Fit the HMM-SSF
##'
##' @param ssf_formula model formula for ssf
##' @param tmp_formula formula for transition probabilities (default = ~1)
##' @param data data with columns for ID, stratum, obs,and covariates in formula
##' @param par0 list of starting values for parameters (betas, alphas)
##' @param n_states how many states in the model
##' @param dist distribution of the controls (default = "gamma", also "uniform")
##'
##' @export
##'
##' @useDynLib hmmSSF

fitHMMSSF <- function(ssf_formula,
                      tpm_formula = ~1,
                      data,
                      par0,
                      n_states,
                      optim_opts = list(trace = 0, maxit = 5e4)) {

  # get vector of parameters
  par <- c(par0$betas, par0$alphas)

  # order data
  data <- data[order(data$ID, data$stratum, -data$obs),]
  obs <- subset(data, obs == 1)

  # get ssf model matrix (without intercept)
  options(na.action = 'na.pass')
  ssf_MM <- model.matrix(ssf_formula, data)
  ssf_MM <- ssf_MM[,!colnames(ssf_MM) == "(Intercept)"]

  # get sampling densities for correction
  sampling_densities <- attr(data, "weights")

  # get transition probabilities model matrix
  options(na.action = 'na.pass')
  tpm_MM <- model.matrix(tpm_formula, obs)

  # get number of covariates for each formula
  n_ssf_cov <- ncol(ssf_MM)
  n_tpm_cov <- ncol(tpm_MM)

  # optimise negative log likelihood
  fit <- optim(par = par,
               fn = nllk,
               ssf_MM = ssf_MM,
               tpm_MM = tpm_MM,
               sampling_densities = sampling_densities,
               stratum = data$stratum,
               ID = data$ID,
               n_states = n_states,
               n_tpm_cov = n_tpm_cov,
               n_ssf_cov = n_ssf_cov,
               n_obs = nrow(obs),
               control = optim_opts,
               hessian = TRUE)

  # unpack, back-transform, and get CIs of fitted parameters
  par_CI <- hessian_CI(fit = fit,
                       n_states = n_states,
                       ssf_MM = ssf_MM,
                       tpm_MM = tpm_MM)

  # Save model formulation in model object
  args <- list(tpm_formula = tpm_formula,
               ssf_formula = ssf_formula,
               data = data,
               n_states = n_states)

  # Returned object
  mod <- list(par_CI = par_CI,
              fit = fit,
              args = args)
  class(mod) <- append("hmmSSF", class(mod))

  return(mod)
}
