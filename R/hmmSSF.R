#' Fit the HMM-SSF
#'
#' @param ssf_formula Model formula for ssf
#' @param tpm_formula Formula for transition probabilities (default = ~1)
#' @param n_states Number of HMM states
#' @param data Data frame with named columns for ID, stratum, obs, and
#' covariates in ssf_formula
#' @param ssf_par0 Matrix of starting values for SSF parameters, with one row
#' for each covariate and one column for each HMM state.
#' @param tpm_par0 Matrix of starting values for transition probability
#' parameters, with one row for each covariate (including the first row for
#' intercepts) and one column for each off-diagonal transition probability.
#' @param method Optimisation method used in optim (by default "Nelder-Mead)
#' but can also take any other option available in optim
#' @param optim_opts Additional options to be passed to optim, as a named list
#' (e.g., maxit, trace)
#'
#' @importFrom stats optim
#'
#' @export
#'
#' @useDynLib hmmSSF

hmmSSF <- function(ssf_formula,
                   tpm_formula = ~1,
                   n_states,
                   data,
                   ssf_par0,
                   tpm_par0 = NULL,
                   optim_opts = list(trace = 0, maxit = 5e4),
                   method = "Nelder-Mead") {

  # order data
  data <- data[order(data$ID, data$stratum, -data$obs),]
  obs <- subset(data, obs == 1)

  # get ssf model matrix (without intercept)
  options(na.action = 'na.pass')
  ssf_MM <- model.matrix(ssf_formula, data)
  ssf_MM <- ssf_MM[,!colnames(ssf_MM) == "(Intercept)"]

  # get transition probabilities model matrix
  options(na.action = 'na.pass')
  tpm_MM <- model.matrix(tpm_formula, obs)

  # get vector of parameters
  if(is.null(tpm_par0)) {
    # if tpm_par0 not provided, default to -2 for intercept and 0 elsewhere
    tpm_par0 <- matrix(0, nrow = ncol(tpm_MM), ncol = n_states * (n_states - 1))
    tpm_par0[1,] <- - 2
  }
  par <- c(ssf_par0, tpm_par0)

  # get sampling densities for correction
  #sampling_densities <- attr(data, "weights")
  sampling_densities <- data$w

  # optimise negative log likelihood
  fit <- optim(par = par,
               fn = nllk,
               ssf_MM = ssf_MM,
               tpm_MM = tpm_MM,
               sampling_densities = sampling_densities,
               stratum = data$stratum,
               ID = data$ID,
               n_states = n_states,
               n_obs = nrow(obs),
               control = optim_opts,
               hessian = TRUE,
               method = method)

  # unpack fitted parameters
  par <-  format_par(par = fit$par,
                     n_states = n_states,
                     ssf_cov = colnames(ssf_MM),
                     tpm_cov = colnames(tpm_MM))

  # save model formulation in model object
  args <- list(tpm_formula = tpm_formula,
               ssf_formula = ssf_formula,
               data = data,
               n_states = n_states,
               ssf_cov = colnames(ssf_MM),
               tpm_cov = colnames(tpm_MM))

  # returned object
  mod <- list(par = par,
              fit = fit,
              args = args)
  class(mod) <- append("hmmSSF", class(mod))

  return(mod)
}
