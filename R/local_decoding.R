
##' Log-sum-exp utility function
logsumexp <- function(x) {
    xmax <- max(x)
    val <- xmax + log(sum(exp(x - xmax)))
    return(val)
}

##' Forward-backward algorithm
##' (modified from hmmTMB function 'HMM$forward_backward')
##'
##' @param ssf_formula SSF formula of the HMM-SSF
##' @param tpm_formula TP formula
##' @param data data used to fit the model
##' @param fit fitted model output from fitHMMSSF
##' @param n_states number of states
##' @param dist control sampling distribution (default = "gamma")
##'
##' @return List with two elements: \code{log_forward}
##' (log forward probabilities) and \code{log_backward}
##' (log backward probabilities)
forward_backward <- function(ssf_formula,
                             tpm_formula,
                             data,
                             fit,
                             n_states,
                             dist = "gamma") {

    # get observed locations
    obs <- subset(data, obs == 1)
    n_obs <- nrow(obs)
    n_by_ID <- as.numeric(table(obs$ID))

    # separate parameters from list
    betas <-  matrix(fit$betas$estimate, ncol = n_states)
    alphas <- matrix(fit$alphas$estimate, ncol = n_states^2 - n_states)

    # get model matrix (without intercept)
    options(na.action = 'na.pass')
    ssf_MM <- model.matrix(ssf_formula, data)
    ssf_MM <- ssf_MM[,!colnames(ssf_MM) == "(Intercept)"]

    # calculate linear predictors
    ssf_LP <- ssf_MM %*% betas

    # get sampling densities
    sampling_densities <- sampling_dens(data, dist = dist, dist_par = NULL)

    # get state-dependent densities
    densities <- state_dens_rcpp(linear_pred = ssf_LP,
                                 stratum = data$stratum,
                                 n_states = n_states,
                                 sampling_densities = sampling_densities,
                                 n_obs = n_obs)
    densities <- ifelse(is.na(densities), 1, densities)

    # get Gamma from TP model matrix
    options(na.action = 'na.pass')
    tpm_MM <- model.matrix(tpm_formula, obs)
    Gamma <-  moveHMM:::trMatrix_rcpp(nbStates = n_states,
                                      beta = alphas,
                                      covs = tpm_MM)

    # get delta from Gamma
    delta <- solve(t(diag(n_states) - Gamma[,,1] + 1), rep(1, n_states))

    # initialise log-forward/backward probabilities
    log_forward <- matrix(0, nrow = n_states, nc = n_obs)
    log_backward <- matrix(0, nrow = n_states, nc = n_obs)

    # loop over ID (tracks)
    k <- 1
    for (ind in 1:length(n_by_ID)) {
        # forward algorithm
        p <- delta * densities[k,]
        psum <- sum(p)
        llk <- log(psum)
        p <- p / psum
        log_forward[, k] <- log(p) + llk
        for (i in 2:n_by_ID[ind]) {
            p <- p %*% Gamma[,, k + i - 2] * densities[k + i - 1,]
            psum <- sum(p)
            llk <- llk + log(psum)
            p <- p / psum
            log_forward[, k + i - 1] <- log(p) + llk
        }

        # backward algorithm
        log_backward[, k + n_by_ID[ind] - 1] <- rep(0, n_states)
        p <- rep(1 / n_states, n_states)
        llk <- log(n_states)
        for (i in (n_by_ID[ind] - 1):1) {
            p <- Gamma[, , k + i - 1] %*% (densities[k + i, ] * p)
            log_backward[, k + i - 1] <- log(p) + llk
            psum <- sum(p)
            p <- p / psum
            llk <- llk + log(psum)
        }

        k <- k + n_by_ID[ind]
    }

    return(list(log_forward = log_forward, log_backward = log_backward))
}

##' Local decoding using the forward-backward algorithm
##' (modified from hmmTMB function 'HMM$state_probs')
##'
##' @param ssf_formula SSF formula of the HMM-SSF
##' @param tpm_formula TP formula
##' @param data data used to fit the model
##' @param fit fitted model output from fitHMMSSF
##' @param n_states number of states
##' @param dist control sampling distribution (default = "gamma")
##'
##' @return Matrix of state probabilities, with one row for each
##' observation time, and one column for each state. The (i, j)-th
##' element is the probability of being in state j at time i.
##'
##' @export
local_decoding <- function(ssf_formula,
                           tpm_formula,
                           data,
                           fit,
                           n_states,
                           dist = "gamma") {

    # get observed locations
    obs <- subset(data, obs == 1)
    n_obs <- nrow(obs)
    n_by_ID <- as.numeric(table(obs$ID))
    cumul_n <- cumsum(n_by_ID)

    # log forward and backward probabilities
    fb <- forward_backward(ssf_formula = ssf_formula, tpm_formula = tpm_formula,
                           data = data, fit = fit, n_states = n_states,
                           dist = dist)
    log_forward <- fb$log_forward
    log_backward <- fb$log_backward

    # initialise matrix of state probabilities
    state_probs <- matrix(0, nrow = n_obs, ncol = n_states)

    # loop over tracks
    k <- 0
    for (ind in 1:length(n_by_ID)) {
        llk <- logsumexp(log_forward[, cumul_n[ind]])
        # loop over time steps
        for (i in 1:n_by_ID[ind]) {
            state_probs[k + i,] <-
                exp(log_forward[, k + i] + log_backward[, k + i] - llk)
        }
        k <- k + n_by_ID[ind]
    }

    colnames(state_probs) <- paste0("state", 1:n_states)
    return(state_probs)
}
