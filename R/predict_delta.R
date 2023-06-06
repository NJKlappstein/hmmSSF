
##' Predict stationary state probabilities (with uncertainty)
##'
##' @param ...
##'
##' @return List with elements 'mle', 'upper', and 'lower'.
##' Each element is an array, where each layer is a transition
##' probability matrix corresponding to a row of 'new_data'.
##'
##' @export
predict_delta <- function(tpm_formula,
                        new_data,
                        fit,
                        n_states,
                        return_CI = TRUE) {

    # separate parameters from list
    tpm_par <- matrix(fit$tpm_par$estimate, ncol = n_states^2 - n_states)

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
        par_ind <- (nrow(fit$ssf_par) + 1) : (nrow(fit$ssf_par) + nrow(fit$tpm_par))
        Sigma <- solve(fit$hessian)[par_ind, par_ind]

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
