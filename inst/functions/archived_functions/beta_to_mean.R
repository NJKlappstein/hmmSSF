##' Get the mean and sd from estimated ssf_par
##'
##' @param fit fitted model output
##' @param n_states number of fitted states
##'
##' @export

beta_to_mean <- function(fit, n_states) {

  # calculate mean from ssf_par
  mean_step <- (fit$ssf_par$estimate[which(fit$ssf_par$cov == "log(step)")] + 2) /
    fit$ssf_par$estimate[which(fit$ssf_par$cov == "step")] * - 1

  # calculate sd from ssf_par
  sd_step <- sqrt(fit$ssf_par$estimate[which(fit$ssf_par$cov == "log(step)")] + 2) /
    fit$ssf_par$estimate[which(fit$ssf_par$cov == "step")] * - 1

  # create dataframe to return
  mean_df <- data.frame(state = rep(1:n_states),
                        mean = mean_step,
                        sd = sd_step)

  return(mean_df)
}

