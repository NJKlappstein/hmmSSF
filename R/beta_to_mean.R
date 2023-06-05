##' Get the mean and sd from estimated betas
##'
##' @param fit fitted model output
##' @param n_states number of fitted states
##'
##' @export

beta_to_mean <- function(fit, n_states) {

  # calculate mean from betas
  mean_step <- (fit$betas$estimate[which(fit$betas$cov == "log(step)")] + 2) /
    fit$betas$estimate[which(fit$betas$cov == "step")] * - 1

  # calculate sd from betas
  sd_step <- sqrt(fit$betas$estimate[which(fit$betas$cov == "log(step)")] + 2) /
    fit$betas$estimate[which(fit$betas$cov == "step")] * - 1

  # create dataframe to return
  mean_df <- data.frame(state = rep(1:n_states),
                        mean = mean_step,
                        sd = sd_step)

  return(mean_df)
}

