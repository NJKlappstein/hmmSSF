##' Get the scale and scale from estimated ssf_par
##'
##' @param fit fitted model output
##' @param n_states number of fitted states
##'
##' @export

get_shape_scale <- function(fit, n_states) {

  # calculate gamma par from ssf_par
  scale <- - 1 / fit$ssf_par$estimate[which(fit$ssf_par$cov == "step")]
  shape <- fit$ssf_par$estimate[which(fit$ssf_par$cov == "log(step)")] + 2

  par_df <- data.frame(state = rep(1:n_states),
                        shape = shape,
                        scale = scale)

  return(par_df)
}
