##' Get the scale and scale from estimated ssf_par
##'
##' @param mod fitted model output
##'
##' @export

get_shape_scale <- function(mod) {

  # calculate gamma par from ssf_par
  scale <- - 1 / mod$par$ssf$estimate[which(mod$par$ssf$cov == "step")]
  shape <- mod$par$ssf$estimate[which(mod$par$ssf$cov == "log(step)")] + 2

  par_df <- data.frame(state = 1:mod$args$n_states,
                        shape = shape,
                        scale = scale)

  return(par_df)
}
