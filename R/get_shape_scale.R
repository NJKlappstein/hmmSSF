##' Get the scale and scale from estimated betas 
##' 
##' @param fit fitted model output
##' @param n_states number of fitted states

get_shape_scale <- function(fit, n_states) {
  
  # calculate gamma par from betas
  scale <- - 1 / fit$betas$estimate[which(fit$betas$cov == "step")]
  shape <- fit$betas$estimate[which(fit$betas$cov == "log(step)")] + 2
  
  par_df <- data.frame(state = rep(1:n_states), 
                        shape = shape, 
                        scale = scale)
  
  return(par_df)
}
