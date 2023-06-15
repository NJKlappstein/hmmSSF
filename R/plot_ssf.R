#' Plot RSS
#'
#' @param mod Fitted model returned by \code{\link{hmmSSF}}
#' @param var Name of variable to plot SSF against
#'
#' @importFrom ggplot2 ggplot geom_line facet_wrap labs aes
#'
#' @export

plot_ssf <- function(mod, var) {
  # all variable names
  vars <- all.vars(mod$args$ssf_formula)

  # new data frame for covariate grid
  new_data <- as.data.frame(matrix(0, nrow = 1000, ncol = length(vars)))
  colnames(new_data) <- vars
  for(i in 1:length(vars)) {
    if(vars[i] != var) {
      # fix to first value
      new_data[,vars[i]] <- mod$args$data[1, vars[i]]
    } else {
      # define range of covariate
      range <- range(mod$args$data[, vars[i]], na.rm = TRUE)
      new_data[,vars[i]] <- seq(range[1], range[2], length = 1000)
    }
  }

  # predict SSF over grid
  ssf <- predict_ssf(mod = mod, new_data = new_data)
  ssf <- apply(ssf, 2, function(col) col / mean(col))

  df <- data.frame(var = new_data[[var]],
                   state = paste0("state ", rep(1:mod$args$n_states, each = 1000)),
                   ssf = as.vector(ssf))

  p <- ggplot(df, aes(var, ssf)) +
    geom_line() +
    facet_wrap("state", scales = "free_y") +
    labs(x = var, y = "RSS")

  return(p)
}
