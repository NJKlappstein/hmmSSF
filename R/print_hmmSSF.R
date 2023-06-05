
#' Print method for hmmSSF model objects
#' @method print hmmSSF
#'
#' @param mod hmmSSF model object
#'
#' @export
print.hmmSSF <- function(mod) {
  print(mod$par_CI)
}
