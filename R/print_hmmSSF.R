
#' Print method for hmmSSF model objects
#' @method print hmmSSF
#'
#' @param mod hmmSSF model object
#'
#' @export
print.hmmSSF <- function(mod) {
  cat("Negative log-likelihood:", mod$fit$value, "\n")
  cat("Convergence code:", mod$fit$convergence, "\n\n")

  cat("SSF model:\n")
  print(mod$args$ssf_formula)
  cat("\n")
  print(mod$par$ssf)

  cat("\n\nTPM model:\n")
  print(mod$args$tpm_formula)
  cat("\n")
  print(mod$par$tpm)
}
