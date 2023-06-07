
#' Print method for hmmSSF model objects
#' @method print hmmSSF
#'
#' @param mod Fitted model object, as returned by \code{\link{fitHMMSSF}}
#'
#' @export
print.hmmSSF <- function(mod) {
  cat("Negative log-likelihood:", mod$fit$value, "\n")
  cat("Convergence code:", mod$fit$convergence, "\n\n")

  # get parameters and 95% confidence intervals
  par <- confint(mod, pretty = TRUE)

  cat("SSF model:\n")
  cat(deparse(mod$args$ssf_formula))
  cat("\n")
  print(round(par$ssf, 3))

  cat("\nTPM model:\n")
  cat(deparse(mod$args$tpm_formula))
  cat("\n")
  print(round(par$tpm, 3))
}
