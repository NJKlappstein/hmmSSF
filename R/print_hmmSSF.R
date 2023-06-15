
#' Print method for hmmSSF model objects
#' @method print hmmSSF
#'
#' @param x Fitted model object, as returned by \code{\link{hmmSSF}}
#'
#' @export
print.hmmSSF <- function(x, ...) {
  mod <- x
  cat("Negative log-likelihood:", mod$fit$minimum, "\n")
  cat("Convergence code:", mod$fit$code, "\n\n")

  # get parameters and 95% confidence intervals
  par <- confint(mod, pretty = TRUE)

  cat("SSF model:\n")
  cat(deparse(mod$args$ssf_formula))
  cat("\n")
  print(round(par$ssf, 3))

  cat("\nTPM model:\n")
  cat(deparse(mod$args$tpm_formula))
  cat("\n")
  if(deparse(mod$args$tpm_formula) == "~1") {
    tpm <- predict_tpm(mod = mod,
                       new_data = mod$args$data[1,],
                       return_CI = FALSE)
    print(round(tpm$mle[,,1], 3))
  } else {
    print(round(par$tpm, 3))
  }
}
