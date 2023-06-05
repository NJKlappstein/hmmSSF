##' Create covariate dataframe
##' 
##' @param formula model formula with covariate names
##' @param cov_data list of covariate rasters
##' @param data location data (either matrix or dataframe with x, y)

cov_df <- function(formula, cov_data, data) { 
  
  if(is.data.frame(data)) {
    data <- matrix(c(data$x, data$y), ncol = 2)
  }
  
  # get non-movement covariates
  pattern <-  paste(c("step", "angle"), collapse="|")
  cov_names <- attr(stats::terms(formula), "term.labels") # covariate names
  cov_names <- cov_names[which(grepl(pattern, cov_names) == F)]
  
  # extract covariates
  covariates <- NULL
  if(length(cov_names) >= 1) { 
    for(j in 1:length(cov_names)){
      cov_layer <- cov_data[[which(names(cov_data) == cov_names[j])]]
      cov_extract <- extract(cov_layer, data, type = "bilinear")
      covariates <- cbind(covariates, cov_extract)
    }
  } 
  
  if(length(covariates) > 0) {colnames(covariates) <- c(cov_names)
  return(covariates)
  }
}
