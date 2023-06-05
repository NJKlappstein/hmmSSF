##' Simulate a movement track from the HMM-SSF
##' 
##' @param ssf_formula model formula for ssf
##' @param tpm_formula model formula for transition probabilities
##' @param ssf_cov named list of covariate rasters for ssf
##' @param tpm_cov dataframe of transition probability covariates (non-spatial)
##' @param par list of parameters:betas (matrix), alphas (matrix)
##' @param n_states how many states in the model
##' @param n_obs how many observations to generate
##' @param n_zeros number of endpoints to choose from
##' @param y1 coordinates of first location
##' @param rmax size of disc radius for generating possible endpoints

simHMMSSF <- function(ssf_formula, tpm_formula = ~1, 
                      ssf_cov = NULL, tpm_cov = NULL, 
                      par, n_states, n_obs, n_zeros, 
                      y1, rmax,
                      levels_cov = NULL, 
                      print = TRUE) { 
  
  # unpack parameters
  betas  <- par$betas
  
  # if there is no tpm covariates, then create a dummy dataframe
  if(tpm_formula == "~1") {
    tpm_cov <- data.frame(cov = rep(1, times = n_obs))
  }
  
  # get array of TPMs (Gamma)
  tpm_MM <- model.matrix(tpm_formula, tpm_cov)
  Gamma <- moveHMM:::trMatrix_rcpp(nbStates = n_states, 
                                   beta = par$alphas, 
                                   covs = tpm_MM)
  
  # get delta from stationary distribution
  delta <- solve(t(diag(n_states) - Gamma[,,1] + 1), rep(1, n_states))
  
  ## Generate state sequence ##
  S <- matrix(ncol = 1, nrow = n_obs)
  # sample initial state with probabilities given by delta
  S[1] <- sample(seq(1:n_states), 1, prob = delta)
  for(i in 2:n_obs) {
    # Sample with probabilities given by the previous states TPs
    S[i] <- sample(1:n_states, 1, prob = Gamma[S[i-1],,i])
  }
  
  ## Generate locations ##
  xy <- matrix(NA, nrow = n_obs, ncol = 2)   # matrix for locations
  step1 <-  sqrt(runif(1, 0, rmax^2))         # first step length
  bearing1 <- runif(1, -pi, pi)               # first bearing
  xy[1,] <- y1 + step1 * cbind(cos(bearing1), sin(bearing1)) # first loc to save
  
  i <- 1
  prev_bearing <- bearing1
  while(i <= n_obs-1) {
    #print simulation completion
    if(i%%10 == 0 & print == TRUE) {
      cat("\rSimulating...", round(i/(n_obs-1), 2)*100, "%")            
    }
    
    # generate step lengths and TAs uniformly from disc with radius rmax
    sim_steps <- sqrt(runif(n_zeros, 0, rmax^2)) 
    sim_bearings <- runif(n_zeros, -pi, pi) 
    sim_angles <- sim_bearings - prev_bearing
    
    # derive possible endpoints from step lengths and bearings
    endpoints <- matrix(rep(xy[i,], each = n_zeros), ncol = 2) +
      sim_steps * cbind(cos(sim_bearings), sin(sim_bearings))
    
    # create dataframe from movement and habitat covariates 
    move_cov <- data.frame(step = sim_steps, angle = sim_angles)
    hab_cov <- as.data.frame(cov_df(ssf_formula, ssf_cov, data = endpoints))
    if(length(hab_cov) > 0) {
      # check if any variables should be formatted as factors
      if(length(levels_cov) > 0) {
        for(k in 1:length(levels_cov)) {
          name <- names(levels_cov)[k]
          hab_cov[[name]] <- factor(hab_cov[[name]])
          levels(hab_cov[[name]]) <- levels_cov[[k]]
        }
      }
      covariates <- cbind(move_cov, hab_cov)
    } else {
      covariates <- move_cov
    }
    
    ## Remove points with NA covariates if necessary
    if(any(is.na(covariates))) {
      missing <- unique(which(is.na(covariates), arr.ind = TRUE)[, "row"])
      covariates <- covariates[-missing,]
      endpoints <- endpoints[-missing,]
      sim_bearings <- sim_bearings[-missing]
    }
    # get model matrix
    mod_matrix <- model.matrix(ssf_formula, covariates)
    mod_matrix <- mod_matrix[,!colnames(mod_matrix) == "(Intercept)"]
    
    # get state-specific linear predictor
    linear_pred <- mod_matrix %*% as.matrix(betas[,S[i+1]])
    
    # calculate ssf for each endpoint
    ssf <- exp(linear_pred) 
    
    # sample next location with probabilities prop to state-specific SSF
    ssf_prob <- ssf/(sum(ssf))
    sample <- sample(1:nrow(endpoints), size=1, prob = ssf_prob)
    xy[i+1,] <- endpoints[sample,]
    
    prev_bearing <- sim_bearings[sample]
    i <- i+1
  }
  return(data.frame(x = xy[,1], 
                    y = xy[,2],
                    state = S))
}
