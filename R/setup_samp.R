
#' Setup importance sampling distributions
setup_samp <- function(dist = "uniform", obs) {
  move_data <- moveHMM::prepData(obs, type = "UTM")
  step <- na.omit(move_data$step)
  angle <- na.omit(move_data$angle)

  if(dist[1] == "uniform") {
    # Uniform distribution in 2D (not uniform step lengths!)
    par_step <- c(0, max(step))
    r_step <- function(n, par) {
      sqrt(runif(n, min = par[1], max = par[2]^2))
    }
    d_step <- function(x, par) {
      rep(1, length(x))
    }
  } else if(dist[1] == "gamma") {
    # Gamma distribution of step lengths
    par_step <- c(mean(step), sd(step))
    r_step <- function(n, par) {
      rgamma(n, shape = par[1]^2/par[2]^2, rate = par[1]/par[2]^2)
    }
    d_step <- function(x, par) {
      dgamma(x, shape = par[1]^2/par[2]^2, rate = par[1]/par[2]^2) / x
    }
  } else if(dist[1] == "exp") {
    # Exponential distribution of step lengths
    par_step <- 1/mean(step)
    r_step <- function(n, par){
      rexp(n, rate = par[1])
    }
    d_step <- function(x, par) {
      dexp(x, rate = par[1]) / x
    }
  } else {
    stop("'dist' is not implemented")
  }

  if(length(dist) == 1) {
    # By default, uniform distribution of turning angles
    r_angle <- function(n, par) runif(n, -pi, pi)
    d_angle <- function(x, par) 1
    par_angle <- NA
  } else if(length(dist) == 2) {
    if(dist[2] == "vm") {
      par_angle <- unlist(CircStats::vm.ml(angle))
      r_angle <- function(n, par) {
        CircStats::rvm(n = n, mean = par[1], k = par[2])
      }
      d_angle <- function(x, par) {
        CircStats::dvm(theta = x, mu = par[1], kappa = par[2])
      }
    } else if(dist[2] == "wrpcauchy") {
      par_angle <- unlist(CircStats::wrpcauchy.ml(angle, mu0 = 0, rho0 = 0.5))
      r_angle <- function(n, par) {
        CircStats::rwrpcauchy(n = n, location = par[1], rho = par[2])
      }
      d_angle <- function(x, par) {
        CircStats::dwrpcauchy(theta = x, mu = par[1], rho = par[2])
      }
    } else {
      stop("'dist' is not implemented")
    }
  } else {
    step("'dist' should have length 1 or 2")
  }

  return(list(r_step = r_step, d_step = d_step, par_step = par_step,
              r_angle = r_angle, d_angle = d_angle, par_angle = par_angle))
}
