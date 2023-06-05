#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat state_dens_rcpp(arma::mat linear_pred, 
                          arma::vec stratum, 
                          int n_states, 
                          arma::vec sampling_densities,
                          int n_obs)
{
  // initialise empty variables to use in loop
  int n = stratum.size();
  arma::mat    densities(n_obs, n_states, fill::zeros);
  arma::rowvec ssf_obs(n_states);
  arma::rowvec ssf_denom(n_states);
  arma::rowvec ssf_control(n_states);
  arma::rowvec state_ssf(n_states);
  
  int strat = 0;
  int n_controls = 0;
  for(int i = 0; i < n; i++) {
    if(i == 0) {
      ssf_obs = exp(linear_pred.row(i));
      
      // set first value of the controls to the obs to include in denom
      ssf_denom = ssf_obs / sampling_densities(i);
      n_controls = 0;
      
    } else if(stratum(i) != stratum(i-1)) { 
      // get ssf for observation
      ssf_obs = exp(linear_pred.row(i));
      
      // set first value of the controls to the obs to include in denom
      ssf_denom = ssf_obs / sampling_densities(i);
      n_controls = 0;
    } else {
      // calculate the ssf of the control location and add it to denominator
      ssf_control = exp(linear_pred.row(i)) / sampling_densities(i);
      
      // add to the number of controls if it was a non-NA control 
      if(!Rcpp::NumericVector::is_na(sum(ssf_control))) {
        ssf_denom = ssf_denom + ssf_control;
        n_controls = n_controls + 1;
      }
    } 
    // if it's the last control of the stratum, calculate the ssf
    if((i == n-1) || (stratum(i) != stratum(i+1))) {
      state_ssf = ssf_obs / (ssf_denom / n_controls); 
      densities.row(strat) = state_ssf; 
      strat = strat + 1;
    }
  }
  return densities;  
    
} 
