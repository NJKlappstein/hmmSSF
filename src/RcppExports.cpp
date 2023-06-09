// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// state_dens_rcpp
arma::mat state_dens_rcpp(arma::mat linear_pred, arma::vec stratum, int n_states, arma::vec sampling_densities, int n_obs);
RcppExport SEXP _hmmSSF_state_dens_rcpp(SEXP linear_predSEXP, SEXP stratumSEXP, SEXP n_statesSEXP, SEXP sampling_densitiesSEXP, SEXP n_obsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type linear_pred(linear_predSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type stratum(stratumSEXP);
    Rcpp::traits::input_parameter< int >::type n_states(n_statesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sampling_densities(sampling_densitiesSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    rcpp_result_gen = Rcpp::wrap(state_dens_rcpp(linear_pred, stratum, n_states, sampling_densities, n_obs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hmmSSF_state_dens_rcpp", (DL_FUNC) &_hmmSSF_state_dens_rcpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_hmmSSF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
