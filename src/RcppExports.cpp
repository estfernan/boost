// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// boost_gp
Rcpp::List boost_gp(arma::mat Y, arma::mat dist, IntegerMatrix nei, NumericVector s, int iter, int burn, double init_b_sigma, double init_h, double update_prop);
RcppExport SEXP _boost_boost_gp(SEXP YSEXP, SEXP distSEXP, SEXP neiSEXP, SEXP sSEXP, SEXP iterSEXP, SEXP burnSEXP, SEXP init_b_sigmaSEXP, SEXP init_hSEXP, SEXP update_propSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist(distSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type nei(neiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< double >::type init_b_sigma(init_b_sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type init_h(init_hSEXP);
    Rcpp::traits::input_parameter< double >::type update_prop(update_propSEXP);
    rcpp_result_gen = Rcpp::wrap(boost_gp(Y, dist, nei, s, iter, burn, init_b_sigma, init_h, update_prop));
    return rcpp_result_gen;
END_RCPP
}
// BOOST_Ising_MCMC_cpp
List BOOST_Ising_MCMC_cpp(arma::vec p, IntegerMatrix A, double mean_omega0, double sigma_omega0, double mean_theta, double sigma_theta, int n_iter, int burn, int M, double tau, double omega0_init, double theta_init);
RcppExport SEXP _boost_BOOST_Ising_MCMC_cpp(SEXP pSEXP, SEXP ASEXP, SEXP mean_omega0SEXP, SEXP sigma_omega0SEXP, SEXP mean_thetaSEXP, SEXP sigma_thetaSEXP, SEXP n_iterSEXP, SEXP burnSEXP, SEXP MSEXP, SEXP tauSEXP, SEXP omega0_initSEXP, SEXP theta_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type mean_omega0(mean_omega0SEXP);
    Rcpp::traits::input_parameter< double >::type sigma_omega0(sigma_omega0SEXP);
    Rcpp::traits::input_parameter< double >::type mean_theta(mean_thetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_theta(sigma_thetaSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type omega0_init(omega0_initSEXP);
    Rcpp::traits::input_parameter< double >::type theta_init(theta_initSEXP);
    rcpp_result_gen = Rcpp::wrap(BOOST_Ising_MCMC_cpp(p, A, mean_omega0, sigma_omega0, mean_theta, sigma_theta, n_iter, burn, M, tau, omega0_init, theta_init));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_boost_boost_gp", (DL_FUNC) &_boost_boost_gp, 9},
    {"_boost_BOOST_Ising_MCMC_cpp", (DL_FUNC) &_boost_BOOST_Ising_MCMC_cpp, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_boost(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
