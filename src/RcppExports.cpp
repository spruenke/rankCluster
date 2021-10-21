// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// f_psi_cpp
arma::vec f_psi_cpp(arma::vec x, Rcpp::List data, arma::vec psi);
RcppExport SEXP _rankCluster_f_psi_cpp(SEXP xSEXP, SEXP dataSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(f_psi_cpp(x, data, psi));
    return rcpp_result_gen;
END_RCPP
}
// f_theta_cpp
arma::vec f_theta_cpp(arma::vec x, Rcpp::List data, arma::vec theta, Rcpp::List psi);
RcppExport SEXP _rankCluster_f_theta_cpp(SEXP xSEXP, SEXP dataSEXP, SEXP thetaSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(f_theta_cpp(x, data, theta, psi));
    return rcpp_result_gen;
END_RCPP
}
// rel_eff_cpp
arma::vec rel_eff_cpp(Rcpp::List data, arma::vec theta, Rcpp::List psi);
RcppExport SEXP _rankCluster_rel_eff_cpp(SEXP dataSEXP, SEXP thetaSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(rel_eff_cpp(data, theta, psi));
    return rcpp_result_gen;
END_RCPP
}
// sigma_est_cpp
arma::mat sigma_est_cpp(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi);
RcppExport SEXP _rankCluster_sigma_est_cpp(SEXP nSEXP, SEXP dataSEXP, SEXP thetaSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_est_cpp(n, data, theta, psi));
    return rcpp_result_gen;
END_RCPP
}
// q_wald_arma
Rcpp::List q_wald_arma(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi, arma::mat cmat);
RcppExport SEXP _rankCluster_q_wald_arma(SEXP nSEXP, SEXP dataSEXP, SEXP thetaSEXP, SEXP psiSEXP, SEXP cmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cmat(cmatSEXP);
    rcpp_result_gen = Rcpp::wrap(q_wald_arma(n, data, theta, psi, cmat));
    return rcpp_result_gen;
END_RCPP
}
// q_anova_arma
Rcpp::List q_anova_arma(arma::vec n, Rcpp::List data, arma::vec theta, Rcpp::List psi, arma::mat cmat);
RcppExport SEXP _rankCluster_q_anova_arma(SEXP nSEXP, SEXP dataSEXP, SEXP thetaSEXP, SEXP psiSEXP, SEXP cmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cmat(cmatSEXP);
    rcpp_result_gen = Rcpp::wrap(q_anova_arma(n, data, theta, psi, cmat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rankCluster_f_psi_cpp", (DL_FUNC) &_rankCluster_f_psi_cpp, 3},
    {"_rankCluster_f_theta_cpp", (DL_FUNC) &_rankCluster_f_theta_cpp, 4},
    {"_rankCluster_rel_eff_cpp", (DL_FUNC) &_rankCluster_rel_eff_cpp, 3},
    {"_rankCluster_sigma_est_cpp", (DL_FUNC) &_rankCluster_sigma_est_cpp, 4},
    {"_rankCluster_q_wald_arma", (DL_FUNC) &_rankCluster_q_wald_arma, 5},
    {"_rankCluster_q_anova_arma", (DL_FUNC) &_rankCluster_q_anova_arma, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_rankCluster(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
