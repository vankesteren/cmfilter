// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// arma_cmf
arma::vec arma_cmf(arma::vec& x, arma::mat& M, arma::vec& y, int& maxIter, int& stableLag, double& critval, int& decfun, int& nCores, int& nStarts, bool& pb);
RcppExport SEXP _cmfilter_arma_cmf(SEXP xSEXP, SEXP MSEXP, SEXP ySEXP, SEXP maxIterSEXP, SEXP stableLagSEXP, SEXP critvalSEXP, SEXP decfunSEXP, SEXP nCoresSEXP, SEXP nStartsSEXP, SEXP pbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int& >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< int& >::type stableLag(stableLagSEXP);
    Rcpp::traits::input_parameter< double& >::type critval(critvalSEXP);
    Rcpp::traits::input_parameter< int& >::type decfun(decfunSEXP);
    Rcpp::traits::input_parameter< int& >::type nCores(nCoresSEXP);
    Rcpp::traits::input_parameter< int& >::type nStarts(nStartsSEXP);
    Rcpp::traits::input_parameter< bool& >::type pb(pbSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_cmf(x, M, y, maxIter, stableLag, critval, decfun, nCores, nStarts, pb));
    return rcpp_result_gen;
END_RCPP
}
// checkOMP
bool checkOMP();
RcppExport SEXP _cmfilter_checkOMP() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(checkOMP());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cmfilter_arma_cmf", (DL_FUNC) &_cmfilter_arma_cmf, 10},
    {"_cmfilter_checkOMP", (DL_FUNC) &_cmfilter_checkOMP, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_cmfilter(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
