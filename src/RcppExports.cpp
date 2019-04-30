// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// csFs
Rcpp::List csFs(std::vector<double> uIn, std::vector<double> vIn, std::vector<double> wIn, const double& ZIn, const double& ustarIn, const double& LinvIn, const double& ZoIn, const double& bwIn, const double& sUustarIn, const double& sVustarIn, const double& kvIn, const double& C0In, const double& alphaIn, const double& MaxFetchIn, int distThresh);
RcppExport SEXP _traxr_csFs(SEXP uInSEXP, SEXP vInSEXP, SEXP wInSEXP, SEXP ZInSEXP, SEXP ustarInSEXP, SEXP LinvInSEXP, SEXP ZoInSEXP, SEXP bwInSEXP, SEXP sUustarInSEXP, SEXP sVustarInSEXP, SEXP kvInSEXP, SEXP C0InSEXP, SEXP alphaInSEXP, SEXP MaxFetchInSEXP, SEXP distThreshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type uIn(uInSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type vIn(vInSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type wIn(wInSEXP);
    Rcpp::traits::input_parameter< const double& >::type ZIn(ZInSEXP);
    Rcpp::traits::input_parameter< const double& >::type ustarIn(ustarInSEXP);
    Rcpp::traits::input_parameter< const double& >::type LinvIn(LinvInSEXP);
    Rcpp::traits::input_parameter< const double& >::type ZoIn(ZoInSEXP);
    Rcpp::traits::input_parameter< const double& >::type bwIn(bwInSEXP);
    Rcpp::traits::input_parameter< const double& >::type sUustarIn(sUustarInSEXP);
    Rcpp::traits::input_parameter< const double& >::type sVustarIn(sVustarInSEXP);
    Rcpp::traits::input_parameter< const double& >::type kvIn(kvInSEXP);
    Rcpp::traits::input_parameter< const double& >::type C0In(C0InSEXP);
    Rcpp::traits::input_parameter< const double& >::type alphaIn(alphaInSEXP);
    Rcpp::traits::input_parameter< const double& >::type MaxFetchIn(MaxFetchInSEXP);
    Rcpp::traits::input_parameter< int >::type distThresh(distThreshSEXP);
    rcpp_result_gen = Rcpp::wrap(csFs(uIn, vIn, wIn, ZIn, ustarIn, LinvIn, ZoIn, bwIn, sUustarIn, sVustarIn, kvIn, C0In, alphaIn, MaxFetchIn, distThresh));
    return rcpp_result_gen;
END_RCPP
}
// csFi
Rcpp::List csFi(std::vector<double> uIn, std::vector<double> vIn, std::vector<double> wIn, const double& ZIn, const double& ustarIn, const double& LinvIn, const double& ZoIn, const double& bwIn, const double& sUustarIn, const double& sVustarIn, const double& kvIn, const double& C0In, const double& alphaIn, const double& MaxFetchIn, int distThresh);
RcppExport SEXP _traxr_csFi(SEXP uInSEXP, SEXP vInSEXP, SEXP wInSEXP, SEXP ZInSEXP, SEXP ustarInSEXP, SEXP LinvInSEXP, SEXP ZoInSEXP, SEXP bwInSEXP, SEXP sUustarInSEXP, SEXP sVustarInSEXP, SEXP kvInSEXP, SEXP C0InSEXP, SEXP alphaInSEXP, SEXP MaxFetchInSEXP, SEXP distThreshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type uIn(uInSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type vIn(vInSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type wIn(wInSEXP);
    Rcpp::traits::input_parameter< const double& >::type ZIn(ZInSEXP);
    Rcpp::traits::input_parameter< const double& >::type ustarIn(ustarInSEXP);
    Rcpp::traits::input_parameter< const double& >::type LinvIn(LinvInSEXP);
    Rcpp::traits::input_parameter< const double& >::type ZoIn(ZoInSEXP);
    Rcpp::traits::input_parameter< const double& >::type bwIn(bwInSEXP);
    Rcpp::traits::input_parameter< const double& >::type sUustarIn(sUustarInSEXP);
    Rcpp::traits::input_parameter< const double& >::type sVustarIn(sVustarInSEXP);
    Rcpp::traits::input_parameter< const double& >::type kvIn(kvInSEXP);
    Rcpp::traits::input_parameter< const double& >::type C0In(C0InSEXP);
    Rcpp::traits::input_parameter< const double& >::type alphaIn(alphaInSEXP);
    Rcpp::traits::input_parameter< const double& >::type MaxFetchIn(MaxFetchInSEXP);
    Rcpp::traits::input_parameter< int >::type distThresh(distThreshSEXP);
    rcpp_result_gen = Rcpp::wrap(csFi(uIn, vIn, wIn, ZIn, ustarIn, LinvIn, ZoIn, bwIn, sUustarIn, sVustarIn, kvIn, C0In, alphaIn, MaxFetchIn, distThresh));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_traxr_csFs", (DL_FUNC) &_traxr_csFs, 15},
    {"_traxr_csFi", (DL_FUNC) &_traxr_csFi, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_traxr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
