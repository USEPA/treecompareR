// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// jaccard_index
double jaccard_index(NumericVector s1, NumericVector s2);
RcppExport SEXP _treecompareR_jaccard_index(SEXP s1SEXP, SEXP s2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2(s2SEXP);
    rcpp_result_gen = Rcpp::wrap(jaccard_index(s1, s2));
    return rcpp_result_gen;
END_RCPP
}
// get_jaccard
NumericVector get_jaccard(List list1, List list2);
RcppExport SEXP _treecompareR_get_jaccard(SEXP list1SEXP, SEXP list2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type list1(list1SEXP);
    Rcpp::traits::input_parameter< List >::type list2(list2SEXP);
    rcpp_result_gen = Rcpp::wrap(get_jaccard(list1, list2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_treecompareR_jaccard_index", (DL_FUNC) &_treecompareR_jaccard_index, 2},
    {"_treecompareR_get_jaccard", (DL_FUNC) &_treecompareR_get_jaccard, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_treecompareR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
