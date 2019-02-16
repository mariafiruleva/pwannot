// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cumulative_matrix
arma::mat cumulative_matrix(const arma::mat& X, const arma::uvec& pathway_sizes, const arma::mat& reference_scores, int sampling_size, bool display_progress);
RcppExport SEXP _pwannot_cumulative_matrix(SEXP XSEXP, SEXP pathway_sizesSEXP, SEXP reference_scoresSEXP, SEXP sampling_sizeSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type pathway_sizes(pathway_sizesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type reference_scores(reference_scoresSEXP);
    Rcpp::traits::input_parameter< int >::type sampling_size(sampling_sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(cumulative_matrix(X, pathway_sizes, reference_scores, sampling_size, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// reference_matrix
arma::vec reference_matrix(const arma::mat& X, const arma::vec& ind);
RcppExport SEXP _pwannot_reference_matrix(SEXP XSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(reference_matrix(X, ind));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pwannot_cumulative_matrix", (DL_FUNC) &_pwannot_cumulative_matrix, 5},
    {"_pwannot_reference_matrix", (DL_FUNC) &_pwannot_reference_matrix, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_pwannot(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
