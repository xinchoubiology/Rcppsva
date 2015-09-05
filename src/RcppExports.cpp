// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// int_eprior
Rcpp::List int_eprior(const Eigen::MatrixXd& sx, const Eigen::VectorXd& ghat, const Eigen::VectorXd& dhat);
RcppExport SEXP Rcppsva_int_eprior(SEXP sxSEXP, SEXP ghatSEXP, SEXP dhatSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type sx(sxSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type ghat(ghatSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type dhat(dhatSEXP);
    __result = Rcpp::wrap(int_eprior(sx, ghat, dhat));
    return __result;
END_RCPP
}
// arma_eigen
arma::vec arma_eigen(const arma::mat& M);
RcppExport SEXP Rcppsva_arma_eigen(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    __result = Rcpp::wrap(arma_eigen(M));
    return __result;
END_RCPP
}
// beta_regress
Rcpp::List beta_regress(const arma::mat& M, const arma::mat& pv, const arma::mat& svs, const int full);
RcppExport SEXP Rcppsva_beta_regress(SEXP MSEXP, SEXP pvSEXP, SEXP svsSEXP, SEXP fullSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pv(pvSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type svs(svsSEXP);
    Rcpp::traits::input_parameter< const int >::type full(fullSEXP);
    __result = Rcpp::wrap(beta_regress(M, pv, svs, full));
    return __result;
END_RCPP
}
// bootstrap_regress
Rcpp::List bootstrap_regress(const arma::mat& M, const arma::mat& mod, const arma::mat& modn, const arma::umat& B);
RcppExport SEXP Rcppsva_bootstrap_regress(SEXP MSEXP, SEXP modSEXP, SEXP modnSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mod(modSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type modn(modnSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type B(BSEXP);
    __result = Rcpp::wrap(bootstrap_regress(M, mod, modn, B));
    return __result;
END_RCPP
}
// single_linkage
arma::mat single_linkage(arma::mat& M);
RcppExport SEXP Rcppsva_single_linkage(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    __result = Rcpp::wrap(single_linkage(M));
    return __result;
END_RCPP
}
// complete_linkage
arma::mat complete_linkage(arma::mat& M);
RcppExport SEXP Rcppsva_complete_linkage(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    __result = Rcpp::wrap(complete_linkage(M));
    return __result;
END_RCPP
}
// average_linkage
arma::mat average_linkage(arma::mat& M);
RcppExport SEXP Rcppsva_average_linkage(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    __result = Rcpp::wrap(average_linkage(M));
    return __result;
END_RCPP
}
// clique_merge
Rcpp::NumericVector clique_merge(arma::mat M);
RcppExport SEXP Rcppsva_clique_merge(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    __result = Rcpp::wrap(clique_merge(M));
    return __result;
END_RCPP
}
// linkage_kinds
Rcpp::CharacterVector linkage_kinds();
RcppExport SEXP Rcppsva_linkage_kinds() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(linkage_kinds());
    return __result;
END_RCPP
}
// distance_kinds
Rcpp::CharacterVector distance_kinds();
RcppExport SEXP Rcppsva_distance_kinds() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(distance_kinds());
    return __result;
END_RCPP
}
// rankm
arma::mat rankm(const arma::mat& M, SEXP byrow);
RcppExport SEXP Rcppsva_rankm(SEXP MSEXP, SEXP byrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< SEXP >::type byrow(byrowSEXP);
    __result = Rcpp::wrap(rankm(M, byrow));
    return __result;
END_RCPP
}
