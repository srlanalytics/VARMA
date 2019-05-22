// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mvrnrm
arma::mat mvrnrm(int n, arma::vec mu, arma::mat Sigma);
RcppExport SEXP _VARMA_mvrnrm(SEXP nSEXP, SEXP muSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnrm(n, mu, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// rinvwish
arma::cube rinvwish(int n, int v, arma::mat S);
RcppExport SEXP _VARMA_rinvwish(SEXP nSEXP, SEXP vSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(rinvwish(n, v, S));
    return rcpp_result_gen;
END_RCPP
}
// invchisq
double invchisq(double nu, double scale);
RcppExport SEXP _VARMA_invchisq(SEXP nuSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(invchisq(nu, scale));
    return rcpp_result_gen;
END_RCPP
}
// BReg
List BReg(arma::mat X, arma::mat Y, bool Int, arma::mat Bp, double lam, double nu, arma::uword reps, arma::uword burn);
RcppExport SEXP _VARMA_BReg(SEXP XSEXP, SEXP YSEXP, SEXP IntSEXP, SEXP BpSEXP, SEXP lamSEXP, SEXP nuSEXP, SEXP repsSEXP, SEXP burnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type Int(IntSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Bp(BpSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    rcpp_result_gen = Rcpp::wrap(BReg(X, Y, Int, Bp, lam, nu, reps, burn));
    return rcpp_result_gen;
END_RCPP
}
// BReg_diag
List BReg_diag(arma::mat X, arma::mat Y, bool Int, arma::mat Bp, double lam, arma::vec nu, arma::uword reps, arma::uword burn);
RcppExport SEXP _VARMA_BReg_diag(SEXP XSEXP, SEXP YSEXP, SEXP IntSEXP, SEXP BpSEXP, SEXP lamSEXP, SEXP nuSEXP, SEXP repsSEXP, SEXP burnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type Int(IntSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Bp(BpSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    rcpp_result_gen = Rcpp::wrap(BReg_diag(X, Y, Int, Bp, lam, nu, reps, burn));
    return rcpp_result_gen;
END_RCPP
}
// comp_form
arma::mat comp_form(arma::mat B);
RcppExport SEXP _VARMA_comp_form(SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(comp_form(B));
    return rcpp_result_gen;
END_RCPP
}
// stack_obs
arma:: mat stack_obs(arma::mat nn, arma::uword p, arma::uword r);
RcppExport SEXP _VARMA_stack_obs(SEXP nnSEXP, SEXP pSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(stack_obs(nn, p, r));
    return rcpp_result_gen;
END_RCPP
}
// advance_vec
arma::vec advance_vec(arma::vec e, arma::vec E);
RcppExport SEXP _VARMA_advance_vec(SEXP eSEXP, SEXP ESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type e(eSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type E(ESEXP);
    rcpp_result_gen = Rcpp::wrap(advance_vec(e, E));
    return rcpp_result_gen;
END_RCPP
}
// VARMA_MSE
List VARMA_MSE(arma::mat B, arma::mat Q, arma::mat Y);
RcppExport SEXP _VARMA_VARMA_MSE(SEXP BSEXP, SEXP QSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(VARMA_MSE(B, Q, Y));
    return rcpp_result_gen;
END_RCPP
}
// KLike
double KLike(arma::mat B, arma::mat q, arma::mat H, arma::vec R, arma::mat Y);
RcppExport SEXP _VARMA_KLike(SEXP BSEXP, SEXP qSEXP, SEXP HSEXP, SEXP RSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type q(qSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(KLike(B, q, H, R, Y));
    return rcpp_result_gen;
END_RCPP
}
// DKsmooth
List DKsmooth(arma::mat B, arma::mat q, arma::mat H, arma::vec R, arma::mat Y);
RcppExport SEXP _VARMA_DKsmooth(SEXP BSEXP, SEXP qSEXP, SEXP HSEXP, SEXP RSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type q(qSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(DKsmooth(B, q, H, R, Y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VARMA_mvrnrm", (DL_FUNC) &_VARMA_mvrnrm, 3},
    {"_VARMA_rinvwish", (DL_FUNC) &_VARMA_rinvwish, 3},
    {"_VARMA_invchisq", (DL_FUNC) &_VARMA_invchisq, 2},
    {"_VARMA_BReg", (DL_FUNC) &_VARMA_BReg, 8},
    {"_VARMA_BReg_diag", (DL_FUNC) &_VARMA_BReg_diag, 8},
    {"_VARMA_comp_form", (DL_FUNC) &_VARMA_comp_form, 1},
    {"_VARMA_stack_obs", (DL_FUNC) &_VARMA_stack_obs, 3},
    {"_VARMA_advance_vec", (DL_FUNC) &_VARMA_advance_vec, 2},
    {"_VARMA_VARMA_MSE", (DL_FUNC) &_VARMA_VARMA_MSE, 3},
    {"_VARMA_KLike", (DL_FUNC) &_VARMA_KLike, 5},
    {"_VARMA_DKsmooth", (DL_FUNC) &_VARMA_DKsmooth, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_VARMA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
