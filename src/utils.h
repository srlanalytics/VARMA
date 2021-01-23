#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
// List BReg(arma::mat X, arma::mat Y, bool Int, arma::mat Bp, double lam, double nu,
//           arma::uword reps = 1000, arma::uword burn = 1000);
// List BReg_diag(arma::mat X,  arma::mat Y, bool Int, arma::mat Bp, arma::vec lam, arma::vec Y_lam, arma::vec nu,
//                arma::uword reps = 1000, arma::uword burn = 1000); 
// arma::sp_mat MakeSparse(arma::mat A);
// arma::sp_mat sp_rows(arma::sp_mat A, arma::uvec r);
// arma::sp_mat sp_cols(arma::sp_mat A, arma::uvec r);
// arma::sp_mat sprow(arma::sp_mat A, arma::mat a, arma::uword r);
arma::mat comp_form(arma::mat B);
// arma::mat mvrnrm(int n, arma::vec mu, arma::mat Sigma);
// arma::cube rinvwish(int n, int v, arma::mat S);
// double invchisq(double nu, double scale);
arma:: mat stack_obs(arma::mat nn, arma::uword p, arma::uword r = 0);
arma::vec advance_vec(arma::vec e, arma::vec E);
arma::mat long_run_var(arma::sp_mat A,
                       arma::sp_mat Q,
                       arma::uword m,
                       arma::uword p);

#endif