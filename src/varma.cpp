// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utils.h"
using namespace arma;
using namespace Rcpp;


//VARMA MSE estimation
// [[Rcpp::export]]
Rcpp::List VARMA_MSE(arma::mat B,     // companion form of transition matrix
                 arma::mat Q,     // measurement equation
                 arma::mat Y){    // mixed frequency weights
  
  
  // preliminaries
  uword T  = Y.n_rows; //number of time peridos
  uword m  = Y.n_cols; //number of factors
  uword p  = B.n_cols/m; //number of lags (must agree with lev/diff structure of data)
  uword q  = Q.n_cols/m; //number of observables
  
  //Declairing variables for the filter
  mat YP(T, m), E(T,m,fill::zeros);
  YP.fill(datum::nan);
  vec z, zp(m*p, fill::zeros), e, ep(m*q, fill::zeros), yp;
  uvec ind;
  
  // -------- Iterating --------------------
  if(p>0 && q==0){
    for(uword t=p; t<T; t++) {
      // Rcpp::Rcout << t << endl;
      z = vectorise(trans(flipud(Y.rows(t-p,t-1))));
      ind = find_nonfinite(z);
      z(ind) = zp(ind);
      yp = B*z;
      YP.row(t) = trans(yp);
      e = trans(Y.row(t)) - yp;
      ind = find_nonfinite(e);
      e(ind) = zeros<vec>(ind.n_elem);
      E.row(t) = trans(e);
      zp = advance_vec(yp,z);
    }
  }else if(p==0 && q>0){
    for(uword t=0; t<T; t++) {
      //Rcpp::Rcout << t << endl;
      yp = Q*ep;
      YP.row(t) = trans(yp);
      e = trans(Y.row(t)) - yp;
      ind = find_nonfinite(e);
      e(ind) = zeros<vec>(ind.n_elem);
      E.row(t) = trans(e);
      ep = advance_vec(e,ep);
    }
  }else{
    for(uword t=p; t<T; t++) {
      // Rcpp::Rcout << "p" << endl;
      z = vectorise(trans(flipud(Y.rows(t-p,t-1))));
      ind = find_nonfinite(z);
      z(ind) = zp(ind);
      yp = B*z + Q*ep;
      YP.row(t) = trans(yp);
      e = trans(Y.row(t)) - yp;
      ind = find_nonfinite(e);
      e(ind) = zeros<vec>(ind.n_elem);
      E.row(t) = trans(e);
      zp = advance_vec(yp,z);
      ep = advance_vec(e,ep);
    }
  }

  double MSE = trace(trans(E)*E);
    
  Rcpp::List Out;
  Out["MSE"]   = MSE;
  Out["YP"] = YP;
  Out["E"] = E;
    
  return(Out);
}

// Lightning fast Uniform frequency filtering
// [[Rcpp::export]]
Rcpp::List Kfilter(arma::mat Y,     // Observations Y
                   arma::mat B,  // transition matrix
                   arma::mat q,  // shocks to factors
                   arma::mat H,
                   arma::vec R){     // shocks to observations

  // preliminaries
  uword T  = Y.n_cols; //number of time peridos
  uword m  = B.n_rows; //number of factors
  uword p  = B.n_cols/m; //number of lags 
  uword sA = m*p; //size of companion matrix A
  
  //Making the A matrix
  sp_mat Bs(B);
  sp_mat tmp_sp(sA-m,m);
  tmp_sp = join_horiz(speye<sp_mat>(sA-m,sA-m), tmp_sp);
  sp_mat A  = join_vert(Bs, tmp_sp);
  
  //Making the Q matrix
  mat qq(sA,sA,fill::zeros);
  qq(span(0,m-1),span(0,m-1)) = q;
  sp_mat Q(qq);
  
  mat P0 = long_run_var(A,Q,m,p);
  mat P1 = P0;
  cube P0str(sA,sA,T+1); P0str.slice(0) = P0;
  mat VarY, ZP(sA,T,fill::zeros), Z(sA,T+1,fill::zeros),  
  K;
  uword i;
  double rn; double yn; double s;
  rowvec hn; double pe; double yp;
  double Lik;
  Lik = 0;
  sp_mat Hn;
  vec PE, Yt, Yn, Yp, Zp;
  uvec ind, indM, ix;
  
  // -------- Filtering --------------------
  for(uword t=0; t!=T; t++) {
    // Prediction for next period
    Zp  = A*Z.col(t); //prediction for Z(t+1) +itcZ
    ZP.col(t) = Zp;
    P1     = A*P0*trans(A)+Q; //variance Z(t+1)|Y(1:t)
    P1     = symmatu((P1+trans(P1))/2); //enforce pos semi def
    // P1str.slice(t) = P1;
    //Allowing for missing Y values
    Yt     = Y.col(t);
    ind    = find_finite(Yt);
    if(ind.is_empty()) {  // if nothing is observed
      Z.col(t+1) = Zp;
      P0       = P1;
      P0str.slice(t+1) = P0;
    } else { //if variables are observed
      for(uword j=0; j != ind.n_elem; j++){
        i = ind(j);
        //Hn        = HJ.rows(ind);
        yn = Yt(i);
        hn = H.row(i); //rows of HJ corresponding to observations
        rn = R(i);  //element of R corresponding to observations
        yp = as_scalar(hn*Zp); //prediction step for Y
        s  = (as_scalar(hn*P1*trans(hn)) + rn); //variance of Yp
        K  = P1*trans(hn)/s;   //trans(solve(S, Hn*Q)); //Kalman gain
        pe = yn-yp;
        Zp = Zp+K*pe; //updating step for Z
        P1 = P1-K*hn*P1; // variance Z(t+1)|Y(1:t+1)
        Lik = Lik - (log(2*datum::pi) + log(s) + pe*pe/s)/2;
      }
      Z.col(t+1) = Zp;
      P0 = symmatu((P1+trans(P1))/2); //enforce pos semi def
      P0str.slice(t+1) = P0;
    }
  }
  Z.shed_col(0);
  P0str.shed_slice(0);
  
  
  Rcpp::List rtrn;
  rtrn["Lik"]  = Lik;
  rtrn["Z"]  = Z; //includes pre-sample value in position 0
  rtrn["var"]  = P0str; // includes pre-sample value in position 0
  rtrn["Yfit"] = H*Z;
  //rtrn(3)  = PS_lag;
  return(rtrn);
}