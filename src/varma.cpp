// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utils.h"
using namespace arma;
using namespace Rcpp;


//VARMA MSE estimation
// [[Rcpp::export]]
List VARMA_MSE(arma::mat B,     // companion form of transition matrix
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
    
  List Out;
  Out["MSE"]   = MSE;
  Out["YP"] = YP;
  Out["E"] = E;
    
  return(Out);
}
  
//Disturbance smoothing --- output is only smoothed factors for simulations
// [[Rcpp::export]]
double KLike(arma::mat B,     // companion form of transition matrix
                arma::mat q,     // covariance matrix of shocks to states
                arma::mat H,     // measurement equation
                arma::vec R,     // covariance matrix of shocks to observables; Y are observations
                arma::mat Y){    // mixed frequency weights
  
  
  // preliminaries
  uword T  = Y.n_rows; //number of time peridos
  uword m  = B.n_rows; //number of factors
  uword p  = B.n_cols/m; //number of lags (must agree with lev/diff structure of data)
  uword k  = H.n_rows; //number of observables
  uword sA = m*p; //size of companion matrix A
  uword sH = H.n_cols;
  
  // Helper matrix
  mat HJ(k,sA, fill::zeros);
  HJ.cols(0,sH-1) = H;
  
  //Making the A matrix
  sp_mat Bs(B);
  sp_mat tmp_sp(sA-m,m);
  tmp_sp = join_horiz(speye<sp_mat>(sA-m,sA-m), tmp_sp);
  sp_mat A = join_vert(Bs, tmp_sp);
  
  //Making the Q matrix
  mat qq(sA,sA,fill::zeros);
  qq(span(0,m-1),span(0,m-1)) = q;
  sp_mat Q(qq);
  
  //Using the long run variance
  mat P0, P1, S, C;
  mat XX(eye<mat>(sA*sA, sA*sA) - kron(A,A));
  mat vQ(reshape(Q, sA*sA, 1));
  mat tmp_P = solve(XX, vQ);
  mat Pi(reshape(tmp_P, sA, sA));
  P0  = Pi; //long run variance
  P1  = P0;
  
  //Declairing variables for the filter
  mat VarY, Z(T,sA,fill::zeros),
  Lik, K, Si, tmp_mat, Ytmp, Hn;
  sp_mat Rn;
  vec PE, Yt, Yn, Yp, Wn;
  uvec ind, indM;
  double tmp;
  double tmpp;
  Lik << 0;
  vec Zp(sA, fill::zeros); //initialize to zero, the long run EV
  
  // -------- Filtering --------------------
  for(uword t=0; t<T; t++) {
    // Rcpp::Rcout << t << endl;
    //Doing the necessary aggregations and allowing for missing Y values
    Yn = trans(Y.row(t)); //returned as a col vec
    ind = find_finite(Yn);
    Yn = Yn(ind);
    // if nothing is observed
    if(Yn.is_empty()){
      Z.row(t) = trans(Zp);
      P0       = P1;
    } else{
      //if variables are observed
      Hn        = HJ.rows(ind); //rows of HJ corresponding to observations
      Rn        = diagmat(R(ind));  //rows of R corresponding to observations
      Yp        = Hn*Zp; //prediction step for Y
      S         = Hn*P1*trans(Hn)+Rn; //variance of Yp
      S         = symmatu((S+trans(S))/2); //enforce pos. semi. def.
      Si        = inv_sympd(S); //invert S
      //K         = trans(solve(S,Hn*P1));
      K         = P1*trans(Hn)*Si; //Kalman gain
      PE        = Yn-Yp; // prediction error
      Z.row(t)  = trans(Zp+K*PE); //updating step for Z
      P0        = P1-P1*trans(Hn)*Si*Hn*P1; // variance Z(t+1)|Y(1:t+1)
      //P0        = P1-P1*trans(Hn)*trans(K); // variance Z(t+1)|Y(1:t+1)
      P0        = symmatu((P0+trans(P0))/2); //enforce pos semi def
      log_det(tmp,tmpp,S); //calculate log determinant of S for the likelihood
      Lik    = -.5*tmp-.5*trans(PE)*Si*PE+Lik; //calculate log likelihood
      //MSE = trans(PE)*PE + MSE; //actual log likelihood gets stuck in local maxima
    }
    // Prediction for next period
    Zp  = A*trans(Z.row(t)); //prediction for Z(t+1) +itcZ
    P1     = A*P0*trans(A)+Q; //variance Z(t+1)|Y(1:t)
    P1     = symmatu((P1+trans(P1))/2); //enforce pos semi def
  }
  return(as_scalar(Lik));
}

// Disturbance smoother. Output is a list.
// [[Rcpp::export]]
List DKsmooth(arma::mat B,     // companion form of transition matrix
              arma::mat q,     // covariance matrix of shocks to states
              arma::mat H,     // measurement equation
              arma::vec R,     // covariance matrix of shocks to observables; Y are observations
              arma::mat Y){    // data
  
  
  // preliminaries
  uword T  = Y.n_rows; //number of time peridos
  uword m  = B.n_rows; //number of factors
  uword p  = B.n_cols/m; //number of lags (must agree with lev/diff structure of data)
  uword k  = H.n_rows; //number of observables
  uword sA = m*p; //size of companion matrix A
  uword sH = H.n_cols;
  
  // Helper matrix
  mat HJ(k,sA, fill::zeros);
  HJ.cols(0,sH-1) = H;
  
  //Making the A matrix
  sp_mat Bs(B);
  sp_mat tmp_sp(sA-m,m);
  tmp_sp = join_horiz(speye<sp_mat>(sA-m,sA-m), tmp_sp);
  sp_mat A  = join_vert(Bs, tmp_sp);
  
  //Making the Q matrix
  mat qq(sA,sA,fill::zeros);
  qq(span(0,m-1),span(0,m-1)) = q;
  sp_mat Q(qq);
  
  //Using the long run variance
  mat P0, P1, S, C;
  mat XX(eye<mat>(sA*sA, sA*sA) - kron(A,A));
  mat vQ(reshape(Q, sA*sA, 1));
  mat tmp_P = solve(XX, vQ);
  mat Pi(reshape(tmp_P, sA, sA));
  P0  = Pi; //long run variance
  P1  = P0;
  
  //Declairing variables for the filter
  field<mat> Kstr(T); //store Kalman gain
  field<vec> PEstr(T); //store prediction error
  field<mat> Hstr(T), Sstr(T); //store H and S^-1
  field<vec> Cmeans;
  mat VarY, ZP(T+1,sA,fill::zeros), Z(T,sA,fill::zeros), Zs(T,sA,fill::zeros),
  Lik, K, Si, tmp_mat, Ytmp, Hn;
  sp_mat Rn;
  vec PE, Yt, Yn, Yp;
  uvec ind, indM;
  double tmp;
  double tmpp;
  Lik << 0;
  vec Zp(sA, fill::zeros); //initialize to zero, the long run EV
  
  mat zippo(1,1,fill::zeros);
  mat zippo_sA(sA,1);
  
  
  
  // -------- Filtering --------------------
  for(uword t=0; t<T; t++) {
    //Rcpp::Rcout << t << endl;
    //Using W to for mixed frequency weights
    Yn = trans(Y.row(t));
    ind = find_finite(Yn);
    Yn = Yn(ind);
    // if nothing is observed
    if(Yn.is_empty()){
      Z.row(t) = trans(Zp);
      P0       = P1;
      Hstr(t)  = trans(zippo_sA);
      Sstr(t)  = zippo;
      PEstr(t) = zippo;
      Kstr(t)  = zippo_sA;
    } else{
      //if variables are observed
      Hn        = HJ.rows(ind); //rows of HJ corresponding to observations
      Hstr(t)   = Hn; //Store for smoothing
      Rn        = diagmat(R(ind));  //rows of R corresponding to observations
      Yp        = Hn*Zp; //prediction step for Y
      S         = Hn*P1*trans(Hn)+Rn; //variance of Yp
      S         = symmatu((S+trans(S))/2); //enforce pos. semi. def.
      Si        = inv_sympd(S); //invert S
      Sstr(t)   = Si; //sotre Si for smoothing
      K         = P1*trans(Hn)*Si; //Kalman gain
      PE        = Yn-Yp; // prediction error
      PEstr(t)  = PE; //store prediction error
      Kstr(t)   = K;  //store Kalman gain
      Z.row(t)  = trans(Zp+K*PE); //updating step for Z
      P0        = P1-P1*trans(Hn)*Si*Hn*P1; // variance Z(t+1)|Y(1:t+1)
      P0        = symmatu((P0+trans(P0))/2); //enforce pos semi def
      log_det(tmp,tmpp,S); //calculate log determinant of S for the likelihood
      Lik    = -.5*tmp-.5*trans(PE)*Si*PE+Lik; //calculate log likelihood
    }
    // Prediction for next period
    Zp  = A*trans(Z.row(t)); //prediction for Z(t+1) +itcZ
    ZP.row(t+1) = trans(Zp);
    P1     = A*P0*trans(A)+Q; //variance Z(t+1)|Y(1:t)
    P1     = symmatu((P1+trans(P1))/2); //enforce pos semi def
  }
  ZP.shed_row(T);
  
  //Smoothing following Durbin Koopman 2001/2012
  mat r(T+1,sA,fill::zeros);
  mat L;
  
  //r is 1 indexed while all other variables are zero indexed
  for(uword t=T; t>0; t--) {
    L     = (A-A*Kstr(t-1)*Hstr(t-1));
    r.row(t-1) = trans(PEstr(t-1))*Sstr(t-1)*Hstr(t-1) + r.row(t)*L;
    // Rcpp::Rcout << t << endl;
  }
  
  Zs.row(0)   = r.row(0)*Pi;
  
  //Forward again
  for(uword t = 0; t<T-1; t++){
    Zs.row(t+1)   = Zs.row(t)*trans(A) + r.row(t+1)*Q; //smoothed values of Z
  }
  
  mat Ys = Zs.cols(0,sH-1)*trans(H); //fitted values of Y
  
  List Out;
  Out["Ys"]   = Ys;
  Out["Lik"]  = Lik;
  // Out["Zz"]   = Z;
  Out["Z"]    = Zs;
  Out["Zp"]   = ZP;
  Out["Kstr"] = Kstr;
  Out["PEstr"]= PEstr;
  Out["A"] = A;
  
  return(Out);
}

// Disturbance smoother. Output is a list.
// [[Rcpp::export]]
List DKsmoothMF(arma::mat B,     // companion form of transition matrix
              arma::mat q,     // covariance matrix of shocks to states
              arma::mat H,     // measurement equation
              arma::vec R,     // covariance matrix of shocks to observables; Y are observations
              arma::mat Y,
              arma::mat W){    // data
  
  
  // preliminaries
  uword T  = Y.n_rows; //number of time peridos
  uword m  = B.n_rows; //number of factors
  uword p  = B.n_cols/m; //number of lags (must agree with lev/diff structure of data)
  uword k  = H.n_rows; //number of observables
  uword sA = m*p; //size of companion matrix A
  uword sH = H.n_cols;
  
  // Helper matrix
  mat HJ(k,sA, fill::zeros);
  HJ.cols(0,sH-1) = H;
  
  //Making the A matrix
  sp_mat Bs(B);
  sp_mat tmp_sp(sA-m,m);
  tmp_sp = join_horiz(speye<sp_mat>(sA-m,sA-m), tmp_sp);
  sp_mat A  = join_vert(Bs, tmp_sp);
  
  //Making the Q matrix
  mat qq(sA,sA,fill::zeros);
  qq(span(0,m-1),span(0,m-1)) = q;
  sp_mat Q(qq);
  
  //Using the long run variance
  mat P0, P1, S, C;
  mat XX(eye<mat>(sA*sA, sA*sA) - kron(A,A));
  mat vQ(reshape(Q, sA*sA, 1));
  mat tmp_P = solve(XX, vQ);
  mat Pi(reshape(tmp_P, sA, sA));
  P0  = Pi; //long run variance
  P1  = P0;
  
  //Declairing variables for the filter
  field<mat> Kstr(T); //store Kalman gain
  field<vec> PEstr(T); //store prediction error
  field<mat> Hstr(T), Sstr(T); //store H and S^-1
  field<vec> Cmeans;
  mat VarY, ZP(T+1,sA,fill::zeros), Z(T,sA,fill::zeros), Zs(T,sA,fill::zeros),
  Lik, K, Si, tmp_mat, Ytmp, Hn;
  sp_mat Rn;
  vec PE, Yt, Yn, Wn, Yp;
  uvec ind, indM;
  double tmp;
  double tmpp;
  Lik << 0;
  vec Zp(sA, fill::zeros); //initialize to zero, the long run EV
  
  mat zippo(1,1,fill::zeros);
  mat zippo_sA(sA,1);
  
  
  
  // -------- Filtering --------------------
  for(uword t=0; t<T; t++) {
    //Rcpp::Rcout << t << endl;
    //Using W to for mixed frequency weights
    Yn = trans(Y.row(t));
    ind = find_finite(Yn);
    Yn = Yn(ind);
    Wn = trans(W.row(t));
    Wn = Wn(ind);
    // if nothing is observed
    if(Yn.is_empty()){
      Z.row(t) = trans(Zp);
      P0       = P1;
      Hstr(t)  = trans(zippo_sA);
      Sstr(t)  = zippo;
      PEstr(t) = zippo;
      Kstr(t)  = zippo_sA;
    } else{
      //if variables are observed
      Hn        = HJ.rows(ind); //rows of HJ corresponding to observations
      Hstr(t)   = Hn; //Store for smoothing
      Rn        = diagmat(R(ind)%Wn);  //rows of R corresponding to observations
      Yp        = Hn*Zp; //prediction step for Y
      S         = Hn*P1*trans(Hn)+Rn; //variance of Yp
      S         = symmatu((S+trans(S))/2); //enforce pos. semi. def.
      Si        = inv_sympd(S); //invert S
      Sstr(t)   = Si; //sotre Si for smoothing
      K         = P1*trans(Hn)*Si; //Kalman gain
      PE        = Yn-Yp; // prediction error
      PEstr(t)  = PE; //store prediction error
      Kstr(t)   = K;  //store Kalman gain
      Z.row(t)  = trans(Zp+K*PE); //updating step for Z
      P0        = P1-P1*trans(Hn)*Si*Hn*P1; // variance Z(t+1)|Y(1:t+1)
      P0        = symmatu((P0+trans(P0))/2); //enforce pos semi def
      log_det(tmp,tmpp,S); //calculate log determinant of S for the likelihood
      Lik    = -.5*tmp-.5*trans(PE)*Si*PE+Lik; //calculate log likelihood
    }
    // Prediction for next period
    Zp  = A*trans(Z.row(t)); //prediction for Z(t+1) +itcZ
    ZP.row(t+1) = trans(Zp);
    P1     = A*P0*trans(A)+Q; //variance Z(t+1)|Y(1:t)
    P1     = symmatu((P1+trans(P1))/2); //enforce pos semi def
  }
  ZP.shed_row(T);
  
  //Smoothing following Durbin Koopman 2001/2012
  mat r(T+1,sA,fill::zeros);
  mat L;
  
  //r is 1 indexed while all other variables are zero indexed
  for(uword t=T; t>0; t--) {
    L     = (A-A*Kstr(t-1)*Hstr(t-1));
    r.row(t-1) = trans(PEstr(t-1))*Sstr(t-1)*Hstr(t-1) + r.row(t)*L;
    // Rcpp::Rcout << t << endl;
  }
  
  Zs.row(0)   = r.row(0)*Pi;
  
  //Forward again
  for(uword t = 0; t<T-1; t++){
    Zs.row(t+1)   = Zs.row(t)*trans(A) + r.row(t+1)*Q; //smoothed values of Z
  }
  
  mat Ys = Zs.cols(0,sH-1)*trans(H); //fitted values of Y
  
  List Out;
  Out["Ys"]   = Ys;
  Out["Lik"]  = Lik;
  // Out["Zz"]   = Z;
  Out["Z"]    = Zs;
  Out["Zp"]   = ZP;
  Out["Kstr"] = Kstr;
  Out["PEstr"]= PEstr;
  Out["A"] = A;
  
  return(Out);
}

// // [[Rcpp::export]]
// List VARMA(arma::mat Y, //multivariate input data
//            arma::mat P, //AR parameters
//            arma::mat Q){ //MA parameters
//   
//   uword T   = Y.n_rows;
//   uword sp  = P.n_cols;
//   uword sq  = Q.n_cols;
//   uword pq  = std::max(sp,sq);
//   double s_AR, s_MA;
//   vec E(T, fill::zeros);
//   vec eps(T, fill::zeros);
//   vec seas(T_long, fill::zeros);
//   vec YP(T, fill::zeros);
//   vec Y_sa(T, fill::zeros);
//   uvec all_ind;
//   uvec Pl(2, fill::zeros);
//   uvec Ql(2, fill::zeros);
//   
//   Z = stack_obs(Y,sp);
// 
//   if(sp==0 && sq==0){
//     Rcpp::stop("Nothing to estimate");
//   }else if(sp>0 && sq==0){
//     for(uword t = srt; t<T; t++){
//       if( !Y(span(t-sp,t)).is_finite() ) continue;
//       YP(t) = as_scalar(P*flipud(Y(span(t-sp,t-1)))); //AR only
//       E(t) = Y(t) - YP(t);
//     }
//     
//   }else if(sp==0 && sq>0){
//     Rcpp::stop("Seasonal adjustment requires at least 1 AR lag");
//   }else{
//     if(P.n_cols==0 && Q.n_cols==0){ //q not used
//       Rcpp::stop("Seasonal adjustment with q>0 requires Q>0");
//     }else if(P.n_cols>0 && Q.n_cols==0){ //q not used
//       Rcpp::stop("Seasonal adjustment with q>0 requires Q>0");
//     }else if(P.n_cols==0 && Q.n_cols>0){
//       for(uword t = srt; t<T; t++){
//         s_MA = as_scalar(Q*eps(Q_lag.row(t))); //seasonal MA covariates
//         if( !Y(span(t-pq,t)).is_finite() || !is_finite(s_MA) ) continue;
//         seas(t) = s_MA; //seasonal component
//         Y_sa(t) = Y(t) - seas(t);
//         YP(t) = as_scalar(p*flipud(Y_sa(span(t-sp,t-1)))); //adding MA components
//         E(t) = Y_sa(t) - YP(t);
//         eps(t) = as_scalar(q*flipud(eps(span(t-sq,t-1)))) + E(t);
//       }
//     }else{
//       for(uword t = srt; t<T; t++){
//         s_MA = as_scalar(Q*eps(Q_lag.row(t))); //seasonal MA covariates
//         s_AR = as_scalar(P*Y(P_lag.row(t))); //seasonal AR covariates
//         if( !Y(span(t-pq,t)).is_finite() || !is_finite(s_MA) || !is_finite(s_AR) ) continue;
//         seas(t) = s_AR + s_MA; //seasonal component
//         Y_sa(t) = Y(t) - seas(t);
//         YP(t) = as_scalar(p*flipud(Y_sa(span(t-sp,t-1)))); //predicting next period SA values
//         E(t) = Y_sa(t) - YP(t);
//         eps(t) = as_scalar(q*flipud(eps(span(t-sq,t-1)))) + E(t);
//       }
//     }
//   }
//   
//   //Get future seasonal adjustments if required
//   if(T_long > T){
//     if(P.n_cols==0 && Q.n_cols>0){
//       for(uword t = T; t<T_long; t++){
//         s_MA = as_scalar(Q*eps(Q_lag.row(t))); //seasonal MA covariates
//         if( !is_finite(s_MA) ) continue;
//         seas(t) = s_MA; //seasonal component
//       }
//     }else if(P.n_cols>0 && Q.n_cols==0){
//       for(uword t = T; t<T_long; t++){
//         s_AR = as_scalar(P*Y(P_lag.row(t))); //seasonal AR covariates
//         if( !is_finite(s_AR) ) continue;
//         seas(t) = s_AR; //seasonal component
//       }
//     }else{
//       for(uword t = T; t<T_long; t++){
//         s_MA = as_scalar(Q*eps(Q_lag.row(t))); //seasonal MA covariates
//         s_AR = as_scalar(P*Y(P_lag.row(t))); //seasonal AR covariates
//         if( !is_finite(s_MA) || !is_finite(s_AR) ) continue;
//         seas(t) = s_AR + s_MA; //seasonal component
//       }
//     }
//   } // if(T_long > T)
//   
//   //double MSE = as_scalar(trans(E(span(srt,T-1)))*E(span(srt,T-1)))/T;
//   
//   uvec ind = find(E); //find non-zero elements of E
//   double MSE = as_scalar(trans(E(ind))*E(ind))/(ind.n_elem);
//   
//   List Out;
//   Out["MSE"] = MSE;
//   Out["seas"]  = seas;
//   Out["Y_sa"] = Y_sa;
//   Out["E"] = E;
//   Out["YP"] = YP;
//   Out["p"]  = p;
//   Out["q"] = q;
//   Out["P"] = P;
//   Out["Q"] = Q;
//   Out["srt"] = srt;
//   
//   return(Out);
// }

