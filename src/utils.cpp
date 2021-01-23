// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;


// arma::sp_mat MakeSparse(arma::mat A){
//   uword n_rows   = A.n_rows;
//   uword n_cols   = A.n_cols;
//   uvec ind       = find(A);
//   umat locations = ind2sub(size(A),ind);
//   vec  values    = A(ind);
//   sp_mat C(locations,values,n_rows,n_cols);
//   return(C);
// }
// 
// arma::sp_mat sp_rows(arma::sp_mat A,
//                      arma::uvec r   ){
//   uword n_rows   = A.n_rows;
//   //  uword n_cols   = A.n_cols;
//   uword n_r      = r.size();
//   uvec  tmp      = regspace<uvec>(0,n_rows-1);
//   tmp      = tmp.elem(r);
//   umat  location = join_vert(trans(regspace<uvec>(0,n_r-1)),trans(tmp));
//   sp_mat J(location,ones<vec>(n_r),n_r,n_rows);
//   sp_mat C       = J*A;
//   return(C);
// }
// 
// arma::sp_mat sp_cols(arma::sp_mat A,
//                      arma::uvec r   ){
//   //  uword n_rows   = A.n_rows;
//   uword n_cols   = A.n_cols;
//   uword n_r      = r.size();
//   uvec  tmp      = regspace<uvec>(0,n_cols-1);
//   tmp            = tmp.elem(r);
//   umat  location = join_vert(trans(tmp),trans(regspace<uvec>(0,n_r-1)));
//   sp_mat J(location,ones<vec>(n_r),n_cols,n_r);
//   sp_mat C       = A*J;
//   return(C);
// }
// 
// 
// //Replace row r of sparse matrix A with the (sparse) vector a.
// //Should be reasonably fast with Armadillo 8 or newer
// arma::sp_mat sprow(arma::sp_mat A,
//                    arma::mat a,
//                    arma::uword r   ){
//   //This intitally used find(a) to inentify non-zero elements of a, but that
//   //did not replace elements that are non-zero in A and zero in a
//   uword n_cols     = A.n_cols;
//   if(n_cols>a.n_elem){
//     a = join_horiz(a, zeros<mat>(1,n_cols-a.n_elem));
//   }
//   for(uword n      = 0; n < n_cols; n++){
//     A(r,n)         = a(n);
//   }
//   return(A);
// }

//Create the companion form of the transition matrix B
// [[Rcpp::export]]
arma::mat comp_form(arma::mat B){
  uword r = B.n_rows;
  uword c = B.n_cols;
  mat A   = join_vert(B, join_horiz(eye<mat>(c-r,c-r), zeros<mat>(c-r,r)));
  return(A);
}

//Stack times series data in VAR format
// [[Rcpp::export]]
arma:: mat stack_obs(arma::mat nn, arma::uword p, arma::uword r = 0){
  uword rr = nn.n_rows;
  uword mn = nn.n_cols;
  if(r == 0){
    r = rr-p+1;
  }else if(rr-p+1 != r){
    Rcpp::stop("Length of input nn and length of data r do not agree.");
  }
  mat N(r,mn*p, fill::zeros);
  uword indx = 0;
  for(uword j = 1; j<=p; j++){
    N.cols(indx,indx+mn-1) = nn.rows(p-j,rr-j);
    indx = indx+mn;
  }
  return(N);
}


// [[Rcpp::export]]
arma::vec advance_vec(arma::vec e, arma::vec E){
  uword s = E.n_elem;
  uword m = e.n_elem;
  vec out(s);
  if(s == m){
    out = e;
  }else{
    out(span(0,m-1)) = e;
    out(span(m,s-1)) = E(span(0,s-m-1));
  }
  return(out);
}

// Getting the long run variance to initiate the filter
// [[Rcpp::export]]
arma::mat long_run_var(arma::sp_mat A,
                       arma::sp_mat Q,
                       arma::uword m,
                       arma::uword p){
  uword sA = A.n_cols;
  uword pp = sA/m;
  uword mp = m*p;
  double mp2 = mp*mp;
  mat B(A(span(0,mp-1), span(0,mp-1))); mat BB; mat b;
  mat XX(eye<mat>(mp2, mp2) - kron(B,B));
  vec vQ(reshape(Q(span(0,m*p-1),span(0,m*p-1)), mp2, 1));
  mat P(reshape(solve(XX, vQ), mp, mp));
  //P = P(span(0,m-1),span(0,m-1));
  mat PP(sA,sA,fill::zeros);
  for(uword j = 0; j<pp; j++){
    BB = B;
    for(uword k=j+1; k<pp; k++){
      b = BB*P;
      PP(span(m*j,m*j+m-1), span(m*k,m*k+m-1)) = b(span(0,m-1),span(0,m-1));
      BB = B*BB;
    }
  }
  PP = PP + trans(PP) + kron(eye<mat>(pp,pp), P(span(0,m-1), span(0,m-1)));
  return(PP);
}

// // [[Rcpp::export]]
// arma::vec vtest(arma::mat A){
//   vec out = vectorise(trans(flipud(A)));
//   return(out);
// }

