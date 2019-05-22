// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat mvrnrm(int n, arma::vec mu, arma::mat Sigma){
  /*-------------------------------------------------------
# Generate draws from a multivariate normal distribution
#--------------------------------------------------------
#  n        number of samples
#  mu       mean vector
#  Sigma    covariance matrix
#-------------------------------------------------------*/
   RNGScope scope;
  int p = Sigma.n_cols;
  mat X = reshape(vec(rnorm(p * n)), p, n);
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, Sigma);
  X = eigvec * diagmat(sqrt(eigval)) * X;
  X.each_col() += mu;
  return(X);
}

/*-------------------------------------------------------
# Generate Draws from an Inverse Wishart Distribution
# via the Bartlett Decomposition
#--------------------------------------------------------
# NOTE: output is identical to riwish from MCMCpack
#       provided the same random seed is used
#--------------------------------------------------------
#   n     number of samples
#   S     scale matrix
#   v     degrees of freedom
#-------------------------------------------------------*/
 // [[Rcpp::export]]
 arma::cube rinvwish(int n, int v, arma::mat S){
   RNGScope scope;
   int p = S.n_rows;
   mat L = chol(inv_sympd(S), "lower");
   cube sims(p, p, n, fill::zeros);
   for(int j = 0; j < n; j++){
     mat A(p,p, fill::zeros);
     for(int i = 0; i < p; i++){
       int df = v - (i + 1) + 1; //zero-indexing
       A(i,i) = sqrt(R::rchisq(df));
     }
     for(int row = 1; row < p; row++){
       for(int col = 0; col < row; col++){
         A(row, col) = R::rnorm(0,1);
       }
     }
     mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
     sims.slice(j) = LA_inv.t() * LA_inv;
   }
   return(sims);
 }

// [[Rcpp::export]]
double invchisq(double nu, double scale){
  /*-------------------------------------------------------
# Generate draws from a scaled inverse chi squared distribution
#--------------------------------------------------------
#  nu       "degrees of freedom"
#  scale    scale parameter
#-------------------------------------------------------*/
   vec    x = randn<vec>(nu)/sqrt(scale);
  double s = 1/sum(square(x));
  return(s);
}


//Bayesian linear regression --- does NOT accept missing variables
// [[Rcpp::export]]
List BReg(arma::mat X,   // RHS variables
          arma::mat Y,   // LHS variables
          bool Int,      // Estimate intercept term?
          arma::mat Bp,  // prior for B
          double lam,    // prior tightness
          double nu,     //prior "deg of freedom"
          arma::uword reps = 1000,
          arma::uword burn = 1000){
  
  uword k    = Y.n_cols;
  uword m    = X.n_cols;
  uword T    = Y.n_rows;
  mat Lam    = lam*eye<mat>(m,m);
  mat tmp;
  
  if(Int){ //if we should estimate an intercept term
    m        = m+1;
    Lam      = lam*eye<mat>(m,m);
    Lam(0,0) = 0;
    tmp      = zeros<mat>(k,1);
    Bp       = join_horiz(tmp, Bp);
    tmp      = ones<mat>(T,1);
    X        = join_horiz(tmp,X);
  }
  
  //declairing variables
  cube Bstore(k,m,reps), Qstore(k,k,reps);
  mat v_1, V_1, Mu, B, Beta, scale, q;
  vec mu;
  
  //Burn Loop
  
  for(uword rep = 0; rep<burn; rep++){
    
    Rcpp::checkUserInterrupt();
    
    v_1   = trans(X)*X+Lam;
    v_1   = (trans(v_1)+v_1)/2;
    v_1   = inv_sympd(v_1);
    Mu    = v_1*(trans(X)*Y+Lam*trans(Bp));
    scale = eye(k,k)+trans(Y-X*Mu)*(Y-X*Mu)+trans(Mu-trans(Bp))*Lam*(Mu-trans(Bp)); // eye(k) is the prior scale parameter for the IW distribution and eye(k)+junk the posterior.
    scale = (scale+trans(scale))/2;
    q     = rinvwish(1,nu+T,scale); // Draw for q
    mu    = vectorise(trans(Mu));   // vectorized posterior mean for beta
    V_1   = kron(v_1,q);            // covariance of vectorized parameters beta
    Beta  = mvrnrm(1,mu,V_1);       //Draw for B
    B     = reshape(Beta,k,m);      //recovering original dimensions
  }
  
  //Sampling Loop
  
  for(uword rep = 0; rep<reps; rep++){
    
    Rcpp::checkUserInterrupt();
    
    v_1   = trans(X)*X+Lam;
    v_1   = (trans(v_1)+v_1)/2;
    v_1   = inv_sympd(v_1);
    Mu    = v_1*(trans(X)*Y+Lam*trans(Bp));
    scale = eye(k,k)+trans(Y-X*Mu)*(Y-X*Mu)+trans(Mu-trans(Bp))*Lam*(Mu-trans(Bp)); // eye(k) is the prior scale parameter for the IW distribution and eye(k)+junk the posterior.
    scale = (scale+trans(scale))/2;
    q     = rinvwish(1,nu+T,scale); // Draw for q
    mu    = vectorise(trans(Mu));   // vectorized posterior mean for beta
    V_1   = kron(v_1,q);            // covariance of vectorized parameters beta
    Beta  = mvrnrm(1,mu,V_1);       //Draw for B
    B     = reshape(Beta,k,m);      //recovering original dimensions
    Qstore.slice(rep) = q;
    Bstore.slice(rep) = B;
  }
  
  //For B
  for(uword rw=0;rw<k;rw++){
    for(uword cl=0;cl<m;cl++){
      B(rw,cl) = as_scalar(median(vectorise(Bstore.tube(rw,cl))));
    }
  }
  //For q
  for(uword rw=0;rw<k;rw++){
    for(uword cl=0;cl<k;cl++){
      q(rw,cl) = as_scalar(median(vectorise(Qstore.tube(rw,cl))));
    }
  }
  
  List Out;
  Out["B"]   = B;
  Out["q"]   = q;
  Out["Bstore"]   = Bstore;
  Out["Qstore"]   = Qstore;
  
  return(Out);
  
}

//Bayesian linear regression with diagonal covariance to shocks. Accepts missing obs.
// [[Rcpp::export]]
List BReg_diag(arma::mat X,   // RHS variables
               arma::mat Y,   // LHS variables
               bool Int,      //Estimate intercept terms?
               arma::mat Bp,  // prior for B
               double lam,    // prior tightness
               arma::vec nu,  //prior "deg of freedom"
               arma::uword reps = 1000, //MCMC sampling iterations
               arma::uword burn = 1000){ //MCMC burn in iterations
  
  uword k    = Y.n_cols;
  uword m    = X.n_cols;
  uword T    = Y.n_rows;
  mat Lam    = lam*eye<mat>(m,m);
  mat tmp;
  
  if(Int){ //if we should estimate an intercept term
    m        = m+1;
    Lam      = lam*eye<mat>(m,m);
    Lam(0,0) = 0;
    tmp      = zeros<mat>(k,1);
    Bp       = join_horiz(tmp, Bp);
    tmp      = ones<mat>(T,1);
    X        = join_horiz(tmp,X);
  }
  
  //declairing variables
  cube Bstore(k,m,reps);
  mat v_1, Beta, Qstore(k,reps), xx;
  vec mu, q(k,fill::zeros);
  double scl;
  uvec ind_rm, ind;
  vec y, x;
  mat B(k, m, fill::zeros);
  uvec indX = find_nonfinite(X.col(0));
  if(m>1){
    for(uword j = 1; j<m; j++){
      x = X.col(j);
      indX = unique( join_cols(indX, find_nonfinite(x))); //index of missing X values
    }
  }
  
  //Burn Loop
  
  for(uword rep = 0; rep<burn; rep++){
    Rcpp::checkUserInterrupt();
    
    for(uword j=0; j<k; j++){
      y       = Y.col(j);
      ind     = unique(join_cols(indX,find_nonfinite(y))); //index of elements to remove
      xx      = X;
      
      //this seems a tedious way to shed non-contiguous indexes
      for(uword n = ind.n_elem; n>0; n--){
        xx.shed_row(ind(n-1));
        y.shed_row(ind(n-1));
      }
      
      v_1   = trans(xx)*xx+Lam;
      v_1   = (trans(v_1)+v_1)/2;
      v_1   = inv_sympd(v_1);
      mu    = v_1*(trans(xx)*y+Lam*trans(Bp.row(j)));
      scl   = as_scalar(trans(y-xx*mu)*(y-xx*mu)+trans(mu-trans(Bp.row(j)))*Lam*(mu-trans(Bp.row(j)))); // prior variance is zero... a little odd but it works
      q(j)  = invchisq(nu(j)+y.n_rows,scl); //Draw for r
      Beta  = mvrnrm(1, mu, v_1*q(j));
      B.row(j) = trans(Beta.col(0));
    }
    
  }
  
  // Sampling loop
  
  for(uword rep = 0; rep<reps; rep++){
    
    Rcpp::checkUserInterrupt();
    
    for(uword j=0; j<k; j++){
      y       = Y.col(j);
      ind     = unique(join_cols(indX,find_nonfinite(y))); //index of elements to remove
      xx      = X;
      //this seems a tedious way to shed non-contiguous indexes
      for(uword n = ind.n_elem; n>0; n--){
        xx.shed_row(ind(n-1));
        y.shed_row(ind(n-1));
      }
      v_1   = trans(xx)*xx+Lam;
      v_1   = (trans(v_1)+v_1)/2;
      v_1   = inv_sympd(v_1);
      mu    = v_1*(trans(xx)*y+Lam*trans(Bp.row(j)));
      scl   = as_scalar(trans(y-xx*mu)*(y-xx*mu)+trans(mu-trans(Bp.row(j)))*Lam*(mu-trans(Bp.row(j)))); // prior variance is zero... a little odd but it works
      q(j)  = invchisq(nu(j)+y.n_rows,scl); //Draw for r
      Beta  = mvrnrm(1, mu, v_1*q(j));
      B.row(j) = trans(Beta.col(0));
    }
    Qstore.col(rep)   = q;
    Bstore.slice(rep) = B;
    
  }
  
  
  //For B
  for(uword rw=0;rw<k;rw++){
    for(uword cl=0;cl<m;cl++){
      B(rw,cl) = as_scalar(median(vectorise(Bstore.tube(rw,cl))));
    }
  }
  //For q
  
  for(uword rw=0;rw<k;rw++){
    q(rw) = as_scalar(median(Qstore.row(rw)));
  }
  
  
  List Out;
  Out["B"]   = B;
  Out["q"]   = q;
  Out["Bstore"]   = Bstore;
  Out["Qstore"]   = Qstore;
  
  return(Out);
}


  
arma::sp_mat MakeSparse(arma::mat A){
  uword n_rows   = A.n_rows;
  uword n_cols   = A.n_cols;
  uvec ind       = find(A);
  umat locations = ind2sub(size(A),ind);
  vec  values    = A(ind);
  sp_mat C(locations,values,n_rows,n_cols);
  return(C);
}

arma::sp_mat sp_rows(arma::sp_mat A,
                     arma::uvec r   ){
  uword n_rows   = A.n_rows;
  //  uword n_cols   = A.n_cols;
  uword n_r      = r.size();
  uvec  tmp      = regspace<uvec>(0,n_rows-1);
  tmp      = tmp.elem(r);
  umat  location = join_vert(trans(regspace<uvec>(0,n_r-1)),trans(tmp));
  sp_mat J(location,ones<vec>(n_r),n_r,n_rows);
  sp_mat C       = J*A;
  return(C);
}

arma::sp_mat sp_cols(arma::sp_mat A,
                     arma::uvec r   ){
  //  uword n_rows   = A.n_rows;
  uword n_cols   = A.n_cols;
  uword n_r      = r.size();
  uvec  tmp      = regspace<uvec>(0,n_cols-1);
  tmp            = tmp.elem(r);
  umat  location = join_vert(trans(tmp),trans(regspace<uvec>(0,n_r-1)));
  sp_mat J(location,ones<vec>(n_r),n_cols,n_r);
  sp_mat C       = A*J;
  return(C);
}


//Replace row r of sparse matrix A with the (sparse) vector a.
//Should be reasonably fast with Armadillo 8 or newer
arma::sp_mat sprow(arma::sp_mat A,
                   arma::mat a,
                   arma::uword r   ){
  //This intitally used find(a) to inentify non-zero elements of a, but that
  //did not replace elements that are non-zero in A and zero in a
  uword n_cols     = A.n_cols;
  if(n_cols>a.n_elem){
    a = join_horiz(a, zeros<mat>(1,n_cols-a.n_elem));
  }
  for(uword n      = 0; n < n_cols; n++){
    A(r,n)         = a(n);
  }
  return(A);
}

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

// // [[Rcpp::export]]
// arma::vec vtest(arma::mat A){
//   vec out = vectorise(trans(flipud(A)));
//   return(out);
// }

