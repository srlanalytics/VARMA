get_from_list <- function(lst, what){
  if(is.character(what)){
    what <- which(names(lst)==what)
  }
  return(lst[[what]])
}

expand_mat <- function(a){
  r <- NROW(a)
  c <- NCOL(a)*r
  A <- matrix(0,r,c)
  idx <- seq(from = 1, to = c-r+1, by = r)
  for(j in seq(1,r)){
    A[j,idx] <- a[j,]
    idx <- idx+1
  }
  return(A)
}

breg <- function(X,Y, Int = TRUE, Bp = NULL, lam = 0, nu = NULL, reps = 2000, burn = 1000, diag = TRUE){
  
  k <- NCOL(X)
  m <- NCOL(Y)
  
  if(is.null(Bp)){
    Bp <- matrix(0,m,k)
  }
  if(is.null(nu)){
    if(diag){
      nu = rep(0,m)
    }else{
      nu = 0
    }
  }
  
  if(diag){
    est <- BReg_diag(X, Y, Int, Bp, lam, nu, reps, burn)
  }else{
    est <- BReg(X, Y, Int, Bp, lam, nu, reps, burn)
  }
  
  return(est)
  
}