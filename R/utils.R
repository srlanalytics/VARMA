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

