
VARMA_loss <- function(x, data, order, shrink_P = 0, shrink_Q = 0){
  
  
  n_obs <- sum(is.finite(data))
  m  <- NCOL(data)
  idx <- m*m*order[1]
  if(order[1]>0){
    P  <- matrix(x[1:idx], m, m*order[1])
    eg <- eigen(comp_form(P))
    if(any(abs(eg$values)>1)){
      return(1e10)
    }
  }else{
    P <- matrix(0,m,0)
  }
  if(order[2]>0){
    Q  <- matrix(x[(idx+1):(idx+m*m*order[2])], m, m*order[2])
  }else{
    Q <- matrix(0,m,0)
  }
  
  MSE <- VARMA_MSE(P,Q,Y = data)
  
  loss <- shrink_P*sum(P^2) + shrink_Q*sum(Q^2) + MSE$MSE/n_obs
  
  return(loss)
}

# VARMA_loss <- function(x, data, order, shrink_P = 0, shrink_Q = 0){
#   
#   m  <- NCOL(data)
#   idx <- m*m*order[1]
#   if(order[1]>0){
#     P  <- matrix(x[1:idx], m, m*order[1])
#   }
#   if(order[2]>0){
#     Q  <- matrix(x[(idx+1):(idx+m*m*order[2])], m, m*order[2])
#   }
#   idx <- idx+m*m*order[2]
#   E <- matrix(0,m,m)
#   E[upper.tri(E, diag = T)] <- x[(idx+1):length(x)]
#   E <- E + t(E) - diag(diag(E),m,m)
# 
#   
#   #Put into state space form
#   r  <- NROW(data)
#   sP <- order[1]*m
#   sQ <- order[2]*m
#   k <- max(sP, sQ+m) 
#   
#   B <- matrix(0,m,k)
#   H <- cbind(diag(1,m,m),matrix(0,m,k-m))
#   
#   if(sP>0){
#     B[,1:sP] <- P
#   }
#   if(sQ>0){
#     H[,(m+1):(sQ+m)] <- Q
#   }
#   
#   lik <- KLike(B,E,H,R = rep(0,m),Y = data)
#   
#   loss <- r*shrink_P*sum(P^2) + r*shrink_Q*sum(Q^2) - lik
#   
#   return(loss)
# }



VARMA_init <- function(data, order, shrink_P = 0, shrink_Q = 0, method = "BFGS"){
  
  init <- apply(data, MARGIN = 2, FUN = VARMA_est, order = order, shrink_P = shrink_P, shrink_Q = shrink_Q)
  
  P <- expand_mat(do.call("rbind", lapply(init, FUN = get_from_list, what = "P")))
  Q <- expand_mat(do.call("rbind", lapply(init, FUN = get_from_list, what = "Q")))
  
  out <- VARMA_est(data, order, shrink_P = shrink_P, shrink_Q = shrink_Q, P_in = P, Q_in = Q, method = method)
  return(out)
}

VARMA_est <- function(data, order, shrink_P = 0, shrink_Q = 0, P_in = NULL, Q_in = NULL, method = "BFGS"){
  
  data <- as.matrix(data)
  m <- NCOL(data)
  
  if(is.null(P_in)){
    P <- matrix(0, m, m*order[1])
  }else{
    P <- P_in
  }
  if(is.null(Q_in)){
    Q <- matrix(0,m,m*order[2])
  }else{
    Q <- Q_in
  }
  
  x <- c(P,Q)
  
  if(length(x)==1){
    parameters <- optim(par = x, fn = VARMA_loss, data = data, order = order, 
                        shrink_P = shrink_P, shrink_Q = shrink_Q, method = "Brent", lower = -1, upper = 1)
  }else{
    parameters <- optim(par = x, fn = VARMA_loss, data = data, order = order, 
                        shrink_P = shrink_P, shrink_Q = shrink_Q, method = method)
  }
  
  
  
  x <- parameters$par
  
  idx <- m*m*order[1]
  if(order[1]>0){
    P  <- matrix(x[1:idx], m, m*order[1])
  }
  if(order[2]>0){
    Q  <- matrix(x[(idx+1):(idx+m*m*order[2])], m, m*order[2])
  }
  
  out <- list(P = P,
              Q = Q,
              MSE = parameters$value)
  
  return(out)
}

# VARMA_est <- function(data, order, shrink_P = 0, shrink_Q = 0){
#   
#   m <- NCOL(data)
#   P <- matrix(0, m, m*order[1])
#   Q <- matrix(0,m,m*order[2])
#   E <- var(data, use = "pairwise.complete.obs")/2
#   
#   x <- c(P,Q,E[upper.tri(E, diag = TRUE)])
#   
#   parameters <- optim(par = x, fn = VARMA_loss, data = data, order = order, 
#                       shrink_P = shrink_P, shrink_Q = shrink_Q)
#   
#   x <- parameters$par
#   idx <- m*m*order[1]
#   if(order[1]>0){
#     P  <- matrix(x[1:idx], m, m*order[1])
#   }
#   if(order[2]>0){
#     Q  <- matrix(x[(idx+1):(idx+m*m*order[2])], m, m*order[2])
#   }
#   idx <- idx+m*m*order[2]
#   E <- matrix(0,m,m)
#   E[upper.tri(E, diag = T)] <- x[(idx+1):length(x)]
#   E <- E + t(E) - diag(diag(E),m,m)
#   
#   out <- list(P = P,
#               Q = Q,
#               E = E,
#               Lik = -parameters$value)
#   
#   return(out)
# }


VARMA_select <- function(data, shrink_P = 0, shrink_Q = 0, pen = 0.1){
 
  m <- NCOL(data)
  n_obs <- sum(is.finite(data))
 
  estAR <- VARMA_est(data, c(1,0), shrink_P = shrink_P, shrink_Q = shrink_Q)
  ucv <- sum(data^2, na.rm = TRUE)
  SICar <- (ucv/n_obs-estAR$MSE)/(ucv/n_obs) - (pen*m^2)/sqrt(n_obs)  #Seth's bogus information criteria
  #BICar <- log(n_obs) * (2*m^2) - 2 * estAR$lik
  
  estMA <- VARMA_est(data, c(0,1), shrink_P = shrink_P, shrink_Q = shrink_Q)
  SICma <- (ucv/n_obs-estMA$MSE)/(ucv/n_obs) - (pen*m^2)/sqrt(n_obs)  #Seth's bogus information criteria
  
  if(SICar>SICma){
    order <- c(1,0)
    SIC <- SICar
    est <- estAR
  }else{
    order <- c(0,1)
    SIC <- SICma
    est <- estMA
  }
  
  for(it in 1:14){
    
    orderAR <- order + c(1,0)
    orderMA <- order + c(0,1)
    
    estAR <- VARMA_est(data, orderAR, shrink_P = shrink_P, shrink_Q = shrink_Q)
    SICar <- (ucv/n_obs-estAR$MSE)/(ucv/n_obs) - (pen*sum(orderAR)*m^2)/sqrt(n_obs)
    
    estMA <- VARMA_est(data, orderMA, shrink_P = shrink_P, shrink_Q = shrink_Q)
    SICma <- (ucv/n_obs-estMA$MSE)/(ucv/n_obs) - (pen*sum(orderMA)*m^2)/sqrt(n_obs)  #Seth's bogus information criteria
    
    if(SIC > max(SICar,SICma)){
      break
    }
    
    if(SICar>SICma){
      order <- orderAR
      SIC <- SICar
      est <- estAR
    }else{
      order <- orderMA
      SIC <- SICma
      est <- estMA
    }
  }
    
    Out <- list(order = order,
                P = est$P,
                Q = est$Q,
                MSE = est$MSE)
    
}


VARMA_select_multi <- function(data, models = 3, shrink_P = 0, shrink_Q = 0, pen = 0.1){
  
  m <- NCOL(data)
  n_obs <- sum(is.finite(data))
  scores <- matrix(0,8,8)
  
  estAR <- VARMA_est(data, c(1,0), shrink_P = shrink_P, shrink_Q = shrink_Q)
  ucv <- sum(data^2, na.rm = TRUE)
  SICar <- (ucv/n_obs-estAR$MSE)/(ucv/n_obs) - (pen*m^2)/sqrt(n_obs)  #Seth's bogus information criteria
  #BICar <- log(n_obs) * (2*m^2) - 2 * estAR$lik
  
  estMA <- VARMA_est(data, c(0,1), shrink_P = shrink_P, shrink_Q = shrink_Q)
  SICma <- (ucv/n_obs-estMA$MSE)/(ucv/n_obs) - (pen*m^2)/sqrt(n_obs)  #Seth's bogus information criteria
  
  scores[2,1] <- SICar
  scores[1,2] <- SICma
  
  if(SICar>SICma){
    order <- c(1,0)
    SIC <- SICar
    est <- estAR
  }else{
    order <- c(0,1)
    SIC <- SICma
    est <- estMA
  }
  
  count <- 0
  it <- 0
  while(count < models && it < 8){
    orderAR <- order + c(1,0)
    orderMA <- order + c(0,1)
    
    estAR <- VARMA_est(data, orderAR, shrink_P = shrink_P, shrink_Q = shrink_Q)
    SICar <- (ucv/n_obs-estAR$MSE)/(ucv/n_obs) - (pen*sum(orderAR)*m^2)/sqrt(n_obs)
    
    estMA <- VARMA_est(data, orderMA, shrink_P = shrink_P, shrink_Q = shrink_Q)
    SICma <- (ucv/n_obs-estMA$MSE)/(ucv/n_obs) - (pen*sum(orderMA)*m^2)/sqrt(n_obs)  #Seth's bogus information criteria
    
    scores[orderAR[1]+1, orderAR[2]+1] <- SICar
    scores[orderMA[1]+1, orderMA[2]+1] <- SICma
    
    if(max(scores) > max(SICar,SICma)){
      count <- count+1
    }
    
    if(SICar>SICma){
      order <- orderAR
      SIC <- SICar
      est <- estAR
    }else{
      order <- orderMA
      SIC <- SICma
      est <- estMA
    }
    it = it+1
  }
  
  
  return(scores)
  
  # Out <- list(order = order,
  #             P = est$P,
  #             Q = est$Q,
  #             MSE = est$MSE)
  
}

VARMA_SS <- function(data, P, Q, Sig = NULL){
  
  if(is.null(Sig)){
    est_basic <- VARMA_MSE(P,Q,data) 
    E <- var((data - est_basic$YP), use = "complete.obs")
  }else{
    E <- Sig
  }
  E <- as.matrix(E)
  
  m <- NCOL(data)
  sP <- NCOL(P)
  sQ <- NCOL(Q)
  k <- max(sP, sQ+m)
  B <- matrix(0,m,k)
  H <- cbind(diag(1,m,m),matrix(0,m,k-m))
  if(sP>0){
    B[,1:sP] <- P
  }
  if(sQ>0){
    H[,(m+1):(sQ+m)] <- Q
  }
  
  SS <- Kfilter(Y = t(data), B,E,H,R = rep(0,m))
  
  Out <- list(Lik = SS$Lik,
              Yfit = t(SS$Yfit),
              var = SS$var,
              Z  = t(SS$Z),
              B = B,
              E = E,
              H = H)
  return(Out)
}



