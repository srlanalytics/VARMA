library(dateutils)

A <- matrix(c(.5, .1, .2, .4), 2, 2)
B <- matrix(c(1, 0, .4, .2, 0, 1, .1, .3),2,4, byrow = TRUE)

r <- 1000

E <- stack_obs(matrix(rnorm((r+1)*2), r+1, 2), 2)

X <- matrix(0, r, 2)

for(j in 2:r){
  X[j, ] = A%*%X[j-1, ] + B%*%E[j, ]
}

out <- VARMA_select(X)
A
out$P

B[,3:4]
out$Q

data <- X
P <- out$P
Q <- out$Q

library(devtools)
load_all("~/VARMA")