library(simdd)
library(fastmatrix)
library(matrixcalc)
library(pracma)

gen <- ""
Y <- list()

for(i in 1:10) {
  
  dat <- read.table(file=paste(c("filename",i-1,".txt"),collapse=""))
  dat <- apply(as.matrix.noquote(dat),2,as.numeric)
  
  d <- length(dat[1,])
  n <- length(dat[,1])
  
  trans_bas <- matrix(0,nrow=d,ncol=d)
  
  for(j in 1:(d-1)){
    vectemp <- rep(0,d)
    vectemp[j] <- 1
    vectemp[d] <- -1
    trans_bas[,j] <- vectemp
  }
  trans_bas[,d] <- c(rep(0,d-1),1)
  
  for(j in 2:d){
    for(k in 1:(j-1)){
      trans_bas[,j] <- trans_bas[,j]-dot(trans_bas[,j],trans_bas[,k])/dot(trans_bas[,k],trans_bas[,k])*trans_bas[,k]
    }
  }
  for(j in 1:d){
    trans_bas[,j] <- trans_bas[,j]/norm(trans_bas[,j],type="2")
  }
  
  trans_bas <- solve(trans_bas)
  
  for(j in 1:n) {
    dat[j,] <- trans_bas%*%dat[j,]
  }
  dat <- dat[,-d]
  
  X <- dat
  d <- d-1
  
  matgen <- function(i,j,d) {
    res <- matrix(0,nrow=d,ncol=d)
    res[i,j] <- 1
    res[j,i] <- 1
    return(res)
  }
  
  ddf2 <- matrix(0,nrow=d*(d+1)/2,ncol=d^2)
  lapf2 <- rep(0,d*(d+1)/2)
  
  indi <- 1
  for(k in 1:d) {
    ddf2[indi,] <- 2*fastmatrix::vec(matgen(k,k,d))
    lapf2[indi] <- 2
    if((d-k)>=1) {
      for(l in 1:(d-k)) {
        ddf2[indi+l,] <- fastmatrix::vec(matgen(k+l,k,d))
      }
    }
    indi <- indi+d-k+1
  }
  
  
  df2_eval <- list()
  
  for(j in 1:n) {
    mat <- matrix(0,ncol=d,nrow=d*(d+1)/2)
    indi <- 1
    for(k in 1:d) {
      mat[indi,k] <- 2*X[j,k]
      if((d-k)>=1) {
        for(l in 1:(d-k)) {
          mat[indi+l,k+l] <- X[j,k]
          mat[indi+l,k] <- X[j,k+l]
        }
      }
      indi <- indi+d-k+1
    }
    df2_eval[[j]] <- mat
  }
  
  dd <- duplication.matrix(d)
  
  M <- matrix(0,nrow=d*(d+1)/2,ncol=d*(d+1)/2)
  for(j in 1:n) {
    M <- M+2*kronecker.prod(df2_eval[[j]]%*%(diag(d)-X[j,]%*%t(X[j,])),t(X[j,]))%*%dd
  }
  M <- M/n
  M <- M[-d*(d+1)/2,]
  M <- M[,-d*(d+1)/2]
  
  D <- rep(0,d*(d+1)/2)
  for(j in 1:n) {
    D <- D + (d-1)*df2_eval[[j]]%*%X[j,]+ddf2%*%kronecker.prod(X[j,],X[j,])-lapf2 
  }
  D <- D/n
  D <- D[-d*(d+1)/2]
  
  
  E <- matrix(0,nrow=d*(d+1)/2,ncol=d)
  for(j in 1:n) {
    E <- E + df2_eval[[j]]%*%(diag(d)-X[j,]%*%t(X[j,]))
  }
  E <- E/n
  E <- E[-d*(d+1)/2,]
  
  G <- matrix(0,nrow=d,ncol=d*(d+1)/2)
  for(j in 1:n) {
    G <- G + 2*kronecker.prod(diag(d)%*%(diag(d)-X[j,]%*%t(X[j,])),t(X[j,]))%*%dd
  }
  G <- G/n
  G <- G[,-d*(d+1)/2]
  
  H <- rep(0,d)
  for(j in 1:n) {
    H <- H + (d-1)*diag(d)%*%X[j,]
  }
  H <- H/n
  
  L <- matrix(0,nrow=d,ncol=d)
  for(j in 1:n) {
    L <- L+diag(d)%*%(diag(d)-X[j,]%*%t(X[j,]))
  }
  L <- L/n
  
  mu_est <- as.numeric(solve(L-G%*%solve(M)%*%E)%*%(H-G%*%solve(M)%*%D))
  
  
  A_est <- solve(M)%*%(D-E%*%mu_est)
  A_est[d*(d+1)/2] <- 0
  A_est <- dd%*%A_est
  A_est <- matrix(A_est,ncol=d)
  
  Y[[i]] <- rFisherBingham(10,mu_est,A_est,mtop=10000000)
  Y[[i]] <- cbind(Y[[i]],rep(0,10))
  
  for(j in 1:10) {
    Y[[i]][j,] <- solve(trans_bas)%*%Y[[i]][j,]
  }
  
  for(j in 1:length(Y[[i]][,1])) {
    gen <- paste(c(gen,Y[[i]][j,]),collapse=" ")
    gen <- paste(gen,"\n")
  }
}

sink("gen_sph_fisherbingham.txt")
cat(gen)
sink()
