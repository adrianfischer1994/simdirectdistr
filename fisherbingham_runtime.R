
library(parallel)
numCores <- 10
cl <- makeCluster(numCores)

clusterEvalQ(cl,{
  
  set.seed(1)

  library(simdd)
  library(fastmatrix)
  library(matrixcalc)

  matgen <- function(i,j,d) {
    res <- matrix(0,nrow=d,ncol=d)
    res[i,j] <- 1
    res[j,i] <- 1
    return(res)
  }

})

sim <- function(param) {
  
  d <- param[1]
  n <- param[2]
  mu <- c(1,rep(0,d-1))
  A <- diag(rep(0,d))
  
  time_vec <- rep(NA,m)
  
  assign("ST_nonex", 0, env=globalenv())
  
  for(i in 1:m) {
    
    X <- rFisherBingham(n,mu,A)
      
    #ST
    
    tryCatch({
      starttime <- Sys.time()
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
      
      if(all(!is.nan(mu_est)) & all(!is.na(mu_est)) & all(!is.nan(A_est)) & all(!is.na(A_est))) {
        endtime <- Sys.time()
        time_vec[i] <- as.double(endtime-starttime, units = "secs")
      } else {
        nonex <- get("ST_nonex", env=globalenv())
        assign("ST_nonex", nonex+1, env=globalenv())
      }
      
      
    },error=function(cond) {
      nonex <- get("ST_nonex", env=globalenv())
      assign("ST_nonex", nonex+1, env=globalenv())
    })
  
  }

  time_res <- mean(time_vec, na.rm=TRUE)
  ex <- ST_nonex
  
  return(c(time_res,ex))
}

m <- 200
clusterExport(cl,c("m"))

#A needs to be symmetric and last entry (d,d) must be equal to zero. mu is vector in R^d.
#Be aware of the sometimes different parametrisation: mu=kappa*mu2 and mu2 a unit vector!

d_grid <- c(3,9,15)
n_grid <- c(50,100,200,500,2000)

param_t <- unname(as.matrix(expand.grid(d_grid,n_grid)))
param <- list()
for(i in 1:length(param_t[,1])) {
  param[[i]] <- param_t[i,] 
}
erg <- parLapply(cl,param,sim)

stopCluster(cl)

saveRDS(erg, file="result.RData")

print <- ""
for(i in 1:length(d_grid)) {
  for(j in 1:length(n_grid)) {
    print <- paste(print,erg[[i+(j-1)*3]][1],collapse=" ")
  }
  print <- paste(print,"\n")
}

sink("result.txt")
cat(print)
sink()



