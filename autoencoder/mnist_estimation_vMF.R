library(pracma)
library(vMF)
library(nleqslv)

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
  
  #MLE
  r <- norm(colSums(X)/n, type="2")
  mu_est <- colSums(X)/n/r
  func_opt <- function(x) {
    res <- besselI(x,d/2)/besselI(x,d/2-1)-r
    return(res)
  }
  kappa_est <- nleqslv(1,func_opt)$x
  
  mu_t <- kappa_est*mu_est
  Y[[i]] <- rvMF(10,mu_t)
  Y[[i]] <- cbind(Y[[i]],rep(0,10))
  
  for(j in 1:10) {
    Y[[i]][j,] <- solve(trans_bas)%*%Y[[i]][j,]
  }
  
  for(j in 1:length(Y[[i]][,1])) {
    gen <- paste(c(gen,Y[[i]][j,]),collapse=" ")
    gen <- paste(gen,"\n")
  }
  
}

sink("gen_vMF.txt")
cat(gen)
sink()
