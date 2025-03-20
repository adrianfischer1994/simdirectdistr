starttime <- Sys.time()

library(parallel)
numCores <- 10
cl <- makeCluster(numCores)

clusterEvalQ(cl,{
  
  set.seed(1)
  
  library(watson)
  library(hypergeo)
  library(fastmatrix)
  
  n <- 100
  
  ll <- function(X,kappa,mu) {
    d <- length(mu)
    c <- n*log(gamma(d/2)/(2*pi^(d/2)*genhypergeo(U=c(0.5), L=c(d/2),z=kappa)))
    tmp <- 0
    for(i in 1:n) {
      tmp <- tmp+(t(mu)%*%X[i,])^2
    }
    return(c+kappa*tmp)
  }
  
  matgen <- function(i,j,d) {
    res <- matrix(0,nrow=d,ncol=d)
    res[i,j] <- 1
    res[j,i] <- 1
    return(res)
  }
  
})

sim <- function(param) {

  mu <- param[[1]]
  mu <- mu/norm(mu, type="2")
  kappa <- param[[2]]
  d <- length(mu)
  
  MLEa <- rep(NA,m)
  EM <- rep(NA,m)
  ST <- rep(NA,m)
  
  assign("MLEa_nonex", 0, env=globalenv())
  assign("EM_nonex", 0, env=globalenv())
  assign("ST_nonex", 0, env=globalenv())
  
  for(i in 1:m) {
    
    X <- rmwat(n,weights=1,kappa=kappa,mu=matrix(mu,nrow=d))
    
    #eigenvalues
    scat <- matrix(0,d,d)
    for(j in 1:n) {
      scat <- scat+X[j,]%*%t(X[j,])
    }
    scat <- scat/n
    eg <- eigen(scat)
    mu_est1 <- eg$vectors[,1]
    mu_est2 <- eg$vectors[,d]
    
    #MLE_approx
    tryCatch({
      r1 <- eigen(scat)$values[1]
      r2 <- eigen(scat)$values[d]
      L1 <- (r1*d/2-0.5)/(r1*(1-r1))*(1+(1-r1)/(d/2-0.5))
      L2 <- (r2*d/2-0.5)/(r2*(1-r2))*(1+(1-r2)/(d/2-0.5))
      U1 <- (r1*d/2-0.5)/(r1*(1-r1))*(1+2*r1)
      U2 <- (r2*d/2-0.5)/(r2*(1-r2))*(1+2*r2)
      MLEa_kappa1 <- 0.5*(L1+U1)
      MLEa_kappa2 <- 0.5*(L2+U2)
      if(MLEa_kappa1>0) {
        if(MLEa_kappa2<0) {
          ll1 <- ll(X,MLEa_kappa1,mu_est1)
          ll2 <- ll(X,MLEa_kappa2,mu_est2)
          if(ll1>ll2) {
            MLEa[i] <- MLEa_kappa1
          } else {
            MLEa[i] <- MLEa_kappa2
          }
        } else {
          MLEa[i] <- MLEa_kappa1
        }
      } else {
        if(MLEa_kappa2<0) {
          MLEa[i] <- MLEa_kappa2
        } else {
          nonex <- get("MLEa_nonex", env=globalenv())
          assign("MLEa_nonex", nonex+1, env=globalenv())
        }
      }
    },error=function(cond) {
      nonex <- get("MLEa_nonex", env=globalenv())
      assign("MLEa_nonex", nonex+1, env=globalenv())
    })
    
    
    #EM
    tryCatch({
      EM_kappa <- watson(X,1)$kappa
      EM[i] <- EM_kappa
    },error=function(cond) {
      nonex <- get("EM_nonex", env=globalenv())
      assign("EM_nonex", nonex+1, env=globalenv())
    })
    
    #ST
    tryCatch({
      
      ddf <- matrix(0,nrow=d*(d+1)/2,ncol=d^2)
      lapf <- rep(0,d*(d+1)/2)
      
      indi <- 1
      for(k in 1:d) {
        ddf[indi,] <- 2*fastmatrix::vec(matgen(k,k,d))
        lapf[indi] <- 2
        if((d-k)>=1) {
          for(l in 1:(d-k)) {
            ddf[indi+l,] <- fastmatrix::vec(matgen(k+l,k,d))
          }
        }
        indi <- indi+d-k+1
      }
      
      
      df_eval <- list()
      
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
        df_eval[[j]] <- mat
      }
      
      
      J1 <- rep(0,d*(d+1)/2)
      for(j in 1:n) {
        J1 <- J1+2*t(X[j,])%*%mu_est1%*%t(mu_est1)%*%(diag(d)-X[j,]%*%t(X[j,]))%*%t(df_eval[[j]])
      }
      J1 <- as.vector(J1/n)
      J1 <- J1[-d*(d+1)/2]
      
      J2 <- rep(0,d*(d+1)/2)
      for(j in 1:n) {
        J2 <- J2+2*t(X[j,])%*%mu_est2%*%t(mu_est2)%*%(diag(d)-X[j,]%*%t(X[j,]))%*%t(df_eval[[j]])
      }
      J2 <- as.vector(J2/n)
      J2 <- J2[-d*(d+1)/2]
      
      V <- rep(0,d*(d+1)/2)
      for(j in 1:n) {
        V <- V + (d-1)*df_eval[[j]]%*%X[j,]+ddf%*%kronecker.prod(X[j,],X[j,])-lapf 
      }
      V <- V/n
      V <- V[-d*(d+1)/2]
      
      ST_kappa1 <- as.numeric(solve(t(J1)%*%J1)%*%t(J1)%*%V)
      ST_kappa2 <- as.numeric(solve(t(J2)%*%J2)%*%t(J2)%*%V)
      
      if(ST_kappa1>0) {
        if(ST_kappa2<0) {
          st1 <- norm(J1*ST_kappa1-V,type="2")
          st2 <- norm(J2*ST_kappa2-V,type="2")
          if(st1>st2) {
            ST[i] <- ST_kappa2
          } else {
            ST[i] <- ST_kappa1
          }
        } else {
          ST[i] <- ST_kappa1
        }
      } else {
        if(ST_kappa2<0) {
          ST[i] <- ST_kappa2
        } else {
          nonex <- get("ST_nonex", env=globalenv())
          assign("ST_nonex", nonex+1, env=globalenv())
        }
      }
      
    },error=function(cond) {
      nonex <- get("ST_nonex", env=globalenv())
      assign("ST_nonex", nonex+1, env=globalenv())
    })
  }
  
  bias <- matrix(0,3,1)
  mse <- matrix(0,3,1)
  ex <- numeric(3)
  
  sd_bias <- matrix(0,3,1)
  sd_mse <- matrix(0,3,1)
  
  bias[1,] <- mean(MLEa-kappa,na.rm=TRUE)
  bias[2,] <- mean(EM-kappa,na.rm=TRUE)
  bias[3,] <- mean(ST-kappa,na.rm=TRUE)
  
  mse[1,] <- mean((MLEa-kappa)^2,na.rm=TRUE)
  mse[2,] <- mean((EM-kappa)^2,na.rm=TRUE)
  mse[3,] <- mean((ST-kappa)^2,na.rm=TRUE)
  
  sd_bias[1,] <- sqrt(var(MLEa-kappa,na.rm=TRUE)/m)
  sd_bias[2,] <- sqrt(var(EM-kappa,na.rm=TRUE)/m)
  sd_bias[3,] <- sqrt(var(ST-kappa,na.rm=TRUE)/m)
  
  sd_mse[1,] <- sqrt(var((MLEa-kappa)^2,na.rm=TRUE)/m)
  sd_mse[2,] <- sqrt(var((EM-kappa)^2,na.rm=TRUE)/m)
  sd_mse[3,] <- sqrt(var((ST-kappa)^2,na.rm=TRUE)/m)
  
  ex[1] <- MLEa_nonex
  ex[2] <- EM_nonex
  ex[3] <- ST_nonex
  
  return(list(bias,mse,ex,sd_bias,sd_mse))
}

m <- 10000
clusterExport(cl,c("m"))

#mu does not need to be normalized

param <- list(list(rep(1,3),-20),list(rep(1,3),-10),list(rep(1,10),-10),list(rep(1,10),-2),
              list(rep(1,20),-2),list(rep(1,3),1),list(rep(1,10),1),
              list(rep(1,20),5),list(rep(1,3),10),list(rep(1,10),20))
erg <- parLapply(cl,param,sim)

stopCluster(cl)

saveRDS(erg, file="result.RData")

nb_est <- length(erg[[1]][[1]][,1]) #number of estimators to compare
nb_par <- length(erg[[1]][[1]][1,]) #number of unknown parameters
nb_parco <- length(erg) #number of parameter constellations

trans1 <- matrix(0,length(erg),nb_est*nb_par*2)
for(i in 1:length(erg)) {
  for(j in 1:nb_par) {
    for(l in 1:2) {
      for(k in 1:nb_est) {
        trans1[i,(j-1)*2*nb_est+(l-1)*nb_est+k] <- erg[[i]][[l]][k,j]
      }    
    }
  }
}

trans2 <- matrix(0,length(erg),nb_est*nb_par*2)
for(i in 1:length(erg)) {
  for(j in 1:nb_par) {
    for(l in 4:5) {
      for(k in 1:nb_est) {
        trans2[i,(j-1)*2*nb_est+(l-4)*nb_est+k] <- erg[[i]][[l]][k,j]
      }    
    }
  }
}

trans_nonex <- matrix(0,nb_parco,nb_est)
for(i in 1:nb_parco) {
  for(j in 1:nb_est) {
    trans_nonex[i,j] <- erg[[i]][[3]][j]
  }
}

rnd <- function(ent) {
  if(is.nan(ent)) {
    return("NaN")
  }
  else { 
    if(abs(ent) < 0.01 | abs(ent) > 10000) {
      ex <- ifelse(ent == 0, 0, floor(log10(abs(ent))))
      mant <- round(ent/10^ex,2)
      ent <- paste(mant,"\\text{e",ex,"}",sep="")
    } else {
      if(1<abs(ent) & abs(ent)<10) {
        ent <- round(ent,2)
      } else {
        if(10<abs(ent) & abs(ent)<100) {
          ent <- round(ent,1)
        } else {
          if(100<abs(ent) & abs(ent)<10000) {
            ent <- round(ent,0)
          } else {
            ent <- round(ent,3)
          }
        }
      }
    }
    return(ent)
  }
}

table <- "\\begin{table} \n\\centering\n\\begin{tabular}{"
for(i in 1:(nb_est*3+2)) {
  table <- paste(table,"c",sep="")
}
table <- paste(table,"}\n $\\theta_0$ & & \\multicolumn{",nb_est,"}{c}{Bias} & \\multicolumn{",nb_est,"}{c}{MSE} & \\multicolumn{",nb_est,"}{c}{NE} \\\\ \n",sep="")
for(i in 1:length(erg)) {
  table <- paste(table,"\\multirow{1}{*}{$(",param[[i]][1],",",param[[i]][2],")$} & $\\kappa$ ",sep="")
  for(j in 1:length(trans1[1,])) { 
    ent <- rnd(trans1[i,j])
    if(j == 2*nb_est) {
      table <- paste(table," & $",ent,"$", sep="")
      for(k in 1:nb_est) {
        table <- paste(table," & \\multirow{",nb_par,"}{*}{$",round(trans_nonex[i,k]/m,2)*100,"$}", sep="")
      }
      table <- paste(table," \\\\", sep="")
    } else {
      table <- paste(table," & $",ent,"$", sep="")
    }
  }
  table <- paste(table," \\hline \n", sep="")
}

table <- paste(table,"\\end{tabular} \n\\end{table}", sep="")

sink("table1.txt")
cat(table)
sink()


table <- "\\begin{table} \n\\centering\n\\begin{tabular}{"
for(i in 1:(nb_est*2+2)) {
  table <- paste(table,"c",sep="")
}
table <- paste(table,"}\n $\\theta_0$ & & \\multicolumn{",nb_est,"}{c}{Bias} & \\multicolumn{",nb_est,"}{c}{MSE}  \\\\ \n",sep="")
for(i in 1:length(erg)) {
  table <- paste(table,"\\multirow{1}{*}{$(",param[[i]][1],",",param[[i]][2],")$} & $\\kappa$ ",sep="")
  for(j in 1:length(trans2[1,])) { 
    ent <- rnd(trans2[i,j])
    if(j == (2*nb_est*nb_par) ) {
      table <- paste(table," & $",ent,"$ \\\\", sep="")
    } else {
      table <- paste(table," & $",ent,"$", sep="")
    }
  }
  table <- paste(table," \\hline \n", sep="")
}

table <- paste(table,"\\end{tabular} \n\\end{table}", sep="")

sink("table2.txt")
cat(table)
sink()

endtime <- Sys.time()
time <- endtime-starttime
