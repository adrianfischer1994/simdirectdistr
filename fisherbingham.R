starttime <- Sys.time()


library(parallel)
numCores <- 10
cl <- makeCluster(numCores)

clusterEvalQ(cl,{
  
  set.seed(1)

  library(simdd)
  library(fastmatrix)
  library(matrixcalc)
  
  n <- 1000

  matgen <- function(i,j,d) {
    res <- matrix(0,nrow=d,ncol=d)
    res[i,j] <- 1
    res[j,i] <- 1
    return(res)
  }

})

sim <- function(param) {
  
  mu <- param[[1]]
  A <- param[[2]]
  
  ST <- list()
  
  assign("ST_nonex", 0, env=globalenv())
  
  d <- length(mu) 
  
  for(i in 1:m) {
    
    X <- rFisherBingham(n,mu,A)
      
    #ST
    
    tryCatch({
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
        ST[[length(ST)+1]] <- list(mu_est,A_est)
      } else {
        nonex <- get("ST_nonex", env=globalenv())
        assign("ST_nonex", nonex+1, env=globalenv())
      }
      
      
    },error=function(cond) {
      nonex <- get("ST_nonex", env=globalenv())
      assign("ST_nonex", nonex+1, env=globalenv())
    })
  
  }
  
  
  l <- length(ST)
  
  STmu_n <- numeric(l)
  STA_n <- numeric(l)
  
  for(i in 1:l) {
    STmu_n[i] <- norm(ST[[i]][[1]]-mu, type="2") 
    STA_n[i] <- spectral.norm(ST[[i]][[2]]-A) 
  }
  
  mse <- c(mean(STmu_n),mean(STA_n))
  
  sd_mse <- c(sqrt(var(STmu_n)/m),sqrt(var(STA_n)/m))
  
  ex <- ST_nonex
  
  return(list(mse,ex,sd_mse))
}

m <- 10000
clusterExport(cl,c("m"))

#A needs to be symmetric and last entry (d,d) must be equal to zero. mu is vector in R^d.
#Be aware of the sometimes different parametrisation: mu=kappa*mu2 and mu2 a unit vector!

param <- list(list(c(10,0,0),diag(c(0,0,0))),
              list(c(5,0,0),diag(c(0,0,0))),
              list(c(0.5,0,0),diag(c(0,0,0))),
              list(c(11,3,10),matrix(c(2,-2,1,-2,12,-2,1,-2,0),ncol=3)),
              list(c(0.05,0.05,0.05),matrix(c(1,2,3,2,6,7,3,7,0),ncol=3)),
              list(c(0,3,3),matrix(c(0,0,0,0,0,-3,0,-3,0),ncol=3)),
              list(mu <- c(0,1,1),matrix(c(-1,-2,-3,-2,5,-3,-3,-3,0),ncol=3)),
              list(c(0,-1,1),matrix(c(-1,-2,-3,-2,1,0,-3,0,0),ncol=3)),
              list(c(0,-1,1),matrix(c(-5,0,-1,0,1,0,-1,0,0),ncol=3)),
              list(mu <- c(0.1,0,0.1), matrix(c(-6,0,0,0,1,0,0,0,0),ncol=3)))
erg <- parLapply(cl,param,sim)

stopCluster(cl)

saveRDS(erg, file="result.RData")

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
  }
  return(ent)
}

nb_par <- length(erg)

table <- "\\begin{table} \n\\centering\n\\begin{tabular}{"
for(i in 1:(nb_par+2)) {
  table <- paste(table,"c",sep="")
}
table <- paste(table,"}\n ",sep="")
for(i in 1:nb_par) {
  table <- paste(table,"& Fig. ",i," ",sep="")
}
table <- paste(table,"\\\\ \n ",sep="")
for(i in 1:nb_par) {
  ent <- rnd(erg[[i]][[1]][1])
  table <- paste(table," & $",ent,"$",sep="")
}
table <- paste(table,"\\\\ \n ",sep="")

for(i in 1:nb_par) {
  ent <- rnd(erg[[i]][[1]][2])
  table <- paste(table," & $",ent,"$",sep="")
}
table <- paste(table,"\\\\ \n ",sep="")

for(i in 1:nb_par) {
  ent <- erg[[i]][[2]]
  table <- paste(table," & $",round(ent/m,2)*100,"$",sep="")
}
table <- paste(table,"\\\\ \\hline \n ",sep="")

table <- paste(table,"\\end{tabular} \n\\end{table}", sep="")

sink("table1.txt")
cat(table)
sink()



table <- "\\begin{table} \n\\centering\n\\begin{tabular}{"
for(i in 1:(nb_par+2)) {
  table <- paste(table,"c",sep="")
}
table <- paste(table,"}\n ",sep="")
for(i in 1:nb_par) {
  table <- paste(table,"& Fig. ",i," ",sep="")
}
table <- paste(table,"\\\\ \n ",sep="")
for(i in 1:nb_par) {
  ent <- rnd(erg[[i]][[3]][1])
  table <- paste(table," & $",ent,"$",sep="")
}
table <- paste(table,"\\\\ \n ",sep="")

for(i in 1:nb_par) {
  ent <- rnd(erg[[i]][[3]][2])
  table <- paste(table," & $",ent,"$",sep="")
}
table <- paste(table,"\\\\ \\hline \n ",sep="")

table <- paste(table,"\\end{tabular} \n\\end{table}", sep="")

sink("table2.txt")
cat(table)
sink()


endtime <- Sys.time()
time <- endtime-starttime

