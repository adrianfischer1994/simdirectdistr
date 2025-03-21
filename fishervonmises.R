starttime <- Sys.time()

library(parallel)
numCores <- 10
cl <- makeCluster(numCores)

clusterEvalQ(cl,{
  
  set.seed(1)
  
  library(vMF)
  library(nleqslv)
  
  m <- 10000
  n <- 100
  
  #GSM functions
  fr <- function(x) {
    return((-1-0.01*x)*log(1+0.01*x))
  }
  dfr <- function(x) {
    return(-0.01*(1+log(1+0.01*x)))
  }
  ddfr <- function(x) {
    return(-0.01^2/(1+0.01*x))
  }
  
})

sim <- function(param) {
  
  mu <- param[[1]]
  mu <- mu/norm(mu, type="2")
  kappa <- param[[2]]
  d <- length(mu)
  
  MLE <- rep(NA,m)
  SM <- rep(NA,m)
  ST <- rep(NA,m)
  GSM <- rep(NA,m)
  # ST2 <- rep(NA,m)
  
  assign("MLE_nonex", 0, env=globalenv())
  assign("SM_nonex", 0, env=globalenv())
  assign("ST_nonex", 0, env=globalenv())
  assign("GSM_nonex", 0, env=globalenv())
  # assign("ST2_nonex", 0, env=globalenv())
  
  for(i in 1:m) {
    
    mu_t <- kappa*mu/norm(mu, type="2")
    X <- rvMF(n,mu_t)
    
    r <- norm(colSums(X)/n, type="2")
    mu_est <- colSums(X)/n*1/r 
    
    #MLE
    func_opt <- function(x) {
      res <- besselI(x,d/2)/besselI(x,d/2-1)-r
      return(res)
    }
    MLE_kappa <- nleqslv(1,func_opt)$x
    if(MLE_kappa>0) {
      MLE[i] <- MLE_kappa
    } else {
      nonex <- get("MLE_nonex", env=globalenv())
      assign("MLE_nonex", nonex+1, env=globalenv())
    }
    
    #SM
    tryCatch({
      sp <- c(1,rep(0,d-1))
      mu_est_o <- mu_est-as.numeric(t(sp)%*%mu_est)*sp
      mu_est_o <- mu_est_o/norm(mu_est_o, type="2")
      theta <- as.numeric(acos(t(sp)%*%mu_est))
      rot_mat <- diag(d)+sin(theta)*(sp%*%t(mu_est_o)-mu_est_o%*%t(sp))+(cos(theta)-1)*(mu_est_o%*%t(mu_est_o)+sp%*%t(sp))
      
      X_rot <- X%*%t(rot_mat)
      
      SM_kappa <- (d-1)*mean(X_rot[,1])/(1-mean(X_rot[,1]^2))
      
      if(SM_kappa>0) {
        SM[i] <- SM_kappa
      } else {
        nonex <- get("SM_nonex", env=globalenv())
        assign("SM_nonex", nonex+1, env=globalenv())
      }
    },error=function(cond) {
      nonex <- get("SM_nonex", env=globalenv())
      assign("SM_nonex", nonex+1, env=globalenv())
    })
    
    #ST
    tryCatch({
      mat <- matrix(0,nrow = d,ncol = d)
      for(j in 1:n) {
        mat <- mat+X[j,]%*%t(X[j,])
      }
      mat <- diag(d)-mat/n
      
      M <- as.vector((t(mu_est)%*%mat))
      D <- as.vector((d-1)*colMeans(X))
      
      ST_kappa <- as.vector(t(M)%*%D)/as.vector(t(M)%*%M)
      
      if(ST_kappa>0) {
        ST[i] <- ST_kappa
      } else {
        nonex <- get("ST_nonex", env=globalenv())
        assign("ST_nonex", nonex+1, env=globalenv())
      }
    },error=function(cond) {
      nonex <- get("ST_nonex", env=globalenv())
      assign("ST_nonex", nonex+1, env=globalenv())
    })
    
    #GSM (generalized score matching)
    tryCatch({

      func_opt <- function(x) {
        temp <- rep(NA,n)
        for(i in 1:n) {
          temp[i] <- t(mu_est)%*%X[i,]
        }
        return(mean(fr(x^2*(1-temp^2))-2*dfr(x^2*(1-temp^2))*(x^2*(1-temp^2)-(d-1)*x*temp)+4*ddfr(x^2*(1-temp^2))*x^3*temp*(1-temp^2)))
      }

      GSM_kappa <- optim(c(1), func_opt, method = "Brent",lower=0,upper=200)$par

      if(GSM_kappa>0) {
        GSM[i] <- GSM_kappa
      } else {
        nonex <- get("GSM_nonex", env=globalenv())
        assign("GSM_nonex", nonex+1, env=globalenv())
      }
    },error=function(cond) {
      nonex <- get("GSM_nonex", env=globalenv())
      assign("GSM_nonex", nonex+1, env=globalenv())
    })

    # #ST2
    # tryCatch({
    #   mat <- matrix(0,nrow = d,ncol = d)
    #   for(j in 1:n) {
    #     mat <- mat+X[j,]%*%t(X[j,])
    #   }
    #   mat <- diag(d)-mat/n
    #   
    #   ST2_kappa <- (d-1)*norm(solve(mat)%*%colSums(X)/n, type="2")
    #   
    #   if(ST2_kappa>0) {
    #     ST2[i] <- ST2_kappa
    #   } else {
    #     nonex <- get("ST2_nonex", env=globalenv())
    #     assign("ST2_nonex", nonex+1, env=globalenv())
    #   }
    # },error=function(cond) {
    #   nonex <- get("ST2_nonex", env=globalenv())
    #   assign("ST2_nonex", nonex+1, env=globalenv())
    # })
  }
  
  bias <- matrix(0,4,1)
  mse <- matrix(0,4,1)
  ex <- numeric(4)
  sd_bias <- matrix(0,4,1)
  sd_mse <- matrix(0,4,1)
  
  bias[1,] <- mean(MLE-kappa,na.rm=TRUE)
  bias[2,] <- mean(SM-kappa,na.rm=TRUE)
  bias[3,] <- mean(ST-kappa,na.rm=TRUE)
  bias[4,] <- mean(GSM-kappa,na.rm=TRUE)
  # bias[5,] <- mean(ST2-kappa,na.rm=TRUE)

  mse[1,] <- mean((MLE-kappa)^2,na.rm=TRUE)
  mse[2,] <- mean((SM-kappa)^2,na.rm=TRUE)
  mse[3,] <- mean((ST-kappa)^2,na.rm=TRUE)
  mse[4,] <- mean((GSM-kappa)^2,na.rm=TRUE)
  # mse[5,] <- mean((ST2-kappa)^2,na.rm=TRUE)
  
  sd_bias[1,] <- sqrt(var(MLE-kappa,na.rm=TRUE)/m)
  sd_bias[2,] <- sqrt(var(SM-kappa,na.rm=TRUE)/m)
  sd_bias[3,] <- sqrt(var(ST-kappa,na.rm=TRUE)/m)
  sd_bias[4,] <- sqrt(var(GSM-kappa,na.rm=TRUE)/m)
  # sd_bias[5,] <- sqrt(var(ST2-kappa,na.rm=TRUE)/m)
  
  sd_mse[1,] <- sqrt(var((MLE-kappa)^2,na.rm=TRUE)/m)
  sd_mse[2,] <- sqrt(var((SM-kappa)^2,na.rm=TRUE)/m)
  sd_mse[3,] <- sqrt(var((ST-kappa)^2,na.rm=TRUE)/m)
  sd_mse[4,] <- sqrt(var((GSM-kappa)^2,na.rm=TRUE)/m)
  # sd_mse[5,] <- sqrt(var((ST2-kappa)^2,na.rm=TRUE)/m)
  
  ex[1] <- MLE_nonex
  ex[2] <- SM_nonex
  ex[3] <- ST_nonex
  ex[4] <- GSM_nonex
  # ex[5] <- ST2_nonex
  
  return(list(bias,mse,ex,sd_bias,sd_mse))
}

#mu does not need to be normalized

param <- list(list(rep(1,3),1),list(rep(1,3),2),list(rep(1,3),10),list(rep(1,3),50),
              list(rep(1,10),1),list(rep(1,10),10),list(rep(1,10),50),
              list(rep(1,20),1),list(rep(1,20),10),list(rep(1,20),50))
erg <- parLapply(cl,param,sim)

stopCluster(cl)

  
saveRDS(erg, file="result.RData")

nb_est <- length(erg[[1]][[1]][,1]) #number of estimators to compare
nb_par <- length(erg[[1]][[1]][1,]) #number of unknown parameters

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
for(i in 1:(nb_est*2+2)) {
  table <- paste(table,"c",sep="")
}
table <- paste(table,"}\n $\\theta_0$ & & \\multicolumn{",nb_est,"}{c}{Bias} & \\multicolumn{",nb_est,"}{c}{MSE} \\\\ \n",sep="")
for(i in 1:length(erg)) {
  table <- paste(table,"\\multirow{1}{*}{$(",param[[i]][1],",",param[[i]][2],")$} & $\\kappa$ ",sep="")
  for(j in 1:length(trans1[1,])) { 
    ent <- rnd(trans1[i,j])
    if(j == (2*nb_est*nb_par) ) {
      table <- paste(table," & $",ent,"$ \\\\", sep="")
    } else {
      table <- paste(table," & $",ent,"$", sep="")
    }
  }
  table <- paste(table," \\hline \n", sep="")
}

table <- paste(table,"\\end{tabular} \n\\end{table}", sep="")

sink(paste("table1.txt", sep=""))
cat(table)
sink()



table <- "\\begin{table} \n\\centering\n\\begin{tabular}{"
for(i in 1:(nb_est*2+2)) {
  table <- paste(table,"c",sep="")
}
table <- paste(table,"}\n $\\theta_0$ & & \\multicolumn{",nb_est,"}{c}{Bias} & \\multicolumn{",nb_est,"}{c}{MSE} \\\\ \n",sep="")
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

sink(paste("table2.txt", sep=""))
cat(table)
sink()
  

endtime <- Sys.time()
time <- endtime-starttime
