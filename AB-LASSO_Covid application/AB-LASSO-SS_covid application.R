#setwd("~/Downloads/Covid-19")
library(mvtnorm)
library(plm)
library(hdm)
library(glmnet)
library(matrixStats)
#library(MASS)
data = load("data_weekly_balanced.Rdata")
attach(sdf_week)
N = length(unique(fips))
T = length(unique(week))

Y = matrix(logdc, nrow = T, ncol = N)  #case level
D = matrix(school, nrow = T, ncol = N) 
ZZ = ZZ.t = list()
ZZ[[1]] = matrix(dlogtests, nrow = T, ncol = N)  
ZZ[[2]] = matrix(college, nrow = T, ncol = N)
ZZ[[3]] = matrix(pmask, nrow = T, ncol = N)
ZZ[[4]] = matrix(pshelter, nrow = T, ncol = N)
ZZ[[5]] = matrix(pgather50, nrow = T, ncol = N)
Y.t = diff(Y) 
D.t = diff(D) 
for(j in 1:length(ZZ)){
  ZZ.t[[j]] = diff(ZZ[[j]])
}
lag = 4

W1 = list()
for(j in 1:lag){
  W1[[j]] = Y.t[(lag-j+1):(nrow(Y.t)-j),]
}
W2 = list()
for(j in 1:2){
  W2[[j]] = D.t[(lag-j+1):(nrow(D.t)-1-j+1),]
}
W2.2 = ZZ.t[[1]][(lag+1):nrow(ZZ.t[[1]]),]
W4 = array(0, dim = c(T-lag-1,N,length(ZZ)-1))
for(j in 1:(length(ZZ)-1)){
  W4[,,j] = ZZ.t[[j+1]][lag:(nrow(ZZ.t[[j+1]])-1),]
}
Q = array(0, dim = c(T,N,T-lag-1))
for(t in (lag+2):T){
  Q[t,,t-lag-1] = rep(1, N)
}
Q.t = array(0, dim = c(T-1,N,T-lag-1))
for(j in 1:(T-lag-1)){
  Q.t[,,j] = diff(Q[,,j])
}
W3 = array(0, dim = c(T-lag-1,N,T-lag-1))
for(j in 1:(T-lag-1)){
  W3[,,j] = Q.t[(lag+1):nrow(Q.t[,,j]),,j]
}

############### AB-LASSO-SS ################
## two folds
Kf = 2
set.seed(202302)
nboot = 100
theta.hat2.all = std.hat2.all = numeric(0)
vcv.all = list()
for(ib in 1:nboot){
  foldid = rep.int(1:Kf, times=ceiling(N/Kf))[sample.int(N)] #fold IDs	
  I = split(1:N, foldid)
  
  Z = list()
  for(t in (lag+2):T){
    Z[[t-lag-1]] = matrix(0, N, (1+lag+dim(W3)[3]+1+length(ZZ)-1))
    y1 = y2.1 = y2.2 = y2 = rep(0, N)
    y3 = matrix(0, N, dim(W3)[3])
    y4 = matrix(0, N, dim(W4)[3])
    zz = matrix(0, N, (length(ZZ)*(t-1)))
    x = matrix(0, N, ((t-2)+(t-1)+dim(zz)[2]))
    for(b in 1:length(I)){
      y1[-I[[b]]] = W1[[1]][t-lag-1,-I[[b]]]  #auxiliary sample
      y2.1[-I[[b]]] = W2[[1]][t-lag-1,-I[[b]]]
      y2[-I[[b]]] = W2.2[t-lag-1,-I[[b]]]
      y3[I[[b]],] = W3[t-lag-1,I[[b]],]       #main sample
      y4[I[[b]],] = W4[t-lag-1,I[[b]],]
      for(j in 1:length(ZZ)){
        zz[I[[b]],((j-1)*(t-1)+1):(j*(t-1))] = t(ZZ[[j]][1:(t-1),I[[b]]])
        zz[-I[[b]],((j-1)*(t-1)+1):(j*(t-1))] = t(ZZ[[j]][1:(t-1),-I[[b]]])
      }
      if(lag==1){if(t>lag+2){x[I[[b]],] = cbind(t(Y[1:(t-2),I[[b]]]), t(D[1:(t-1),I[[b]]]), zz[I[[b]],])}else{x[I[[b]],] = cbind(Y[1:(t-2),I[[b]]], t(D[1:(t-1),I[[b]]]), zz[I[[b]],])}}
      if(lag>1){x[I[[b]],] = cbind(t(Y[1:(t-2),I[[b]]]), t(D[1:(t-1),I[[b]]]), zz[I[[b]],])}
      if(lag==1){if(t>lag+2){x[-I[[b]],] = cbind(t(Y[1:(t-2),-I[[b]]]), t(D[1:(t-1),-I[[b]]]), zz[-I[[b]],])}else{x[-I[[b]],] = cbind(Y[1:(t-2),-I[[b]]], t(D[1:(t-1),-I[[b]]]), zz[-I[[b]],])}}
      if(lag>1){x[-I[[b]],] = cbind(t(Y[1:(t-2),-I[[b]]]), t(D[1:(t-1),-I[[b]]]), zz[-I[[b]],])}
      fit1_s2 = rlasso(x[-I[[b]],], y1[-I[[b]]])
      Z[[t-lag-1]][I[[b]],1] = predict(fit1_s2, x[I[[b]],])
      j = 1
      while(j<lag){
        Z[[t-lag-1]][I[[b]],j+1] = W1[[j+1]][t-lag-1,I[[b]]]
        j = j + 1
      }
      fit2_s2 = rlasso(x[-I[[b]],], y2[-I[[b]]])
      Z[[t-lag-1]][I[[b]],lag+1] = W2[[1]][t-lag-1,I[[b]]]
      Z[[t-lag-1]][I[[b]],lag+1+1] = predict(fit2_s2, x[I[[b]],])
      Z[[t-lag-1]][I[[b]],(lag+2+1):(3+lag+dim(W4)[3]-1)] = y4[I[[b]],]
      Z[[t-lag-1]][I[[b]],(lag+3+dim(W4)[3]):(3+lag+dim(W4)[3]+dim(W3)[3]-1)] = y3[I[[b]],]
    }
  }
  theta.hat2 = matrix(0, (1+lag+dim(W3)[3]+1+length(ZZ)-1), length(I))
  for(b in 1:length(I)){
    sum1 = matrix(0, T-lag-1+lag+1+1+length(ZZ)-1, T-lag-1+lag+1+1+length(ZZ)-1)
    sum2 = rep(0, T-lag-1+lag+1+1+length(ZZ)-1)
    for(i in I[[b]]){
      for(t in (lag+2):T){
        j = 1
        XX = numeric(0)
        while(j<=lag){
          XX = c(XX, W1[[j]][t-lag-1,i])
          j = j + 1
        }
        sum1 = sum1 + Z[[t-lag-1]][i,]%*%t(c(XX, W2[[1]][t-lag-1,i], W2.2[t-lag-1,i], W4[t-lag-1,i,], W3[t-lag-1,i,]))
        sum2 = sum2 + Z[[t-lag-1]][i,]*Y.t[t-1,i]
      }
    }
    sum1.inv = try(solve(sum1))
    theta.hat2[,b] = sum1.inv%*%sum2
  }
  theta.hat2 = rowMeans(theta.hat2)
  
  W4.all = numeric()
  for(j in 1:dim(W4)[3]){
    W4.all = cbind(W4.all, as.vector(W4[,,j]))
  }
  W3.all = numeric()
  for(j in 1:dim(W3)[3]){
    W3.all = cbind(W3.all, as.vector(W3[,,j]))
  }
  j = 1
  XX = numeric(0)
  while(j<=lag){
    XX = cbind(XX, as.vector(W1[[j]]))
    j = j + 1
  }
  vps.t = matrix(as.vector(Y.t[(lag+1):nrow(Y.t),]) - cbind(XX, as.vector(W2[[1]]), as.vector(W2.2), W4.all, W3.all)%*%theta.hat2, nrow = T-lag-1, ncol = N)
  mu = matrix(0, N, T-lag-1+lag+1+1+length(ZZ)-1)
  for(i in 1:N){
    for(t in (lag+2):T){
      mu[i,] = mu[i,] + Z[[t-lag-1]][i,]*vps.t[t-lag-1,i]
    }
    mu[i,] = mu[i,]/(T-lag-1)
  }
  sum3 = sum4 = matrix(0, T-lag-1+lag+1+1+length(ZZ)-1, T-lag-1+lag+1+1+length(ZZ)-1)
  for(i in 1:N){
    for(t in (lag+2):T){
      sum3 = sum3 + (Z[[t-lag-1]][i,]*vps.t[t-lag-1,i] - mu[i,])%*%t(Z[[t-lag-1]][i,]*vps.t[t-lag-1,i] - mu[i,])
    }
  }
  for(i in 1:N){
    for(t in (lag+2):(T-1)){
      sum4 = sum4 + (Z[[t-lag-1]][i,]*vps.t[t-lag-1,i] - mu[i,])%*%t(Z[[t-lag]][i,]*vps.t[t-lag,i] - mu[i,])
    }
  }
  Sigma = sum3 + sum4*(T-lag-2)/(T-lag-1) + t(sum4)*(T-lag-2)/(T-lag-1)
  sum1 = matrix(0, T-lag-1+lag+1+1+length(ZZ)-1, T-lag-1+lag+1+1+length(ZZ)-1)
  for(i in 1:N){
    for(t in (lag+2):T){
      j = 1
      XX = numeric(0)
      while(j<=lag){
        XX = c(XX, W1[[j]][t-lag-1,i])
        j = j + 1
      }
      sum1 = sum1 + Z[[t-lag-1]][i,]%*%t(c(XX, W2[[1]][t-lag-1,i], W2.2[t-lag-1,i], W4[t-lag-1,i,], W3[t-lag-1,i,]))
    }
  }
  sum1.inv = try(solve(sum1))
  std.hat2 = sqrt(diag(sum1.inv%*%Sigma%*%t(sum1.inv)))
  vcv = sum1.inv%*%Sigma%*%t(sum1.inv)
  
  theta.hat2.all = rbind(theta.hat2.all, theta.hat2)
  std.hat2.all = rbind(std.hat2.all, std.hat2)
  vcv.all[[ib]] = vcv
  
  print(paste(ib,"/",nboot))
}
results = list(theta.hat2.all = theta.hat2.all, std.hat2.all = std.hat2.all, vcv.all = vcv.all)
save(results, file = "ablasso_covid_K2.dat")

## five folds
Kf = 5
set.seed(202302)
nboot = 100
theta.hat5.all = std.hat5.all = numeric(0)
vcv.all = list()
for(ib in 1:nboot){
  foldid = rep.int(1:Kf, times=ceiling(N/Kf))[sample.int(N)] #fold IDs	
  I = split(1:N, foldid)
  
  Z = list()
  for(t in (lag+2):T){
    Z[[t-lag-1]] = matrix(0, N, (1+lag+dim(W3)[3]+1+length(ZZ)-1))
    y1 = y2.1 = y2.2 = y2 = rep(0, N)
    y3 = matrix(0, N, dim(W3)[3])
    y4 = matrix(0, N, dim(W4)[3])
    zz = matrix(0, N, (length(ZZ)*(t-1)))
    x = matrix(0, N, ((t-2)+(t-1)+dim(zz)[2]))
    for(b in 1:length(I)){
      y1[-I[[b]]] = W1[[1]][t-lag-1,-I[[b]]]  #auxiliary sample
      y2.1[-I[[b]]] = W2[[1]][t-lag-1,-I[[b]]]
      y2[-I[[b]]] = W2.2[t-lag-1,-I[[b]]]
      y3[I[[b]],] = W3[t-lag-1,I[[b]],]       #main sample
      y4[I[[b]],] = W4[t-lag-1,I[[b]],]
      for(j in 1:length(ZZ)){
        zz[I[[b]],((j-1)*(t-1)+1):(j*(t-1))] = t(ZZ[[j]][1:(t-1),I[[b]]])
        zz[-I[[b]],((j-1)*(t-1)+1):(j*(t-1))] = t(ZZ[[j]][1:(t-1),-I[[b]]])
      }
      if(lag==1){if(t>lag+2){x[I[[b]],] = cbind(t(Y[1:(t-2),I[[b]]]), t(D[1:(t-1),I[[b]]]), zz[I[[b]],])}else{x[I[[b]],] = cbind(Y[1:(t-2),I[[b]]], t(D[1:(t-1),I[[b]]]), zz[I[[b]],])}}
      if(lag>1){x[I[[b]],] = cbind(t(Y[1:(t-2),I[[b]]]), t(D[1:(t-1),I[[b]]]), zz[I[[b]],])}
      if(lag==1){if(t>lag+2){x[-I[[b]],] = cbind(t(Y[1:(t-2),-I[[b]]]), t(D[1:(t-1),-I[[b]]]), zz[-I[[b]],])}else{x[-I[[b]],] = cbind(Y[1:(t-2),-I[[b]]], t(D[1:(t-1),-I[[b]]]), zz[-I[[b]],])}}
      if(lag>1){x[-I[[b]],] = cbind(t(Y[1:(t-2),-I[[b]]]), t(D[1:(t-1),-I[[b]]]), zz[-I[[b]],])}      
      fit1_s2 = rlasso(x[-I[[b]],], y1[-I[[b]]])
      Z[[t-lag-1]][I[[b]],1] = predict(fit1_s2, x[I[[b]],])
      j = 1
      while(j<lag){
        Z[[t-lag-1]][I[[b]],j+1] = W1[[j+1]][t-lag-1,I[[b]]]
        j = j + 1
      }
      fit2_s2 = rlasso(x[-I[[b]],], y2[-I[[b]]])
      Z[[t-lag-1]][I[[b]],lag+1] = W2[[1]][t-lag-1,I[[b]]]
      Z[[t-lag-1]][I[[b]],lag+1+1] = predict(fit2_s2, x[I[[b]],])
      Z[[t-lag-1]][I[[b]],(lag+2+1):(3+lag+dim(W4)[3]-1)] = y4[I[[b]],]
      Z[[t-lag-1]][I[[b]],(lag+3+dim(W4)[3]):(3+lag+dim(W4)[3]+dim(W3)[3]-1)] = y3[I[[b]],]
    }
  }
  theta.hat5 = matrix(0, (1+lag+dim(W3)[3]+1+length(ZZ)-1), length(I))
  for(b in 1:length(I)){
    sum1 = matrix(0, T-lag-1+lag+1+1+length(ZZ)-1, T-lag-1+lag+1+1+length(ZZ)-1)
    sum2 = rep(0, T-lag-1+lag+1+1+length(ZZ)-1)
    for(i in I[[b]]){
      for(t in (lag+2):T){
        j = 1
        XX = numeric(0)
        while(j<=lag){
          XX = c(XX, W1[[j]][t-lag-1,i])
          j = j + 1
        }
        sum1 = sum1 + Z[[t-lag-1]][i,]%*%t(c(XX, W2[[1]][t-lag-1,i], W2.2[t-lag-1,i], W4[t-lag-1,i,], W3[t-lag-1,i,]))
        sum2 = sum2 + Z[[t-lag-1]][i,]*Y.t[t-1,i]
      }
    }
    sum1.inv = try(solve(sum1))
    theta.hat5[,b] = sum1.inv%*%sum2
  }
  theta.hat5 = rowMeans(theta.hat5)
  
  W4.all = numeric()
  for(j in 1:dim(W4)[3]){
    W4.all = cbind(W4.all, as.vector(W4[,,j]))
  }
  W3.all = numeric()
  for(j in 1:dim(W3)[3]){
    W3.all = cbind(W3.all, as.vector(W3[,,j]))
  }
  j = 1
  XX = numeric(0)
  while(j<=lag){
    XX = cbind(XX, as.vector(W1[[j]]))
    j = j + 1
  }
  vps.t = matrix(as.vector(Y.t[(lag+1):nrow(Y.t),]) - cbind(XX, as.vector(W2[[1]]), as.vector(W2.2), W4.all, W3.all)%*%theta.hat5, nrow = T-lag-1, ncol = N)
  mu = matrix(0, N, T-lag-1+lag+1+1+length(ZZ)-1)
  for(i in 1:N){
    for(t in (lag+2):T){
      mu[i,] = mu[i,] + Z[[t-lag-1]][i,]*vps.t[t-lag-1,i]
    }
    mu[i,] = mu[i,]/(T-lag-1)
  }
  sum3 = sum4 = matrix(0, T-lag-1+lag+1+1+length(ZZ)-1, T-lag-1+lag+1+1+length(ZZ)-1)
  for(i in 1:N){
    for(t in (lag+2):T){
      sum3 = sum3 + (Z[[t-lag-1]][i,]*vps.t[t-lag-1,i] - mu[i,])%*%t(Z[[t-lag-1]][i,]*vps.t[t-lag-1,i] - mu[i,])
    }
  }
  for(i in 1:N){
    for(t in (lag+2):(T-1)){
      sum4 = sum4 + (Z[[t-lag-1]][i,]*vps.t[t-lag-1,i] - mu[i,])%*%t(Z[[t-lag]][i,]*vps.t[t-lag,i] - mu[i,])
    }
  }
  Sigma = sum3 + sum4*(T-lag-2)/(T-lag-1) + t(sum4)*(T-lag-2)/(T-lag-1)
  sum1 = matrix(0, T-lag-1+lag+1+1+length(ZZ)-1, T-lag-1+lag+1+1+length(ZZ)-1)
  for(i in 1:N){
    for(t in (lag+2):T){
      j = 1
      XX = numeric(0)
      while(j<=lag){
        XX = c(XX, W1[[j]][t-lag-1,i])
        j = j + 1
      }
      sum1 = sum1 + Z[[t-lag-1]][i,]%*%t(c(XX, W2[[1]][t-lag-1,i], W2.2[t-lag-1,i], W4[t-lag-1,i,], W3[t-lag-1,i,]))
    }
  }
  sum1.inv = try(solve(sum1))
  std.hat5 = sqrt(diag(sum1.inv%*%Sigma%*%t(sum1.inv)))
  vcv = sum1.inv%*%Sigma%*%t(sum1.inv)
  
  theta.hat5.all = rbind(theta.hat5.all, theta.hat5)
  std.hat5.all = rbind(std.hat5.all, std.hat5)
  vcv.all[[ib]] = vcv
  
  print(paste(ib,"/",nboot))
}
results = list(theta.hat5.all = theta.hat5.all, std.hat5.all = std.hat5.all, vcv.all = vcv.all)
save(results, file = "ablasso_covid_K5.dat")

load("ablasso_covid_K2.dat")  
nboot = 100

# short-run effects 
round(colMedians(results$theta.hat2.all)[1:10],2)  
round(colMedians(results$std.hat2.all)[1:10],3)
round(colMedians(results$theta.hat2.all)[1:10]/colMedians(results$std.hat2.all)[1:10], 2) # T-stat

# long-run effects 
lr = cse.lr = rep(0, nboot)
for(ib in 1:nboot){
  coefs  = results$theta.hat2.all[ib,]
  HCV.coefs = results$vcv.all[[ib]]
  lr[ib] = coefs[5]/(1-sum(coefs[1:4]))  
  jac.lr = c(rep(lr[ib],4),1)/(1-sum(coefs[1:4]))
  cse.lr[ib] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,5),c(1:4,5)] %*% jac.lr)  
} 
round(median(lr),2)
round(median(cse.lr),3)
round(median(lr)/median(cse.lr),2) # T-stat

lr = cse.lr = rep(0, nboot)
for(ib in 1:nboot){
  coefs  = results$theta.hat2.all[ib,]
  HCV.coefs = results$vcv.all[[ib]]
  lr[ib] = coefs[7]/(1-sum(coefs[1:4]))  
  jac.lr = c(rep(lr[ib],4),1)/(1-sum(coefs[1:4]))
  cse.lr[ib] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,7),c(1:4,7)] %*% jac.lr)  
} 
round(median(lr),2)
round(median(cse.lr),3)
round(median(lr)/median(cse.lr),2) # T-stat

lr = cse.lr = rep(0, nboot)
for(ib in 1:nboot){
  coefs  = results$theta.hat2.all[ib,]
  HCV.coefs = results$vcv.all[[ib]]
  lr[ib] = coefs[8]/(1-sum(coefs[1:4]))  
  jac.lr = c(rep(lr[ib],4),1)/(1-sum(coefs[1:4]))
  cse.lr[ib] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,8),c(1:4,8)] %*% jac.lr)  
} 
round(median(lr),2)
round(median(cse.lr),3)
round(median(lr)/median(cse.lr),2) # T-stat

lr = cse.lr = rep(0, nboot)
for(ib in 1:nboot){
  coefs  = results$theta.hat2.all[ib,]
  HCV.coefs = results$vcv.all[[ib]]
  lr[ib] = coefs[9]/(1-sum(coefs[1:4]))  
  jac.lr = c(rep(lr[ib],4),1)/(1-sum(coefs[1:4]))
  cse.lr[ib] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,9),c(1:4,9)] %*% jac.lr)  
} 
round(median(lr),2)
round(median(cse.lr),3)
round(median(lr)/median(cse.lr),2) # T-stat

lr = cse.lr = rep(0, nboot)
for(ib in 1:nboot){
  coefs  = results$theta.hat2.all[ib,]
  HCV.coefs = results$vcv.all[[ib]]
  lr[ib] = coefs[10]/(1-sum(coefs[1:4]))  
  jac.lr = c(rep(lr[ib],4),1)/(1-sum(coefs[1:4]))
  cse.lr[ib] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,10),c(1:4,10)] %*% jac.lr)  
} 
round(median(lr),2)
round(median(cse.lr),3)
round(median(lr)/median(cse.lr),2) # T-stat

load("ablasso_covid_K5.dat")
nboot = 100

# short-run effects 
round(colMedians(results$theta.hat5.all)[1:10],2)  
round(colMedians(results$std.hat5.all)[1:10],3)
round(colMedians(results$theta.hat5.all)[1:10]/colMedians(results$std.hat5.all)[1:10], 2) # T-stat

# long-run effects 
lr = cse.lr = rep(0, nboot)
for(ib in 1:nboot){
  coefs  = results$theta.hat5.all[ib,]
  HCV.coefs = results$vcv.all[[ib]]
  lr[ib] = coefs[5]/(1-sum(coefs[1:4]))  
  jac.lr = c(rep(lr[ib],4),1)/(1-sum(coefs[1:4]))
  cse.lr[ib] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,5),c(1:4,5)] %*% jac.lr)  
} 
round(median(lr),2)
round(median(cse.lr),3)
round(median(lr)/median(cse.lr),2) # T-stat

lr = cse.lr = rep(0, nboot)
for(ib in 1:nboot){
  coefs  = results$theta.hat5.all[ib,]
  HCV.coefs = results$vcv.all[[ib]]
  lr[ib] = coefs[7]/(1-sum(coefs[1:4]))  
  jac.lr = c(rep(lr[ib],4),1)/(1-sum(coefs[1:4]))
  cse.lr[ib] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,7),c(1:4,7)] %*% jac.lr)  
} 
round(median(lr),2)
round(median(cse.lr),3)
round(median(lr)/median(cse.lr),2) # T-stat

lr = cse.lr = rep(0, nboot)
for(ib in 1:nboot){
  coefs  = results$theta.hat5.all[ib,]
  HCV.coefs = results$vcv.all[[ib]]
  lr[ib] = coefs[8]/(1-sum(coefs[1:4]))  
  jac.lr = c(rep(lr[ib],4),1)/(1-sum(coefs[1:4]))
  cse.lr[ib] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,8),c(1:4,8)] %*% jac.lr)  
} 
round(median(lr),2)
round(median(cse.lr),3)
round(median(lr)/median(cse.lr),2) # T-stat

lr = cse.lr = rep(0, nboot)
for(ib in 1:nboot){
  coefs  = results$theta.hat5.all[ib,]
  HCV.coefs = results$vcv.all[[ib]]
  lr[ib] = coefs[9]/(1-sum(coefs[1:4]))  
  jac.lr = c(rep(lr[ib],4),1)/(1-sum(coefs[1:4]))
  cse.lr[ib] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,9),c(1:4,9)] %*% jac.lr)  
} 
round(median(lr),2)
round(median(cse.lr),3)
round(median(lr)/median(cse.lr),2) # T-stat

lr = cse.lr = rep(0, nboot)
for(ib in 1:nboot){
  coefs  = results$theta.hat5.all[ib,]
  HCV.coefs = results$vcv.all[[ib]]
  lr[ib] = coefs[10]/(1-sum(coefs[1:4]))  
  jac.lr = c(rep(lr[ib],4),1)/(1-sum(coefs[1:4]))
  cse.lr[ib] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,10),c(1:4,10)] %*% jac.lr)  
} 
round(median(lr),2)
round(median(cse.lr),3)
round(median(lr)/median(cse.lr),2) # T-stat
