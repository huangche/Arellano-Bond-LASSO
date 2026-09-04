setwd("~/Documents/GitHub/Arellano-Bond-LASSO/AB-LASSO_Covid application")
library(mvtnorm)
library(plm)
library(hdm)
library(glmnet)
library(matrixStats)
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

Y.t = D.t = matrix(0, T-1, N)
for(i in 1:N){
  for(t in 1:(T-1)){
    Y.t[t,i] = (Y[t,i] - mean(Y[(t+1):T,i]))*sqrt((T-t)/(T-t+1))
    D.t[t,i] = (D[t,i] - mean(D[(t+1):T,i]))*sqrt((T-t)/(T-t+1))
  }
}
Y.t = Y.t - rowMeans(Y.t)
D.t = D.t - rowMeans(D.t)
for(j in 1:length(ZZ)){
  ZZ.t[[j]] = matrix(0, T-1, N)
  for(i in 1:N){
    for(t in 1:(T-1)){
      ZZ.t[[j]][t,i] = (ZZ[[j]][t,i] - mean(ZZ[[j]][(t+1):T,i]))*sqrt((T-t)/(T-t+1))
    }
  }
  ZZ.t[[j]] = ZZ.t[[j]] - rowMeans(ZZ.t[[j]])
}
lag = 4

W1 = list()
for(j in 1:lag){
  W1[[j]] = Y.t[(lag-j+1):(nrow(Y.t)-j),]
}
W2 = D.t[(lag-1+1):(nrow(D.t)-1),]
W2.2 = ZZ.t[[1]][(lag+1):nrow(ZZ.t[[1]]),]
W4 = array(0, dim = c(T-lag-1,N,length(ZZ)-1))
for(j in 1:(length(ZZ)-1)){
  W4[,,j] = ZZ.t[[j+1]][lag:(nrow(ZZ.t[[j+1]])-1),]
}

############### AB-LASSO ################
Z = list()
for(t in (lag+1):(T-1)){
  Z[[t-lag]] = matrix(0, N, (1+lag+1+length(ZZ)-1))
  y1 = list()
  y2 = y2.2 = rep(0, N)
  y4 = matrix(0, N, dim(W4)[3])
  zz = matrix(0, N, (length(ZZ)*(t)))
  x = matrix(0, N, ((t-1)+(t)+dim(zz)[2]))
  for(j in 1:lag){
    y1[[j]] = W1[[j]][t-lag,]  
  }
  y2 = W2[t-lag,]
  y2.2 = W2.2[t-lag,]
  y4= W4[t-lag,,]
  for(j in 1:length(ZZ)){
    zz[,((j-1)*(t)+1):(j*(t))] = t(ZZ[[j]][1:(t),])
  }
  if(lag==1){if(t>lag+1){x = cbind(t(Y[1:(t-1),]), t(D[1:(t),]), zz)}else{x = cbind(Y[1:(t-1),], t(D[1:(t),]), zz)}}
  if(lag>1){x = cbind(t(Y[1:(t-1),]), t(D[1:(t),]), zz)}
  for(j in 1:lag){
    fit1 = rlasso(x, y1[[j]], penalty = list(homoscedastic = "none", lambda.start = 1.5*sqrt(N)*qnorm(1-0.1/(2*dim(x)[2]))))
    Z[[t-lag]][,j] = predict(fit1, x)
  }
  fit2_1 = rlasso(x, y2, penalty = list(homoscedastic = "none", lambda.start = 1.5*sqrt(N)*qnorm(1-0.1/(2*dim(x)[2]))))
  Z[[t-lag]][,lag+1] = predict(fit2_1, x)
  fit2_2 = rlasso(x, y2.2, penalty = list(homoscedastic = "none", lambda.start = 1.5*sqrt(N)*qnorm(1-0.1/(2*dim(x)[2]))))
  Z[[t-lag]][,lag+1+1] = predict(fit2_2, x)
  for(j in 1:ncol(y4)){
    fit2_3 = rlasso(x, y4[,j], penalty = list(homoscedastic = "none", lambda.start = 1.5*sqrt(N)*qnorm(1-0.1/(2*dim(x)[2]))))
    Z[[t-lag]][,lag+2+j] =  predict(fit2_3, x)
  }
}
sum1 = matrix(0, lag+1+1+length(ZZ)-1, lag+1+1+length(ZZ)-1)
sum2 = rep(0, lag+1+1+length(ZZ)-1)
for(i in 1:N){
  for(t in (lag+1):(T-1)){
    j = 1
    XX = numeric(0)
    while(j<=lag){
      XX = c(XX, W1[[j]][t-lag,i])
      j = j + 1
    }
    sum1 = sum1 + Z[[t-lag]][i,]%*%t(c(XX, W2[t-lag,i], W2.2[t-lag,i], W4[t-lag,i,]))
    sum2 = sum2 + Z[[t-lag]][i,]*Y.t[t,i]
  }
}
sum1.inv = try(solve(sum1))
theta.hat = sum1.inv%*%sum2

W4.all = numeric()
for(j in 1:dim(W4)[3]){
  W4.all = cbind(W4.all, as.vector(W4[,,j]))
}
j = 1
XX = numeric(0)
while(j<=lag){
  XX = cbind(XX, as.vector(W1[[j]]))
  j = j + 1
}
vps.t = matrix(as.vector(Y.t[(lag+1):nrow(Y.t),]) - cbind(XX, as.vector(W2), as.vector(W2.2), W4.all)%*%theta.hat, nrow = T-lag-1, ncol = N)
mu = matrix(0, N, lag+1+1+length(ZZ)-1)
for(i in 1:N){
  for(t in (lag+1):(T-1)){
    mu[i,] = mu[i,] + Z[[t-lag]][i,]*vps.t[t-lag,i]
  }
  mu[i,] = mu[i,]/(T-lag-1)
}
sum3 = sum4 = matrix(0, lag+1+1+length(ZZ)-1, lag+1+1+length(ZZ)-1)
for(i in 1:N){
  for(t in (lag+1):(T-1)){
    sum3 = sum3 + (Z[[t-lag]][i,]*vps.t[t-lag,i] - mu[i,])%*%t(Z[[t-lag]][i,]*vps.t[t-lag,i] - mu[i,])
  }
}
Sigma = sum3 
std.hat = sqrt(diag(sum1.inv%*%Sigma%*%t(sum1.inv)))
vcv = sum1.inv%*%Sigma%*%t(sum1.inv)
results = list(theta.hat = theta.hat, std.hat = std.hat, vcv = vcv)
save(results, file = "ablasso_covid.dat")

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
  
  Y.t = D.t = matrix(0, T-1, N)
  for(i in 1:N){
    for(t in 1:(T-1)){
      Y.t[t,i] = (Y[t,i] - mean(Y[(t+1):T,i]))*sqrt((T-t)/(T-t+1))
      D.t[t,i] = (D[t,i] - mean(D[(t+1):T,i]))*sqrt((T-t)/(T-t+1))
    }
  }
  for(j in 1:length(ZZ)){
    ZZ.t[[j]] = matrix(0, T-1, N)
    for(i in 1:N){
      for(t in 1:(T-1)){
        ZZ.t[[j]][t,i] = (ZZ[[j]][t,i] - mean(ZZ[[j]][(t+1):T,i]))*sqrt((T-t)/(T-t+1))
      }
    }
  }
  for(b in 1:length(I)){
    Y.t[,I[[b]]] = Y.t[,I[[b]]] - rowMeans(Y.t[,I[[b]]])
    D.t[,I[[b]]] = D.t[,I[[b]]] - rowMeans(D.t[,I[[b]]])
    for(j in 1:length(ZZ)){
      ZZ.t[[j]][,I[[b]]] = ZZ.t[[j]][,I[[b]]] - rowMeans(ZZ.t[[j]][,I[[b]]])
    }
  }
  W1 = list()
  for(j in 1:lag){
    W1[[j]] = Y.t[(lag-j+1):(nrow(Y.t)-j),]
  }
  W2 = D.t[(lag-1+1):(nrow(D.t)-1),]
  W2.2 = ZZ.t[[1]][(lag+1):nrow(ZZ.t[[1]]),]
  W4 = array(0, dim = c(T-lag-1,N,length(ZZ)-1))
  for(j in 1:(length(ZZ)-1)){
    W4[,,j] = ZZ.t[[j+1]][lag:(nrow(ZZ.t[[j+1]])-1),]
  }
  
  Z = list()
  for(t in (lag+1):(T-1)){
    Z[[t-lag]] = matrix(0, N, (1+lag+1+length(ZZ)-1))
    y1 = list()
    y2 = y2.2 = rep(0, N)
    y4 = matrix(0, N, dim(W4)[3])
    zz = matrix(0, N, (length(ZZ)*(t)))
    x = matrix(0, N, ((t-1)+(t)+dim(zz)[2]))
    for(b in 1:length(I)){
      for(j in 1:lag){
        y1[[j]] = rep(0, N)
        y1[[j]][-I[[b]]] = W1[[j]][t-lag,-I[[b]]]  #auxiliary sample
      }
      y2[-I[[b]]] = W2[t-lag,-I[[b]]]
      y2.2[-I[[b]]] = W2.2[t-lag,-I[[b]]]
      y4[-I[[b]],] = W4[t-lag,-I[[b]],]
      y2[I[[b]]] = W2[t-lag,I[[b]]]         #main sample
      y2.2[I[[b]]] = W2.2[t-lag,I[[b]]]
      y4[I[[b]],] = W4[t-lag,I[[b]],]
      for(j in 1:length(ZZ)){
        zz[I[[b]],((j-1)*(t)+1):(j*(t))] = t(ZZ[[j]][1:(t),I[[b]]])
        zz[-I[[b]],((j-1)*(t)+1):(j*(t))] = t(ZZ[[j]][1:(t),-I[[b]]])
      }
      if(lag==1){if(t>lag+1){x[I[[b]],] = cbind(t(Y[1:(t-1),I[[b]]]), t(D[1:(t),I[[b]]]), zz[I[[b]],])}else{x[I[[b]],] = cbind(Y[1:(t-1),I[[b]]], t(D[1:(t),I[[b]]]), zz[I[[b]],])}}
      if(lag>1){x[I[[b]],] = cbind(t(Y[1:(t-1),I[[b]]]), t(D[1:(t),I[[b]]]), zz[I[[b]],])}
      if(lag==1){if(t>lag+1){x[-I[[b]],] = cbind(t(Y[1:(t-1),-I[[b]]]), t(D[1:(t),-I[[b]]]), zz[-I[[b]],])}else{x[-I[[b]],] = cbind(Y[1:(t-1),-I[[b]]], t(D[1:(t),-I[[b]]]), zz[-I[[b]],])}}
      if(lag>1){x[-I[[b]],] = cbind(t(Y[1:(t-1),-I[[b]]]), t(D[1:(t),-I[[b]]]), zz[-I[[b]],])}
      for(j in 1:lag){
        fit1_s2 = rlasso(x[-I[[b]],], y1[[j]][-I[[b]]], penalty = list(homoscedastic = "none", lambda.start = 1.5*sqrt(N-length(I[[b]]))*qnorm(1-0.1/(2*dim(x)[2]))))
        Z[[t-lag]][I[[b]],j] = predict(fit1_s2, x[I[[b]],])
      }
      fit2_1_s2 = rlasso(x[-I[[b]],], y2[-I[[b]]], penalty = list(homoscedastic = "none", lambda.start = 1.5*sqrt(N-length(I[[b]]))*qnorm(1-0.1/(2*dim(x)[2]))))
      Z[[t-lag]][I[[b]],lag+1] = predict(fit2_1_s2, x[I[[b]],])
      fit2_2_s2 = rlasso(x[-I[[b]],], y2.2[-I[[b]]], penalty = list(homoscedastic = "none", lambda.start = 1.5*sqrt(N-length(I[[b]]))*qnorm(1-0.1/(2*dim(x)[2]))))
      Z[[t-lag]][I[[b]],lag+1+1] = predict(fit2_2_s2, x[I[[b]],])
      for(j in 1:ncol(y4)){
        fit2_3_s2 = rlasso(x[-I[[b]],], y4[-I[[b]],j], penalty = list(homoscedastic = "none", lambda.start = 1.5*sqrt(N-length(I[[b]]))*qnorm(1-0.1/(2*dim(x)[2]))))
        Z[[t-lag]][I[[b]],lag+2+j] =  predict(fit2_3_s2, x[I[[b]],])
      }
    }
  }
  theta.hat2 = matrix(0, (1+lag+1+length(ZZ)-1), length(I))
  for(b in 1:length(I)){
    sum1 = matrix(0, lag+1+1+length(ZZ)-1, lag+1+1+length(ZZ)-1)
    sum2 = rep(0, lag+1+1+length(ZZ)-1)
    for(i in I[[b]]){
      for(t in (lag+1):(T-1)){
        j = 1
        XX = numeric(0)
        while(j<=lag){
          XX = c(XX, W1[[j]][t-lag,i])
          j = j + 1
        }
        sum1 = sum1 + Z[[t-lag]][i,]%*%t(c(XX, W2[t-lag,i], W2.2[t-lag,i], W4[t-lag,i,]))
        sum2 = sum2 + Z[[t-lag]][i,]*Y.t[t,i]
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
  j = 1
  XX = numeric(0)
  while(j<=lag){
    XX = cbind(XX, as.vector(W1[[j]]))
    j = j + 1
  }
  vps.t = matrix(as.vector(Y.t[(lag+1):nrow(Y.t),]) - cbind(XX, as.vector(W2), as.vector(W2.2), W4.all)%*%theta.hat2, nrow = T-lag-1, ncol = N)
  mu = matrix(0, N, lag+1+1+length(ZZ)-1)
  for(i in 1:N){
    for(t in (lag+1):(T-1)){
      mu[i,] = mu[i,] + Z[[t-lag]][i,]*vps.t[t-lag,i]
    }
    mu[i,] = mu[i,]/(T-lag-1)
  }
  sum3 = sum4 = matrix(0, lag+1+1+length(ZZ)-1, lag+1+1+length(ZZ)-1)
  for(i in 1:N){
    for(t in (lag+1):(T-1)){
      sum3 = sum3 + (Z[[t-lag]][i,]*vps.t[t-lag,i] - mu[i,])%*%t(Z[[t-lag]][i,]*vps.t[t-lag,i] - mu[i,])
    }
  }
  Sigma = sum3 
  sum1 = matrix(0, lag+1+1+length(ZZ)-1, lag+1+1+length(ZZ)-1)
  for(i in 1:N){
    for(t in (lag+1):(T-1)){
      j = 1
      XX = numeric(0)
      while(j<=lag){
        XX = c(XX, W1[[j]][t-lag,i])
        j = j + 1
      }
      sum1 = sum1 + Z[[t-lag]][i,]%*%t(c(XX, W2[t-lag,i], W2.2[t-lag,i], W4[t-lag,i,]))
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

load("ablasso_covid.dat") 

# short-run effects 
round(results$theta.hat[1:10],2)  
round(results$std.hat[1:10],3)
round(results$theta.hat[1:10]/results$std.hat[1:10],2) # T-stat

# long-run effects
coefs  = results$theta.hat
HCV.coefs = results$vcv
lr = coefs[5]/(1-sum(coefs[1:4]))  
jac.lr = c(rep(lr,4),1)/(1-sum(coefs[1:4]))
cse.lr = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,5),c(1:4,5)] %*% jac.lr)  
round(lr,2)
round(cse.lr,3)
round(lr/cse.lr,2) # T-stat

lr = coefs[7]/(1-sum(coefs[1:4]))  
jac.lr = c(rep(lr,4),1)/(1-sum(coefs[1:4]))
cse.lr = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,7),c(1:4,7)] %*% jac.lr)  
round(lr,2)
round(cse.lr,3)
round(lr/cse.lr,2) # T-stat

lr = coefs[8]/(1-sum(coefs[1:4]))  
jac.lr = c(rep(lr,4),1)/(1-sum(coefs[1:4]))
cse.lr = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,8),c(1:4,8)] %*% jac.lr)  
round(lr,2)
round(cse.lr,3)
round(lr/cse.lr,2) # T-stat

lr = coefs[9]/(1-sum(coefs[1:4]))  
jac.lr = c(rep(lr,4),1)/(1-sum(coefs[1:4]))
cse.lr = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,9),c(1:4,9)] %*% jac.lr)  
round(lr,2)
round(cse.lr,3)
round(lr/cse.lr,2) # T-stat

lr = coefs[10]/(1-sum(coefs[1:4]))  
jac.lr = c(rep(lr,4),1)/(1-sum(coefs[1:4]))
cse.lr = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,10),c(1:4,10)] %*% jac.lr)  
round(lr,2)
round(cse.lr,3)
round(lr/cse.lr,2) # T-stat

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

############### DFE-A ################
dataset = cbind(as.factor(fips), as.factor(week), as.vector(dlogdc), as.vector(school), as.vector(Y), as.vector(pmask), as.vector(pgather50), as.vector(college), as.vector(pshelter), as.vector(dlogtests))
colnames(dataset) = c("fips", "week", "dlogdc", "school", "logdc", "pmask", "pgather50", "college", "pshelter", "dlogtests")

data.fe = pdata.frame(dataset, index = c("fips","week"))
form.fe = logdc ~ lag(logdc, 1:lag) + lag(school, 1) + lag(college, 1) + lag(pmask, 1)  + lag(pshelter, 1) + lag(pgather50, 1) + dlogtests - 1
fit.fe = plm(form.fe, data.fe, model = "within", effect = "twoways", index = c("fips","week"))
fit.fe = summary(fit.fe)
theta.hat.fe = fit.fe$coefficients[,"Estimate"]
HCV.coefs = vcovHC(fit.fe, cluster = 'group')
se.fe = sqrt(diag(HCV.coefs)) 

form = character(0)
for(j in 1:lag){
  name = paste("logdc.lag", j, sep="")
  data.fe[[name]] = lag(data.fe$logdc, j)
  form = paste(form, name, " + ", sep="")
}
data.fe$school.lag1 = lag(data.fe$school, 1)
data.fe$college.lag1 = lag(data.fe$college, 1)
data.fe$pmask.lag1 = lag(data.fe$pmask, 1)
data.fe$pshelter.lag1 = lag(data.fe$pshelter, 1)
data.fe$pgather50.lag1 = lag(data.fe$pgather50, 1)
form.feabc = paste("logdc ~ ", form, " school.lag1 + college.lag1 + pmask.lag1 + pshelter.lag1 + pgather50.lag1 + dlogtests + factor(fips) + factor(week)", sep="")
fit.feabc = lm(form.feabc, data.fe, x = TRUE, na.action = na.omit)
res.feabc = fit.feabc$residuals
jac = solve(t(fit.feabc$x)%*%fit.feabc$x/length(res.feabc))[2:(1+lag+length(ZZ)+1),2:(1+lag+length(ZZ)+1)]
indexes = c(1:length(res.feabc))
indexes = indexes[-c(1+c(0:(N-1))*length(res.feabc)/N)]
bscore = t(fit.feabc$x[indexes, 2:(1+lag+length(ZZ)+1)])%*%res.feabc[indexes-1]/length(indexes)
bias = -jac%*%bscore*N/length(res.feabc)
theta.hat.feabc = theta.hat.fe - bias
se.feabc = se.fe

results = list(theta.hat.feabc = theta.hat.feabc, se.feabc = se.feabc, HCV.coefs = HCV.coefs)
save(results, file = "dfe_covid.dat")

load("dfe_covid.dat") 

# short-run effects 
round(results$theta.hat.feabc[1:10],2)  
round(results$se.feabc[1:10],3)
round(results$theta.hat.feabc[1:10]/results$se.feabc[1:10],2) # T-stat

# long-run effects
coefs  = results$theta.hat.feabc
HCV.coefs = results$HCV.coefs
lr = coefs[5]/(1-sum(coefs[1:4]))  
jac.lr = c(rep(lr,4),1)/(1-sum(coefs[1:4]))
cse.lr = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,5),c(1:4,5)] %*% jac.lr)  
round(lr,2)
round(cse.lr,3)
round(lr/cse.lr,2) # T-stat

lr = coefs[6]/(1-sum(coefs[1:4]))  
jac.lr = c(rep(lr,4),1)/(1-sum(coefs[1:4]))
cse.lr = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,6),c(1:4,6)] %*% jac.lr)  
round(lr,2)
round(cse.lr,3)
round(lr/cse.lr,2) # T-stat

lr = coefs[7]/(1-sum(coefs[1:4]))  
jac.lr = c(rep(lr,4),1)/(1-sum(coefs[1:4]))
cse.lr = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,7),c(1:4,7)] %*% jac.lr)  
round(lr,2)
round(cse.lr,3)
round(lr/cse.lr,2) # T-stat

lr = coefs[8]/(1-sum(coefs[1:4]))  
jac.lr = c(rep(lr,4),1)/(1-sum(coefs[1:4]))
cse.lr = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,8),c(1:4,8)] %*% jac.lr)  
round(lr,2)
round(cse.lr,3)
round(lr/cse.lr,2) # T-stat

lr = coefs[9]/(1-sum(coefs[1:4]))  
jac.lr = c(rep(lr,4),1)/(1-sum(coefs[1:4]))
cse.lr = sqrt(t(jac.lr) %*% HCV.coefs[c(1:4,9),c(1:4,9)] %*% jac.lr)  
round(lr,2)
round(cse.lr,3)
round(lr/cse.lr,2) # T-stat
