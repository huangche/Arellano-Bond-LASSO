setwd("~/Documents/GitHub/Arellano-Bond-LASSO/AB-LASSO_calibration")
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

l = 1
dataset = cbind(as.factor(fips), as.factor(week), as.vector(dlogdc), as.vector(school), as.vector(logdc), as.vector(pmask), as.vector(pgather50), as.vector(college), as.vector(pshelter), as.vector(dlogtests))
colnames(dataset) = c("fips", "week", "dlogdc", "school", "logdc", "pmask", "pgather50", "college", "pshelter", "dlogtests")
# FE
data.fe = pdata.frame(dataset, index = c("fips","week"))
form.fe = logdc ~ lag(logdc, 1:l) + lag(school, 1) + lag(college, 1) + lag(pmask, 1)  + lag(pshelter, 1) + lag(pgather50, 1) + dlogtests - 1
fit.fe = plm(form.fe, data.fe, model = "within", effect = "twoways", index = c("fips","week"))
theta = coef(fit.fe)
residual = fit.fe$residuals
effects = fixef(fit.fe, effect = "twoways", type = "level")

prepara = list(theta = theta, effects = effects, sigma = sd(residual))
save(prepara, file = "application_calibration2.dat")

####### for a single replication
load("application_calibration2.dat")
effects = matrix(prepara$effects, nrow = T-lag, ncol = N)
theta = prepara$theta[c(1:(lag+length(ZZ)),length(prepara$theta))]
sigma2 = prepara$sigma^2
theta.lr = theta[(lag+1):(lag+length(ZZ))]/(1-sum(theta[1:lag])) 

errors = matrix(0, T-lag, N)
for(tt in 1:(T-lag)){
  errors[tt,] = rnorm(N, 0, sigma2)
}
Y = matrix(logdc, nrow = T, ncol = N)
for(t in (lag+1):(T)){
  Y[t,] = theta[1:lag]%*%Y[(t-1):(t-lag),] + theta[(lag+1):(lag+length(ZZ))]%*%rbind(D[t-1,], ZZ[[2]][t-1,], ZZ[[3]][t-1,], ZZ[[4]][t-1,], ZZ[[5]][t-1,]) + tail(theta, 1)*ZZ[[1]][t,] + effects[t-lag,] + errors[t-lag,]
}

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
cat("Starting AB-LASSO...\n")
t1 = Sys.time()
Z = list()
for(t in (lag+1):(T-1)){
  Z[[t-lag]] = matrix(0, N, (1+lag+1+length(ZZ)-1))
  y1 = list()
  y2 = rep(0, N)
  y4 = matrix(0, N, dim(W4)[3])
  zz = matrix(0, N, (length(ZZ)*(t)))
  x = matrix(0, N, ((t-1)+(t)+dim(zz)[2]))
  for(j in 1:lag){
    y1[[j]] = W1[[j]][t-lag,]  
  }
  y2 = W2.2[t-lag,]
  y4= W4[t-lag,,]
  for(j in 1:length(ZZ)){
    zz[,((j-1)*(t)+1):(j*(t))] = t(ZZ[[j]][1:(t),])
  }
  if(lag==1){if(t>lag+1){x = cbind(t(Y[1:(t-1),]), t(D[1:(t),]), zz)}else{x = cbind(Y[1:(t-1),], t(D[1:(t),]), zz)}}
  if(lag>1){x = cbind(t(Y[1:(t-1),]), t(D[1:(t),]), zz)}
  for(j in 1:lag){
    fit1_s2 = rlasso(x, y1[[j]], penalty = list(homoscedastic = "none", lambda.start = 1.5*sqrt(N)*qnorm(1-0.1/(2*dim(x)[2]))))
    Z[[t-lag]][,j] = predict(fit1_s2, x)
  }
  Z[[t-lag]][,lag+1] = W2[t-lag,]
  Z[[t-lag]][,lag+1+1] = y2
  Z[[t-lag]][,(lag+2+1):(3+lag+dim(W4)[3]-1)] = y4
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
se = sqrt(diag(sum1.inv%*%Sigma%*%t(sum1.inv)))[c(1:(lag+1),(lag+3):(lag+3+length(ZZ)-2),lag+2)]
vcv = sum1.inv%*%Sigma%*%t(sum1.inv)

theta.hat = theta.hat[c(1:(lag+1),(lag+3):(lag+3+length(ZZ)-2),lag+2)]
cover = as.numeric(theta>=(theta.hat-qnorm(0.975)*se)&theta<=(theta.hat+qnorm(0.975)*se))
stat = theta.hat/se

lr = cse.lr = rep(0, length(ZZ))
coefs = theta.hat
HCV.coefs = vcv[c(1:(lag+1),(lag+3):(lag+3+length(ZZ)-2)),c(1:(lag+1),(lag+3):(lag+3+length(ZZ)-2))]
for(j in 1:length(ZZ)){
  lr[j] = coefs[j+lag]/(1-sum(coefs[1:lag])) 
  jac.lr = c(rep(lr[j],lag),1)/(1-sum(coefs[1:lag]))
  cse.lr[j] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:lag,j+lag),c(1:lag,j+lag)] %*% jac.lr)
}
cover.lr = as.numeric(theta.lr>=(lr-qnorm(0.975)*cse.lr)&theta.lr<=(lr+qnorm(0.975)*cse.lr))
stat.lr = lr/cse.lr
c1 = difftime(Sys.time(), t1, units = "secs")

############### AB-LASSO-SS (2 folds) ################
cat("Starting AB-LASSO-SS...\n")
t2 = Sys.time()
Kf = 2
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
    y2 = rep(0, N)
    y4 = matrix(0, N, dim(W4)[3])
    zz = matrix(0, N, (length(ZZ)*(t)))
    x = matrix(0, N, ((t-1)+(t)+dim(zz)[2]))
    for(b in 1:length(I)){
      for(j in 1:lag){
        y1[[j]] = rep(0, N)
        y1[[j]][-I[[b]]] = W1[[j]][t-lag,-I[[b]]]  #auxiliary sample
      }
      y2[I[[b]]] = W2.2[t-lag,I[[b]]]     #main sample
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
      Z[[t-lag]][I[[b]],lag+1] = W2[t-lag,I[[b]]]
      Z[[t-lag]][I[[b]],lag+1+1] = y2[I[[b]]]
      Z[[t-lag]][I[[b]],(lag+2+1):(3+lag+dim(W4)[3]-1)] = y4[I[[b]],]
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

theta.hat2 = colMedians(theta.hat2.all)[c(1:(lag+1),(lag+3):(lag+3+length(ZZ)-2),lag+2)]
se2 = colMedians(std.hat2.all)[c(1:(lag+1),(lag+3):(lag+3+length(ZZ)-2),lag+2)]
sd2 = colSds(theta.hat2.all)[c(1:(lag+1),(lag+3):(lag+3+length(ZZ)-2),lag+2)]
cover2 = as.numeric(theta>=(theta.hat2-qnorm(0.975)*se2)&theta<=(theta.hat2+qnorm(0.975)*se2))
stat2 = theta.hat2/se2

lr.all = cse.lr.all = matrix(0, nboot, length(ZZ))
for(ib in 1:nboot){
  coefs  = theta.hat2.all[ib,c(1:(lag+1),(lag+3):(lag+3+length(ZZ)-2))]
  HCV.coefs = vcv.all[[ib]][c(1:(lag+1),(lag+3):(lag+3+length(ZZ)-2)),c(1:(lag+1),(lag+3):(lag+3+length(ZZ)-2))]
  for(j in 1:length(ZZ)){
    lr.all[ib,j] = coefs[j+lag]/(1-sum(coefs[1:lag])) 
    jac.lr = c(rep(lr.all[ib,j],lag),1)/(1-sum(coefs[1:lag]))
    cse.lr.all[ib,j] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:lag,j+lag),c(1:lag,j+lag)] %*% jac.lr)
  }
} 
lr2 = colMedians(lr.all)
cse.lr2 = colMedians(cse.lr.all)
cover.lr2 = as.numeric(theta.lr>=(lr2-qnorm(0.975)*cse.lr2)&theta.lr<=(lr2+qnorm(0.975)*cse.lr2))
stat.lr2 = lr2/cse.lr2
c2 = difftime(Sys.time(), t2, units = "secs")

dataset = cbind(as.factor(fips), as.factor(week), as.vector(dlogdc), as.vector(school), as.vector(Y), as.vector(pmask), as.vector(pgather50), as.vector(college), as.vector(pshelter), as.vector(dlogtests))
colnames(dataset) = c("fips", "week", "dlogdc", "school", "logdc", "pmask", "pgather50", "college", "pshelter", "dlogtests")

############### FE ################
cat("Starting FE...\n")
t3 = Sys.time()
data.fe = pdata.frame(dataset, index = c("fips","week"))
form.fe = logdc ~ lag(logdc, 1:lag) + lag(school, 1) + lag(college, 1) + lag(pmask, 1)  + lag(pshelter, 1) + lag(pgather50, 1) + dlogtests - 1
fit.fe = plm(form.fe, data.fe, model = "within", effect = "twoways", index = c("fips","week"))
fit.fe = summary(fit.fe)
theta.hat.fe = fit.fe$coefficients[,"Estimate"]
HCV.coefs = vcovHC(fit.fe, cluster = 'group')
se.fe = sqrt(diag(HCV.coefs)) #Clustered std errors
cover.fe = as.numeric(theta>=(theta.hat.fe-qnorm(0.975)*se.fe)&theta<=(theta.hat.fe+qnorm(0.975)*se.fe))
stat.fe = theta.hat.fe/se.fe

lr.fe = cse.lr.fe = rep(0, length(ZZ))
coefs = theta.hat.fe
for(j in 1:length(ZZ)){
  lr.fe[j] = coefs[j+lag]/(1-sum(coefs[1:lag])) 
  jac.lr = c(rep(lr.fe[j],lag),1)/(1-sum(coefs[1:lag]))
  cse.lr.fe[j] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:lag,j+lag),c(1:lag,j+lag)] %*% jac.lr)
}
cover.lr.fe = as.numeric(theta.lr>=(lr.fe-qnorm(0.975)*cse.lr.fe)&theta.lr<=(lr.fe+qnorm(0.975)*cse.lr.fe))
stat.lr.fe = lr.fe/cse.lr.fe
c3 = difftime(Sys.time(), t3, units = "secs")

############### DFE-A ################
cat("Starting DFE...\n")
t4 = Sys.time()
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
cover.feabc = as.numeric(theta>=(theta.hat.feabc-qnorm(0.975)*se.feabc)&theta<=(theta.hat.feabc+qnorm(0.975)*se.feabc))
stat.feabc = theta.hat.feabc/se.feabc

lr.feabc = cse.lr.feabc = rep(0, length(ZZ))
coefs = theta.hat.feabc
for(j in 1:length(ZZ)){
  lr.feabc[j] = coefs[j+lag]/(1-sum(coefs[1:lag])) 
  jac.lr = c(rep(lr.feabc[j],lag),1)/(1-sum(coefs[1:lag]))
  cse.lr.feabc[j] = sqrt(t(jac.lr) %*% HCV.coefs[c(1:lag,j+lag),c(1:lag,j+lag)] %*% jac.lr)
}
cover.lr.feabc = as.numeric(theta.lr>=(lr.feabc-qnorm(0.975)*cse.lr.feabc)&theta.lr<=(lr.feabc+qnorm(0.975)*cse.lr.feabc))
stat.lr.feabc = lr.feabc/cse.lr.feabc
c4 = difftime(Sys.time(), t4, units = "secs")

time = c(c1, c2, c3, c4)
result = list(theta.hat = theta.hat, se = se, cover = cover, stat = stat, 
              lr = lr, cse.lr = cse.lr, cover.lr = cover.lr, stat.lr = stat.lr,
              theta.hat2 = theta.hat2, se2 = se2, cover2 = cover2, stat2 = stat2, 
              lr2 = lr2, cse.lr2 = cse.lr2, cover.lr2 = cover.lr2, stat.lr2 = stat.lr2,
              theta.hat.fe = theta.hat.fe, se.fe = se.fe, cover.fe = cover.fe, stat.fe = stat.fe, 
              lr.fe = lr.fe, cse.lr.fe = cse.lr.fe, cover.lr.fe = cover.lr.fe, stat.lr.fe = stat.lr.fe,
              theta.hat.feabc = theta.hat.feabc, se.feabc = se.feabc, cover.feabc = cover.feabc, stat.feabc = stat.feabc, 
              lr.feabc = lr.feabc, cse.lr.feabc = cse.lr.feabc, cover.lr.feabc = cover.lr.feabc, stat.lr.feabc = stat.lr.feabc,
              sd2 = sd2, time = time)

out_dir = paste("application_calibraion_K", Kf, sep="")
if(!dir.exists(out_dir)){dir.create(out_dir)}
filename = file.path(out_dir, paste0("result_", rep_id, ".dat"))
save(result, file = filename)
