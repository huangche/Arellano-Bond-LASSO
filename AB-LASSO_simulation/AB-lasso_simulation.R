library(mvtnorm)
library(plm)
library(hdm)
library(glmnet)
library(matrixStats)
rep = 500
NN = c(200,300) 
TT = seq(30,50,10) 
sigma_alpha = sigma_gamma = sigma_eps.d = sigma_eps.y = 1
cov_eps = 0.5
theta = c(0.8,1)
rho = 0.5
nboot = 100
S2 = 1
data = list()
ii = 1
for(nn in 1:length(NN)){
  for(tt in 1:length(TT)){
    N = NN[nn]
    T = TT[tt]
    alpha = matrix(0, nrow = rep*2, ncol = N)
    gamma = matrix(0, nrow = rep*2, ncol = T+10)
    eps.y = eps.d = array(0, dim = c(rep*2, T+10, N))
    eps = array(0, dim = c(rep*2, T+10, N, 2))
    for(r in 1:(rep*2)){
      alpha[r,] = rnorm(N, 0, sigma_alpha)
      gamma[r,] = rnorm(T+10, 0, sigma_gamma)
      eps.y[r,,] = rmvnorm(T+10, rep(0,N), diag(N)*sigma_eps.y)
      eps.d[r,,] = rmvnorm(T+10, rep(0,N), diag(N)*sigma_eps.d)
      for(i in 1:N){
        eps[r,,i,] = rmvnorm(T+10, rep(0,2), cbind(c(sigma_eps.y, cov_eps), c(cov_eps, sigma_eps.d)))
      }
    }
    data[[ii]] = list(alpha = alpha, gamma = gamma, eps.y = eps.y, eps.d = eps.d, eps = eps)
    ii = ii + 1
  }
}
save(data, file = "datasample.dat")
load("datasample.dat")
results = numeric(0)
output = list()
ii = 1
#options(warn=-1)
for(nn in 1:length(NN)){
  for(tt in 1:length(TT)){
    N = NN[nn]
    T = TT[tt]
    theta.hat = theta.hat2 = theta.hat.ab = theta.hat5 = theta.hat.dab = matrix(0, rep, 2)
    std.hat = std.hat2 = std.hat.ab = std.hat5 = std.hat.dab = matrix(0, rep, 2)
    sd2 = sd5 = matrix(0, rep, 2)
    cover = cover2 = cover.ab = cover5 = cover.dab = matrix(0, rep, 2)
    stat = stat2 = stat.ab = stat5 = stat.dab = matrix(0, rep, 2)
    c1 = c2 = c3 = c4 = c5 = rep(0, rep)
    r = 1
    r.data = 0
    while(r<=rep){
      r.data = r.data + 1
      print(paste(N,"-",T,":",r.data))
      alpha = data[[ii]]$alpha[r.data,]
      gamma = data[[ii]]$gamma[r.data,]
      eps.y = data[[ii]]$eps.y[r.data,,]
      eps.d = data[[ii]]$eps.d[r.data,,]
      eps = data[[ii]]$eps[r.data,,,]
      for(i in 1:N){
        eps.y[,i] = eps[,i,1] 
        eps.d[2:(T+10),i] = eps[1:(T+10-1),i,2] 
      }
      Y = D = matrix(0, nrow = T+10, ncol = N)
      for(t in 2:(T+10)){
        D[t,] = rho*D[t-1,] + eps.d[t,]
        Y[t,] = theta[1]*Y[t-1,] + theta[2]*D[t,] + alpha + eps.y[t,] + gamma[t]
      }
      Y = Y[-(1:10),]
      D = D[-(1:10),]
      ############### AB-LASSO ################
      ## without sample splitting
      t1 = Sys.time()
      Y.t = diff(Y) - rowMeans(diff(Y))
      D.t = diff(D) - rowMeans(diff(D))
      W1 = Y.t[1:(nrow(Y.t)-1),]
      W2 = D.t[2:nrow(D.t),]
      Z = list()
      for(t in 3:T){
        y1 = W1[t-2,]
        y2 = W2[t-2,]
        if(t>3){x = cbind(t(Y[1:(t-2),]), t(D[1:(t-1),]))}else{x = cbind(Y[1:(t-2),], t(D[1:(t-1),]))}
        fit1 = rlasso(x, y1)
        fit2 = rlasso(x, y2)
        Z[[t-2]] = cbind(predict(fit1), predict(fit2))
      }
      sum1 = matrix(0, 2, 2)
      sum2 = rep(0, 2)
      for(i in 1:N){
        for(t in 3:T){
          sum1 = sum1 + Z[[t-2]][i,]%*%t(c(W1[t-2,i],W2[t-2,i]))
          sum2 = sum2 + Z[[t-2]][i,]*Y.t[t-1,i]
        }
      }
      theta.hat[r,] = solve(sum1)%*%sum2
      vps.t = matrix(as.vector(Y.t[2:nrow(Y.t),]) - cbind(as.vector(W1),as.vector(W2))%*%theta.hat[r,], nrow = T-2, ncol = N)
      mu = matrix(0, N, 2)
      for(i in 1:N){
        for(t in 3:T){
          mu[i,] = mu[i,] + Z[[t-2]][i,]*vps.t[t-2,i]
        }
        mu[i,] = mu[i,]/(T-2)
      }
      sum3 = sum4 = matrix(0, 2, 2)
      for(i in 1:N){
        for(t in 3:T){
          sum3 = sum3 + (Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])%*%t(Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])
        }
      }
      for(i in 1:N){
        for(t in 3:(T-1)){
          sum4 = sum4 + (Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])%*%t(Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])
        }
      }
      Sigma = sum3 + sum4*(T-3)/(T-2) + t(sum4)*(T-3)/(T-2)
      std.hat[r,] = sqrt(diag(solve(sum1)%*%Sigma%*%t(solve(sum1))*N*(T-2)))
      cover[r,1] = as.numeric(theta[1]>=(theta.hat[r,1]-qnorm(0.975)*std.hat[r,1]/sqrt(N*(T-2)))&&theta[1]<=(theta.hat[r,1]+qnorm(0.975)*std.hat[r,1]/sqrt(N*(T-2))))
      cover[r,2] = as.numeric(theta[2]>=(theta.hat[r,2]-qnorm(0.975)*std.hat[r,2]/sqrt(N*(T-2)))&&theta[2]<=(theta.hat[r,2]+qnorm(0.975)*std.hat[r,2]/sqrt(N*(T-2))))
      stat[r,] = sqrt(N*(T-2))*theta.hat[r,]/std.hat[r,]
      c1[r] = difftime(Sys.time(), t1, units = "secs")
      
      ## 2 folds random sample splitting
      t2 = Sys.time()
      Kf = 2
      set.seed(202304)
      theta.hat.ss = std.hat.ss = numeric(0)
      for(ib in 1:nboot){
        foldid = rep.int(1:Kf, times=ceiling(N/Kf))[sample.int(N)] #fold IDs
        I = split(1:N, foldid)
        
        Y.t = diff(Y) 
        D.t = diff(D)
        for(b in 1:length(I)){
          Y.t[,I[[b]]] = Y.t[,I[[b]]] - rowMeans(Y.t[,I[[b]]]) 
          D.t[,I[[b]]] = D.t[,I[[b]]] - rowMeans(D.t[,I[[b]]]) 
        }
        W1 = Y.t[1:(nrow(Y.t)-1),]
        W2 = D.t[2:nrow(D.t),]
        
        Z = list()
        for(t in 3:T){
          Z[[t-2]] = matrix(0, N, 2)
          y1 = y2 = rep(0, N)
          x = matrix(0, N, ((t-2)+(t-1)))
          for(b in 1:length(I)){
            y1[-I[[b]]] = W1[t-2,-I[[b]]]  #s2 - auxiliary sample
            y2[-I[[b]]] = W2[t-2,-I[[b]]]
            if(t>3){x[I[[b]],] = cbind(t(Y[1:(t-2),I[[b]]]), t(D[1:(t-1),I[[b]]]))}else{x[I[[b]],]  = cbind(Y[1:(t-2),I[[b]]], t(D[1:(t-1),I[[b]]]))}
            if(t>3){x[-I[[b]],] = cbind(t(Y[1:(t-2),-I[[b]]]), t(D[1:(t-1),-I[[b]]]))}else{x[-I[[b]],] = cbind(Y[1:(t-2),-I[[b]]], t(D[1:(t-1),-I[[b]]]))}
            fit1_s2 = rlasso(x[-I[[b]],], y1[-I[[b]]])
            fit2_s2 = rlasso(x[-I[[b]],], y2[-I[[b]]])
            Z[[t-2]][I[[b]],] = cbind(predict(fit1_s2, x[I[[b]],]), predict(fit2_s2, x[I[[b]],])) 
          }
        }
        theta.est = matrix(0, 2, length(I))
        for(b in 1:length(I)){
          sum1 = matrix(0, 2, 2)
          sum2 = rep(0, 2)
          for(i in I[[b]]){
            for(t in 3:T){
              sum1 = sum1 + Z[[t-2]][i,]%*%t(c(W1[t-2,i],W2[t-2,i]))
              sum2 = sum2 + Z[[t-2]][i,]*Y.t[t-1,i]
            }
          }
          theta.est[,b] = solve(sum1)%*%sum2
        }
        theta.hat.ss = rbind(theta.hat.ss, rowMeans(theta.est))
        vps.t = matrix(as.vector(Y.t[2:nrow(Y.t),]) - cbind(as.vector(W1[,]),as.vector(W2[,]))%*%rowMeans(theta.est), nrow = T-2, ncol = N)
        mu = matrix(0, N, 2)
        for(i in 1:N){
          for(t in 3:T){
            mu[i,] = mu[i,] + Z[[t-2]][i,]*vps.t[t-2,i]
          }
          mu[i,] = mu[i,]/(T-2)
        }
        sum3 = sum4 = matrix(0, 2, 2)
        for(i in 1:N){
          for(t in 3:T){
            sum3 = sum3 + (Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])%*%t(Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])
          }
        }
        for(i in 1:N){
          for(t in 3:(T-1)){
            sum4 = sum4 + (Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])%*%t(Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])
          }
        }
        Sigma = sum3 + sum4*(T-3)/(T-2) + t(sum4)*(T-3)/(T-2)
        sum1 = matrix(0, 2, 2)
        for(i in 1:N){
          for(t in 3:T){
            sum1 = sum1 + Z[[t-2]][i,]%*%t(c(W1[t-2,i],W2[t-2,i]))
          }
        }
        std.hat.ss = rbind(std.hat.ss, sqrt(diag(solve(sum1)%*%Sigma%*%t(solve(sum1))*N*(T-2))))
      }
      theta.hat2[r,] = colMedians(theta.hat.ss)
      std.hat2[r,] = colMedians(std.hat.ss)
      sd2[r,] = colSds(theta.hat.ss)
      cover2[r,1] = as.numeric(theta[1]>=(theta.hat2[r,1]-qnorm(0.975)*std.hat2[r,1]/sqrt(N*(T-2)))&&theta[1]<=(theta.hat2[r,1]+qnorm(0.975)*std.hat2[r,1]/sqrt(N*(T-2))))
      cover2[r,2] = as.numeric(theta[2]>=(theta.hat2[r,2]-qnorm(0.975)*std.hat2[r,2]/sqrt(N*(T-2)))&&theta[2]<=(theta.hat2[r,2]+qnorm(0.975)*std.hat2[r,2]/sqrt(N*(T-2))))
      stat2[r,] = sqrt(N*(T-2))*theta.hat2[r,]/std.hat2[r,]
      c2[r] = difftime(Sys.time(), t2, units = "secs")
      
      ## 5 folds random sample splitting
      t3 = Sys.time()
      Kf = 5
      set.seed(202304)
      theta.hat.ss = std.hat.ss = numeric(0)
      for(ib in 1:nboot){
        foldid = rep.int(1:Kf, times=ceiling(N/Kf))[sample.int(N)] #fold IDs
        I = split(1:N, foldid)
        
        W1.all = W2.all = Y.t.all = list()
        for(b in 1:length(I)){
          Y.t = diff(Y) 
          D.t = diff(D)
          Y.t[,I[[b]]] = Y.t[,I[[b]]] - rowMeans(Y.t[,I[[b]]]) 
          D.t[,I[[b]]] = D.t[,I[[b]]] - rowMeans(D.t[,I[[b]]])
          Y.t[,-I[[b]]] = Y.t[,-I[[b]]] - rowMeans(Y.t[,-I[[b]]]) 
          D.t[,-I[[b]]] = D.t[,-I[[b]]] - rowMeans(D.t[,-I[[b]]])
          W1 = Y.t[1:(nrow(Y.t)-1),]
          W2 = D.t[2:nrow(D.t),]
          W1.all[[b]] = W1
          W2.all[[b]] = W2
          Y.t.all[[b]] = Y.t
        }
        
        Z = list()
        for(t in 3:T){
          Z[[t-2]] = matrix(0, N, 2)
          y1 = y2 = rep(0, N)
          x = matrix(0, N, ((t-2)+(t-1)))
          for(b in 1:length(I)){
            W1 = W1.all[[b]] 
            W2 = W2.all[[b]] 
            Y.t = Y.t.all[[b]]
            y1[-I[[b]]] = W1[t-2,-I[[b]]]  #s2 - auxiliary sample
            y2[-I[[b]]] = W2[t-2,-I[[b]]]
            if(t>3){x[I[[b]],] = cbind(t(Y[1:(t-2),I[[b]]]), t(D[1:(t-1),I[[b]]]))}else{x[I[[b]],] = cbind(Y[1:(t-2),I[[b]]], t(D[1:(t-1),I[[b]]]))}
            if(t>3){x[-I[[b]],] = cbind(t(Y[1:(t-2),-I[[b]]]), t(D[1:(t-1),-I[[b]]]))}else{x[-I[[b]],] = cbind(Y[1:(t-2),-I[[b]]], t(D[1:(t-1),-I[[b]]]))}
            fit1_s2 = rlasso(x[-I[[b]],], y1[-I[[b]]])
            fit2_s2 = rlasso(x[-I[[b]],], y2[-I[[b]]])
            Z[[t-2]][I[[b]],] = cbind(predict(fit1_s2, x[I[[b]],]), predict(fit2_s2, x[I[[b]],])) 
          }
        }
        theta.est = matrix(0, 2, length(I))
        for(b in 1:length(I)){
          W1 = W1.all[[b]] 
          W2 = W2.all[[b]] 
          Y.t = Y.t.all[[b]]
          sum1 = matrix(0, 2, 2)
          sum2 = rep(0, 2)
          for(i in I[[b]]){
            for(t in 3:T){
              sum1 = sum1 + Z[[t-2]][i,]%*%t(c(W1[t-2,i],W2[t-2,i]))
              sum2 = sum2 + Z[[t-2]][i,]*Y.t[t-1,i]
            }
          }
          theta.est[,b] = solve(sum1)%*%sum2
        }
        theta.hat.ss = rbind(theta.hat.ss, rowMeans(theta.est))
        vps.t = matrix(as.vector(Y.t[2:nrow(Y.t),]) - cbind(as.vector(W1[,]),as.vector(W2[,]))%*%rowMeans(theta.est), nrow = T-2, ncol = N)
        mu = matrix(0, N, 2)
        for(i in 1:N){
          for(t in 3:T){
            mu[i,] = mu[i,] + Z[[t-2]][i,]*vps.t[t-2,i]
          }
          mu[i,] = mu[i,]/(T-2)
        }
        sum3 = sum4 = matrix(0, 2, 2)
        for(i in 1:N){
          for(t in 3:T){
            sum3 = sum3 + (Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])%*%t(Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])
          }
        }
        for(i in 1:N){
          for(t in 3:(T-1)){
            sum4 = sum4 + (Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])%*%t(Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])
          }
        }
        Sigma = sum3 + sum4*(T-3)/(T-2) + t(sum4)*(T-3)/(T-2)
        sum1 = matrix(0, 2, 2)
        for(i in 1:N){
          for(t in 3:T){
            sum1 = sum1 + Z[[t-2]][i,]%*%t(c(W1[t-2,i],W2[t-2,i]))
          }
        }
        std.hat.ss = rbind(std.hat.ss, sqrt(diag(solve(sum1)%*%Sigma%*%t(solve(sum1))*N*(T-2))))
      }
      theta.hat5[r,] = colMedians(theta.hat.ss)
      std.hat5[r,] = colMedians(std.hat.ss)
      sd5[r,] = colSds(theta.hat.ss)
      cover5[r,1] = as.numeric(theta[1]>=(theta.hat5[r,1]-qnorm(0.975)*std.hat5[r,1]/sqrt(N*(T-2)))&&theta[1]<=(theta.hat5[r,1]+qnorm(0.975)*std.hat5[r,1]/sqrt(N*(T-2))))
      cover5[r,2] = as.numeric(theta[2]>=(theta.hat5[r,2]-qnorm(0.975)*std.hat5[r,2]/sqrt(N*(T-2)))&&theta[2]<=(theta.hat5[r,2]+qnorm(0.975)*std.hat5[r,2]/sqrt(N*(T-2))))
      stat5[r,] = sqrt(N*(T-2))*theta.hat5[r,]/std.hat5[r,]
      c3[r] = difftime(Sys.time(), t3, units = "secs")
      
      ############### AB ################
      index.i = index.t = numeric(0)
      for(i in 1:N){
        index.i = c(index.i, rep(i,T))
        index.t = c(index.t, c(2000:(2000+T-1)))
      }
      dataset.ab = cbind(as.factor(index.i), as.factor(index.t), as.vector(Y), as.vector(D))
      colnames(dataset.ab) = c("firm", "year", "y.ab", "d.ab")
      ## two-step AB 
      t4 = Sys.time()
      fit.ab = try(pgmm(y.ab ~ lag(y.ab,1) + d.ab | lag(y.ab,2:9999) + lag(d.ab,1:9999), dataset.ab, model = "twosteps", effect = "twoways"))#"individual")
      if(class(fit.ab)[1]!="try-error"){
        try.ab = try(summary(fit.ab), silent = TRUE)
        c4[r] = difftime(Sys.time(), t4, units = "secs")
      }
      
      ## debiased two-step AB with 1 partition
      if(class(fit.ab)[1]!="try-error"&&class(try.ab)[1]!="try-error"){
        t5 = Sys.time()
        coeff.dab = se.dab = 0*fit.ab$coefficients[[2]][1:2]
        dab.fit1 = dab.fit2 = try.dab1 = try.dab2 = list()
        s = 1
        set.seed(202304)
        while(s<=S2){
          sample1 = sample(N, ceiling(N/2), replace = FALSE)
          dab.fit1[[s]] = try(pgmm(y.ab ~ lag(y.ab,1) + d.ab | lag(y.ab,2:9999) + lag(d.ab,1:9999), dataset.ab[as.double(index.i) %in% sample1, ], model = "twosteps", effect = "twoways"))#"individual")
          dab.fit2[[s]] = try(pgmm(y.ab ~ lag(y.ab,1) + d.ab | lag(y.ab,2:9999) + lag(d.ab,1:9999), dataset.ab[!as.double(index.i) %in% sample1, ], model = "twosteps", effect = "twoways"))#"individual")
          if(class(dab.fit1[[s]])[1]!="try-error"&&class(dab.fit2[[s]])[1]!="try-error"){
            try.dab1[[s]] = try(summary(dab.fit1[[s]]), silent = TRUE)
            try.dab2[[s]] = try(summary(dab.fit2[[s]]), silent = TRUE)
            if(class(try.dab1[[s]])[1]!="try-error"&&class(try.dab2[[s]])[1]!="try-error"){
              coeff.dab = coeff.dab + ((dab.fit1[[s]]$coefficients[[2]][1:2] + dab.fit2[[s]]$coefficients[[2]][1:2])/2)/S2
              s = s + 1
            }else{s = S2 + 2}
          }else{s = S2 + 2}
        }
        c5[r] = difftime(Sys.time(), t5, units = "secs")
        if(s==(S2+1)){
          theta.hat.ab[r,] = fit.ab$coefficients[[2]][1:2]
          std.hat.ab[r,] = try.ab$coefficients[,"Std. Error"] #summary(fit.ab)$coefficients[,"Std. Error"]
          cover.ab[r,1] = as.numeric(theta[1]>=(theta.hat.ab[r,1]-qnorm(0.975)*std.hat.ab[r,1])&&theta[1]<=(theta.hat.ab[r,1]+qnorm(0.975)*std.hat.ab[r,1]))
          cover.ab[r,2] = as.numeric(theta[2]>=(theta.hat.ab[r,2]-qnorm(0.975)*std.hat.ab[r,2])&&theta[2]<=(theta.hat.ab[r,2]+qnorm(0.975)*std.hat.ab[r,2]))
          stat.ab[r,] = theta.hat.ab[r,]/std.hat.ab[r,]
          theta.hat.dab[r,] = 2*fit.ab$coefficients[[2]][1:2] - coeff.dab
          std.hat.dab[r,] = std.hat.ab[r,]
          cover.dab[r,1] = as.numeric(theta[1]>=(theta.hat.dab[r,1]-qnorm(0.975)*std.hat.dab[r,1])&&theta[1]<=(theta.hat.dab[r,1]+qnorm(0.975)*std.hat.dab[r,1]))
          cover.dab[r,2] = as.numeric(theta[2]>=(theta.hat.dab[r,2]-qnorm(0.975)*std.hat.dab[r,2])&&theta[2]<=(theta.hat.dab[r,2]+qnorm(0.975)*std.hat.dab[r,2]))
          stat.dab[r,] = theta.hat.dab[r,]/std.hat.dab[r,]
          print(paste(N,"-",T,":",r,"/",rep))
          r = r + 1
        }
      }
    }
    res = c(N = N, T = T, rmse.ablasso = sqrt(colMeans((theta.hat - matrix(theta, rep, 2, byrow = TRUE))^2)),
            rmse.ablasso2 = sqrt(colMeans((theta.hat2 - matrix(theta, rep, 2, byrow = TRUE))^2)),
            rmse.ablasso5 = sqrt(colMeans((theta.hat5 - matrix(theta, rep, 2, byrow = TRUE))^2)),
            rmse.ab = sqrt(colMeans((theta.hat.ab - matrix(theta, rep, 2, byrow = TRUE))^2)),
            rmse.dab = sqrt(colMeans((theta.hat.dab - matrix(theta, rep, 2, byrow = TRUE))^2)),
            std.ablasso = colSds(theta.hat), std.ablasso2 = colSds(theta.hat2), std.ablasso5 = colSds(theta.hat5), 
            std.ab = colSds(theta.hat.ab), std.dab = colSds(theta.hat.dab),
            bias.ablasso = colMeans(theta.hat) - theta, bias.ablasso2 = colMeans(theta.hat2) - theta, bias.ablasso5 = colMeans(theta.hat5) - theta, 
            bias.ab = colMeans(theta.hat.ab) - theta, bias.dab = colMeans(theta.hat.dab) - theta,
            ct.ablasso = mean(c1), ct.ablasso2 = mean(c2), ct.ablasso5 = mean(c3), ct.ab = mean(c4), ct.dab = mean(c5),
            length = colMeans(2*std.hat*qnorm(0.975)/sqrt(N*(T-2))), length2 = colMeans(2*std.hat2*qnorm(0.975)/sqrt(N*(T-2))), length5 = colMeans(2*std.hat5*qnorm(0.975)/sqrt(N*(T-2))), 
            length.ab = colMeans(2*std.hat.ab*qnorm(0.975)), length.dab = colMeans(2*std.hat.dab*qnorm(0.975)),
            coverage = colMeans(cover), coverage2 = colMeans(cover2), coverage5 = colMeans(cover5),
            coverage.ab = colMeans(cover.ab), coverage.dab = colMeans(cover.dab)
    )
    
    results = rbind(results, res)
    
    output[[ii]] = list(N = N, T = T, theta.hat = theta.hat, theta.hat2 = theta.hat2, theta.hat5 = theta.hat5, theta.hat.ab = theta.hat.ab, theta.hat.dab = theta.hat.dab,
                        std.hat = std.hat, std.hat2 = std.hat2, std.hat5 = std.hat5, std.hat.ab = std.hat.ab, std.hat.dab = std.hat.dab,
                        stat = stat, stat2 = stat2, stat5 = stat5, stat.ab = stat.ab, stat.dab = stat.dab)
    ii = ii + 1
  }
}

results.all = list(results = results, output = output)
save(results.all, file = "ablasso_simulation.dat")

load("ablasso_simulation.dat")
results = results.all$results
t(round(results[,c("N","T","rmse.ablasso51","rmse.ablasso52","rmse.ab1","rmse.ab2",
                   "std.ablasso51","std.ablasso52","std.ab1","std.ab2",
                   "bias.ablasso51","bias.ablasso52","bias.ab1","bias.ab2",
                   "length51","length52","length.ab1","length.ab2",
                   "coverage51","coverage52","coverage.ab1","coverage.ab2",
                   "ct.ablasso5","ct.ab")],2))
t(round(results[,c("N","T","rmse.ablasso1","rmse.ablasso2","rmse.ab1","rmse.ab2",
                   "std.ablasso1","std.ablasso2","std.ab1","std.ab2",
                   "bias.ablasso1","bias.ablasso2","bias.ab1","bias.ab2",
                   "length1","length2","length.ab1","length.ab2",
                   "coverage1","coverage2","coverage.ab1","coverage.ab2",
                   "ct.ablasso","ct.ab")],2))
t(round(results[,c("N","T","rmse.ablasso21","rmse.ablasso22","rmse.dab1","rmse.dab2",
                   "std.ablasso21","std.ablasso22","std.dab1","std.dab2",
                   "bias.ablasso21","bias.ablasso22","bias.dab1","bias.dab2",
                   "length21","length22","length.dab1","length.dab2",
                   "coverage21","coverage22","coverage.dab1","coverage.dab2",
                   "ct.ablasso2","ct.dab")],2))
