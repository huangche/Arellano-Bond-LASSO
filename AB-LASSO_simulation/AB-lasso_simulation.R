####### for a single replication
setwd("~/Documents/GitHub/Arellano-Bond-LASSO/AB-LASSO_simulation")

libraries = c("mvtnorm", "plm", "hdm", "glmnet", "matrixStats", "dpm")
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

NN = c(100,200)
TT = c(20,30,40,50,60)
sigmas = c(sqrt(2.96),sqrt(6.58),1)
theta = c(0.75,0.25) 
theta.lr = theta[2]/(1-theta[1])
rho = 0.5
phi = -0.17
pi = 0.67
df = 4
nboot = 100
S2 = 1

sigma_alpha = sigmas[1]
sigma_eps.d = sigmas[2]
sigma_eps.y = sigmas[3]

for(nn in 1:length(NN)){
  for(tt in 1:length(TT)){
    N = NN[nn]
    T = TT[tt]
    con = 1.1
    r = 1
    while(r<2){
      alpha = rnorm(N, 0, sigma_alpha)
      eps.y = eps.d = Y = D = matrix(0, nrow = T+1, ncol = N)
      for(t in 1:(T+1)){
        eps.y[t,] = rt(N, df) 
        eps.d[t,] = rt(N, df) 
      }
      #eps.y[which(eps.d>0)] = 1.5*eps.y[which(eps.d>0)]  #for heteroskedastic case
      D[1,] = ((phi+pi*(1-theta[1]))/((1-rho)*(1-theta[1])-theta[2]*phi))*alpha + eps.d[1,]
      Y[1,] = ((theta[2]*pi+(1-rho))/((1-rho)*(1-theta[1])-theta[2]*phi))*alpha + eps.y[1,]
      for(t in 2:(T+1)){
        D[t,] = rho*D[t-1,] + phi*Y[t-1,] + pi*alpha + eps.d[t,]
        Y[t,] = theta[1]*Y[t-1,] + theta[2]*D[t,] + alpha + eps.y[t,]
      }
      Y = Y[-1,]
      D = D[-1,]
      
      ############### AB-LASSO ################
      ## without sample splitting
      cat("Starting AB-LASSO...\n")
      t1 = Sys.time()
      Y.t = D.t = matrix(0, T-1, N)
      for(i in 1:N){
        for(t in 1:(T-1)){
          Y.t[t,i] = (Y[t,i] - mean(Y[(t+1):T,i]))*sqrt((T-t)/(T-t+1))
          D.t[t,i] = (D[t,i] - mean(D[(t+1):T,i]))*sqrt((T-t)/(T-t+1))
        }
      } 
      W1 = Y.t[1:(nrow(Y.t)-1),]
      W2 = D.t[2:nrow(D.t),]
      Z = list()
      for(t in 2:(T-1)){
        y1 = W1[t-1,]
        y2 = W2[t-1,]
        if(t>2){x = cbind(t(Y[1:(t-1),]), t(D[1:(t),]))}else{x = cbind(Y[1:(t-1),], t(D[1:(t),]))}
        fit1 = rlasso(x, y1, penalty = list(homoscedastic = "none", lambda.start = con*sqrt(N)*qnorm(1-0.1/(2*dim(x)[2]))))
        fit2 = rlasso(x, y2, penalty = list(homoscedastic = "none", lambda.start = con*sqrt(N)*qnorm(1-0.1/(2*dim(x)[2]))))
        Z[[t-1]] = cbind(predict(fit1), predict(fit2))
      }
      sum1 = matrix(0, 2, 2)
      sum2 = rep(0, 2)
      for(i in 1:N){
        for(t in 2:(T-1)){
          sum1 = sum1 + Z[[t-1]][i,]%*%t(c(W1[t-1,i],W2[t-1,i]))
          sum2 = sum2 + Z[[t-1]][i,]*Y.t[t,i]
        }
      }
      theta.hat = solve(sum1)%*%sum2
      vps.t = matrix(as.vector(Y.t[2:nrow(Y.t),]) - cbind(as.vector(W1),as.vector(W2))%*%theta.hat, nrow = T-2, ncol = N)
      mu = matrix(0, N, 2)
      for(i in 1:N){
        for(t in 2:(T-1)){
          mu[i,] = mu[i,] + Z[[t-1]][i,]*vps.t[t-1,i]
        }
        mu[i,] = mu[i,]/(T-2)
      }
      sum3 = sum4 = matrix(0, 2, 2)
      for(i in 1:N){
        for(t in 2:(T-1)){
          sum3 = sum3 + (Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])%*%t(Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])
        }
      }
      Sigma = sum3 
      se = sqrt(diag(solve(sum1)%*%Sigma%*%t(solve(sum1))))
      vcv = solve(sum1)%*%Sigma%*%t(solve(sum1))
      cover = as.numeric(theta>=(theta.hat-qnorm(0.975)*se)&theta<=(theta.hat+qnorm(0.975)*se))
      stat = theta.hat/se
      #long-run effect
      lr = theta.hat[2]/(1-theta.hat[1])
      jac.lr = c(lr,1)/(1-theta.hat[1])
      cse.lr = sqrt(t(jac.lr) %*% vcv %*% jac.lr)
      cover.lr = as.numeric(theta.lr>=(lr-qnorm(0.975)*cse.lr)&theta.lr<=(lr+qnorm(0.975)*cse.lr))
      stat.lr = lr/cse.lr
      c1 = difftime(Sys.time(), t1, units = "secs")
      
      ## 2 folds random sample splitting
      cat("Starting 2 folds SS...\n")
      t2 = Sys.time()
      Kf = 2
      theta.hat.ss = se.ss = numeric(0)
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
        W1 = Y.t[1:(nrow(Y.t)-1),]
        W2 = D.t[2:nrow(D.t),]
        
        Z = list()
        for(t in 2:(T-1)){
          Z[[t-1]] = matrix(0, N, 2)
          y1 = y2 = rep(0, N)
          x = matrix(0, N, ((t-1)+(t)))
          for(b in 1:length(I)){
            y1[-I[[b]]] = W1[t-1,-I[[b]]]  #s2 - auxiliary sample
            y2[-I[[b]]] = W2[t-1,-I[[b]]]
            if(t>2){x[I[[b]],] = cbind(t(Y[1:(t-1),I[[b]]]), t(D[1:(t),I[[b]]]))}else{x[I[[b]],]  = cbind(Y[1:(t-1),I[[b]]], t(D[1:(t),I[[b]]]))}
            if(t>2){x[-I[[b]],] = cbind(t(Y[1:(t-1),-I[[b]]]), t(D[1:(t),-I[[b]]]))}else{x[-I[[b]],] = cbind(Y[1:(t-1),-I[[b]]], t(D[1:(t),-I[[b]]]))}
            fit1_s2 = rlasso(x[-I[[b]],], y1[-I[[b]]], penalty = list(homoscedastic = "none", lambda.start = con*sqrt(N-length(I[[b]]))*qnorm(1-0.1/(2*dim(x)[2]))))
            fit2_s2 = rlasso(x[-I[[b]],], y2[-I[[b]]], penalty = list(homoscedastic = "none", lambda.start = con*sqrt(N-length(I[[b]]))*qnorm(1-0.1/(2*dim(x)[2]))))
            Z[[t-1]][I[[b]],] = cbind(predict(fit1_s2, x[I[[b]],]), predict(fit2_s2, x[I[[b]],])) 
          }
        }
        theta.est = matrix(0, 2, length(I))
        for(b in 1:length(I)){
          sum1 = matrix(0, 2, 2)
          sum2 = rep(0, 2)
          for(i in I[[b]]){
            for(t in 2:(T-1)){
              sum1 = sum1 + Z[[t-1]][i,]%*%t(c(W1[t-1,i],W2[t-1,i]))
              sum2 = sum2 + Z[[t-1]][i,]*Y.t[t,i]
            }
          }
          theta.est[,b] = solve(sum1)%*%sum2
        }
        theta.hat.ss = rbind(theta.hat.ss, rowMeans(theta.est))
        vps.t = matrix(as.vector(Y.t[2:nrow(Y.t),]) - cbind(as.vector(W1),as.vector(W2))%*%rowMeans(theta.est), nrow = T-2, ncol = N)
        mu = matrix(0, N, 2)
        for(i in 1:N){
          for(t in 2:(T-1)){
            mu[i,] = mu[i,] + Z[[t-1]][i,]*vps.t[t-1,i]
          }
          mu[i,] = mu[i,]/(T-2)
        }
        sum3 = sum4 = matrix(0, 2, 2)
        for(i in 1:N){
          for(t in 2:(T-1)){
            sum3 = sum3 + (Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])%*%t(Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])
          }
        }
        Sigma = sum3 
        sum1 = matrix(0, 2, 2)
        for(i in 1:N){
          for(t in 2:(T-1)){
            sum1 = sum1 + Z[[t-1]][i,]%*%t(c(W1[t-1,i],W2[t-1,i]))
          }
        }
        se.ss = rbind(se.ss, sqrt(diag(solve(sum1)%*%Sigma%*%t(solve(sum1)))))
        vcv.all[[ib]] = solve(sum1)%*%Sigma%*%t(solve(sum1))
      }
      theta.hat2 = colMedians(theta.hat.ss)
      se2 = colMedians(se.ss)
      sd2 = colSds(theta.hat.ss)
      cover2 = as.numeric(theta>=(theta.hat2-qnorm(0.975)*se2)&theta<=(theta.hat2+qnorm(0.975)*se2))
      stat2 = theta.hat2/se2
      #long-run effect
      lr.all = cse.lr.all = rep(0, nboot)
      for(ib in 1:nboot){
        lr.all[ib] = theta.hat.ss[ib,2]/(1-theta.hat.ss[ib,1])
        jac.lr = c(lr.all[ib],1)/(1-theta.hat.ss[ib,1])
        cse.lr.all[ib] = sqrt(t(jac.lr) %*% vcv.all[[ib]] %*% jac.lr)  
      }
      lr2 = median(lr.all)
      cse.lr2 = median(cse.lr.all)
      cover.lr2 = as.numeric(theta.lr>=(lr2-qnorm(0.975)*cse.lr2)&theta.lr<=(lr2+qnorm(0.975)*cse.lr2))
      stat.lr2 = lr2/cse.lr2
      c2 = difftime(Sys.time(), t2, units = "secs")
      
      ## 5 folds random sample splitting
      cat("Starting 5 folds SS...\n")
      t3 = Sys.time()
      Kf = 5
      theta.hat.ss = se.ss = numeric(0)
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
        W1 = Y.t[1:(nrow(Y.t)-1),]
        W2 = D.t[2:nrow(D.t),]
        
        Z = list()
        for(t in 2:(T-1)){
          Z[[t-1]] = matrix(0, N, 2)
          y1 = y2 = rep(0, N)
          x = matrix(0, N, ((t-1)+(t)))
          for(b in 1:length(I)){
            y1[-I[[b]]] = W1[t-1,-I[[b]]]  #s2 - auxiliary sample
            y2[-I[[b]]] = W2[t-1,-I[[b]]]
            if(t>2){x[I[[b]],] = cbind(t(Y[1:(t-1),I[[b]]]), t(D[1:(t),I[[b]]]))}else{x[I[[b]],]  = cbind(Y[1:(t-1),I[[b]]], t(D[1:(t),I[[b]]]))}
            if(t>2){x[-I[[b]],] = cbind(t(Y[1:(t-1),-I[[b]]]), t(D[1:(t),-I[[b]]]))}else{x[-I[[b]],] = cbind(Y[1:(t-1),-I[[b]]], t(D[1:(t),-I[[b]]]))}
            fit1_s2 = rlasso(x[-I[[b]],], y1[-I[[b]]], penalty = list(homoscedastic = "none", lambda.start = con*sqrt(N-length(I[[b]]))*qnorm(1-0.1/(2*dim(x)[2]))))
            fit2_s2 = rlasso(x[-I[[b]],], y2[-I[[b]]], penalty = list(homoscedastic = "none", lambda.start = con*sqrt(N-length(I[[b]]))*qnorm(1-0.1/(2*dim(x)[2]))))
            Z[[t-1]][I[[b]],] = cbind(predict(fit1_s2, x[I[[b]],]), predict(fit2_s2, x[I[[b]],])) 
          }
        }
        theta.est = matrix(0, 2, length(I))
        for(b in 1:length(I)){
          sum1 = matrix(0, 2, 2)
          sum2 = rep(0, 2)
          for(i in I[[b]]){
            for(t in 2:(T-1)){
              sum1 = sum1 + Z[[t-1]][i,]%*%t(c(W1[t-1,i],W2[t-1,i]))
              sum2 = sum2 + Z[[t-1]][i,]*Y.t[t,i]
            }
          }
          theta.est[,b] = solve(sum1)%*%sum2
        }
        theta.hat.ss = rbind(theta.hat.ss, rowMeans(theta.est))
        vps.t = matrix(as.vector(Y.t[2:nrow(Y.t),]) - cbind(as.vector(W1),as.vector(W2))%*%rowMeans(theta.est), nrow = T-2, ncol = N)
        mu = matrix(0, N, 2)
        for(i in 1:N){
          for(t in 2:(T-1)){
            mu[i,] = mu[i,] + Z[[t-1]][i,]*vps.t[t-1,i]
          }
          mu[i,] = mu[i,]/(T-2)
        }
        sum3 = sum4 = matrix(0, 2, 2)
        for(i in 1:N){
          for(t in 2:(T-1)){
            sum3 = sum3 + (Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])%*%t(Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])
          }
        }
        Sigma = sum3 
        sum1 = matrix(0, 2, 2)
        for(i in 1:N){
          for(t in 2:(T-1)){
            sum1 = sum1 + Z[[t-1]][i,]%*%t(c(W1[t-1,i],W2[t-1,i]))
          }
        }
        se.ss = rbind(se.ss, sqrt(diag(solve(sum1)%*%Sigma%*%t(solve(sum1)))))
        vcv.all[[ib]] = solve(sum1)%*%Sigma%*%t(solve(sum1))
      }
      theta.hat5 = colMedians(theta.hat.ss)
      se5 = colMedians(se.ss)
      sd5 = colSds(theta.hat.ss)
      cover5 = as.numeric(theta>=(theta.hat5-qnorm(0.975)*se5)&theta<=(theta.hat5+qnorm(0.975)*se5))
      stat5 = theta.hat5/se5
      #long-run effect
      lr.all = cse.lr.all = rep(0, nboot)
      for(ib in 1:nboot){
        lr.all[ib] = theta.hat.ss[ib,2]/(1-theta.hat.ss[ib,1])
        jac.lr = c(lr.all[ib],1)/(1-theta.hat.ss[ib,1])
        cse.lr.all[ib] = sqrt(t(jac.lr) %*% vcv.all[[ib]] %*% jac.lr)  
      }
      lr5 = median(lr.all)
      cse.lr5 = median(cse.lr.all)
      cover.lr5 = as.numeric(theta.lr>=(lr5-qnorm(0.975)*cse.lr5)&theta.lr<=(lr5+qnorm(0.975)*cse.lr5))
      stat.lr5 = lr5/cse.lr5
      c3 = difftime(Sys.time(), t3, units = "secs")
      
      ############### AB and others ################
      index.i = index.t = numeric(0)
      for(i in 1:N){
        index.i = c(index.i, rep(i,T))
        index.t = c(index.t, c(2000:(2000+T-1)))
      }
      dataset.ab = cbind(as.factor(index.i), as.factor(index.t), as.vector(Y), as.vector(D))
      colnames(dataset.ab) = c("firm", "year", "y.ab", "d.ab")
      ## DFE-A
      t4 = Sys.time()
      data.fe = pdata.frame(dataset.ab, index = c("firm","year"))
      data.fe$y.ab.lag1 = lag(data.fe$y.ab, 1)
      form.feabc = y.ab ~ y.ab.lag1 + d.ab + firm
      fit.feabc = lm(form.feabc, data.fe, x = TRUE, na.action = na.omit)
      res.feabc = fit.feabc$residuals
      jac = solve(t(fit.feabc$x)%*%fit.feabc$x/length(res.feabc))[2:3,2:3]
      indexes = c(1:length(res.feabc))
      indexes = indexes[-c(1+c(0:(N-1))*length(res.feabc)/N)]
      bscore = t(fit.feabc$x[indexes, 2:3])%*%res.feabc[indexes-1]/length(indexes)
      bias = -jac%*%bscore*N/length(res.feabc)
      theta.hat.feabc = theta.hat.fe - bias
      se.feabc = se.fe
      cover.feabc = as.numeric(theta>=(theta.hat.feabc-qnorm(0.975)*se.feabc)&theta<=(theta.hat.feabc+qnorm(0.975)*se.feabc))
      stat.feabc = theta.hat.feabc/se.feabc
      #long-run effect
      lr.feabc = theta.hat.feabc[2]/(1-theta.hat.feabc[1])
      jac.lr = c(lr.feabc,1)/(1-theta.hat.feabc[1])
      cse.lr.feabc = sqrt(t(jac.lr) %*% HCV.coefs %*% jac.lr)
      cover.lr.feabc = as.numeric(theta.lr>=(lr.feabc-qnorm(0.975)*cse.lr.feabc)&theta.lr<=(lr.feabc+qnorm(0.975)*cse.lr.feabc))
      stat.lr.feabc = lr.feabc/cse.lr.feabc
      c4 = difftime(Sys.time(), t4, units = "secs")
      
      ## Likelihood method
      cat("Starting ML...\n")
      if(N==100&T==60){
        theta.hat.ml = se.ml = cover.ml = stat.ml = rep(NA, 2)
        lr.ml = cse.lr.ml = cover.lr.ml = stat.lr.ml = NA
        c5 = 0
      }else{
        if(N==100&T==50){
          theta.hat.ml = se.ml = cover.ml = stat.ml = rep(NA, 2)
          lr.ml = cse.lr.ml = cover.lr.ml = stat.lr.ml = NA
          c5 = 0
        }else{
          t5 = Sys.time()
          dataset.ml = panel_data(data.frame(dataset.ab), id = firm, wave = year)
          fit.ml = dpm(y.ab ~ pre(d.ab), data = dataset.ml, error.inv = FALSE, information = "observed")
          try.ml = try(summary(fit.ml), silent = TRUE)
          if(class(try.ml)[1]!="try-error"){
            theta.hat.ml = rev(try.ml$coefficients[,"Est."])
            se.ml = rev(try.ml$coefficients[,"S.E."])
          }else{
            par_table = as.data.frame(fit.ml@ParTable)
            coef_table = subset(par_table, op == "~")
            idx = c(1:nrow(coef_table))
            theta.hat.ml = c(mean(coef_table[idx[idx %% 2 == 0], "est"]), mean(coef_table[idx[idx %% 2 != 0], "est"]))
            se.ml = c(mean(coef_table[idx[idx %% 2 == 0], "se"]), mean(coef_table[idx[idx %% 2 != 0], "se"]))
          }
          cover.ml = as.numeric(theta>=(theta.hat.ml-qnorm(0.975)*se.ml)&theta<=(theta.hat.ml+qnorm(0.975)*se.ml))
          stat.ml = theta.hat.ml/se.ml
          #long-run effect
          lr.ml = theta.hat.ml[2]/(1-theta.hat.ml[1])
          jac.lr = c(lr.ml,1)/(1-theta.hat.ml[1])
          HCV.coefs = fit.ml@vcov[[3]][1:2,1:2]
          if(is.matrix(HCV.coefs)){
            diag(HCV.coefs) = rev(diag(fit.ml@vcov[[3]][1:2,1:2]))
            cse.lr.ml = sqrt(t(jac.lr) %*% HCV.coefs %*% jac.lr)
            cover.lr.ml = as.numeric(theta.lr>=(lr.ml-qnorm(0.975)*cse.lr.ml)&theta.lr<=(lr.ml+qnorm(0.975)*cse.lr.ml))
            stat.lr.ml = lr.ml/cse.lr.ml
          }else{
            cse.lr.ml = cover.lr.ml = stat.lr.ml = NA
          }
          c5 = difftime(Sys.time(), t5, units = "secs")
        }
      }
      cat("Finishing ML...\n")
      
      ## two-step AB 
      t6 = Sys.time()
      fit.ab = try(pgmm(y.ab ~ lag(y.ab,1) + d.ab | lag(y.ab,2:9999) + lag(d.ab,1:9999), dataset.ab, model = "twosteps", effect = "individual"))
      if(class(fit.ab)[1]!="try-error"){
        try.ab = try(summary(fit.ab), silent = TRUE)
        c6 = difftime(Sys.time(), t6, units = "secs")
      }
      
      ## debiased two-step AB with 1 partition
      if(class(fit.ab)[1]!="try-error"&&class(try.ab)[1]!="try-error"){
        t7 = Sys.time()
        coeff.dab = se.dab = 0*fit.ab$coefficients[[2]][1:2]
        dab.fit1 = dab.fit2 = try.dab1 = try.dab2 = list()
        s = 1
        #set.seed(202304)
        while(s<=S2){
          sample1 = sample(N, ceiling(N/2), replace = FALSE)
          dab.fit1[[s]] = try(pgmm(y.ab ~ lag(y.ab,1) + d.ab | lag(y.ab,2:9999) + lag(d.ab,1:9999), dataset.ab[as.double(index.i) %in% sample1, ], model = "twosteps", effect = "individual"))
          dab.fit2[[s]] = try(pgmm(y.ab ~ lag(y.ab,1) + d.ab | lag(y.ab,2:9999) + lag(d.ab,1:9999), dataset.ab[!as.double(index.i) %in% sample1, ], model = "twosteps", effect = "individual"))
          if(class(dab.fit1[[s]])[1]!="try-error"&&class(dab.fit2[[s]])[1]!="try-error"){
            try.dab1[[s]] = try(summary(dab.fit1[[s]]), silent = TRUE)
            try.dab2[[s]] = try(summary(dab.fit2[[s]]), silent = TRUE)
            if(class(try.dab1[[s]])[1]!="try-error"&&class(try.dab2[[s]])[1]!="try-error"){
              coeff.dab = coeff.dab + ((dab.fit1[[s]]$coefficients[[2]][1:2] + dab.fit2[[s]]$coefficients[[2]][1:2])/2)/S2
              s = s + 1
            }else{s = S2 + 2}
          }else{s = S2 + 2}
        }
        c7 = difftime(Sys.time(), t7, units = "secs")
        if(s==(S2+1)){
          theta.hat.ab = fit.ab$coefficients[[2]][1:2]
          se.ab = try.ab$coefficients[,"Std. Error"]
          #HCV.coefs = vcovHC(fit.ab, cluster = 'group')
          #se.ab = sqrt(diag(HCV.coefs))[1:2] #Clustered std error (get same results for AB, just use the direct output)
          cover.ab = as.numeric(theta>=(theta.hat.ab-qnorm(0.975)*se.ab)&theta<=(theta.hat.ab+qnorm(0.975)*se.ab))
          stat.ab = theta.hat.ab/se.ab
          theta.hat.dab = 2*fit.ab$coefficients[[2]][1:2] - coeff.dab
          se.dab = se.ab
          cover.dab = as.numeric(theta>=(theta.hat.dab-qnorm(0.975)*se.dab)&theta<=(theta.hat.dab+qnorm(0.975)*se.dab))
          stat.dab = theta.hat.dab/se.dab
          #long-run effect
          lr.ab = theta.hat.ab[2]/(1-theta.hat.ab[1])
          jac.lr = c(lr.ab,1)/(1-theta.hat.ab[1])
          HCV.coefs = vcovHC(fit.ab, cluster = 'group')
          cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs %*% jac.lr)
          cover.lr.ab = as.numeric(theta.lr>=(lr.ab-qnorm(0.975)*cse.lr.ab)&theta.lr<=(lr.ab+qnorm(0.975)*cse.lr.ab))
          stat.lr.ab = lr.ab/cse.lr.ab
          lr.dab = theta.hat.dab[2]/(1-theta.hat.dab[1])
          jac.lr = c(lr.dab,1)/(1-theta.hat.dab[1])
          cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs %*% jac.lr)
          cover.lr.dab = as.numeric(theta.lr>=(lr.dab-qnorm(0.975)*cse.lr.dab)&theta.lr<=(lr.dab+qnorm(0.975)*cse.lr.dab))
          stat.lr.dab = lr.dab/cse.lr.dab
          r = r + 1
        }
      }
    }
    time = c(c1,c2,c3,c4,c5,c6,c7)
    result = list(theta.hat = theta.hat, theta.hat2 = theta.hat2, theta.hat5 = theta.hat5, theta.hat.feabc = theta.hat.feabc, theta.hat.ml = theta.hat.ml, 
                  theta.hat.ab = theta.hat.ab, theta.hat.dab = theta.hat.dab,
                  se = se, se2 = se2, se5 = se5, se.feabc = se.feabc, se.ml = se.ml, 
                  se.ab = se.ab, se.dab = se.dab,
                  cover = cover, cover2 = cover2, cover5 = cover5, cover.feabc = cover.feabc, cover.ml = cover.ml, 
                  cover.ab = cover.ab, cover.dab = cover.dab,
                  stat = stat, stat2 = stat2, stat5 = stat5, stat.feabc = stat.feabc, stat.ml = stat.ml, 
                  stat.ab = stat.ab, stat.dab = stat.dab, 
                  lr = lr, cse.lr = cse.lr, cover.lr = cover.lr, lr2 = lr2, cse.lr2 = cse.lr2, cover.lr2 = cover.lr2, lr5 = lr5, cse.lr5 = cse.lr5, cover.lr5 = cover.lr5, 
                  lr.feabc = lr.feabc, cse.lr.feabc = cse.lr.feabc, cover.lr.feabc = cover.lr.feabc,  
                  lr.ml = lr.ml, cse.lr.ml = cse.lr.ml, cover.lr.ml = cover.lr.ml, lr.ab = lr.ab, cse.lr.ab = cse.lr.ab, cover.lr.ab = cover.lr.ab, lr.dab = lr.dab, cse.lr.dab = cse.lr.dab, cover.lr.dab = cover.lr.dab, 
                  stat.lr = stat.lr, stat.lr2 = stat.lr2, stat.lr5 = stat.lr5, stat.lr.feabc = stat.lr.feabc, stat.lr.ml = stat.lr.ml, 
                  stat.lr.ab = stat.lr.ab, stat.lr.dab = stat.lr.dab, 
                  sd2 = sd2, sd5 = sd5, time = time)
    
    out_dir = paste("ablasso_N",N,"_T",T,"QE2",sep="")
    if(!dir.exists(out_dir)){dir.create(out_dir)}
    filename = file.path(out_dir, paste0("result_", rep_id, ".dat"))
    save(result, file = filename)
  }
}

