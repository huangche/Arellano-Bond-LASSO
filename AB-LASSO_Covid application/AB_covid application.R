setwd("~/Downloads/Covid-19")
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

############### AB ################
l = 4
dataset.ab = cbind(as.factor(fips), as.factor(week), as.vector(dlogdc), as.vector(school), as.vector(logdc), as.vector(pmask), as.vector(pgather50), as.vector(college), as.vector(pshelter), as.vector(dlogtests))
colnames(dataset.ab) = c("fips", "week", "dlogdc", "school", "logdc", "pmask", "pgather50", "college", "pshelter", "dlogtests")

form.ab = logdc ~ lag(school, 1) + lag(logdc, 1:l) + dlogtests + lag(pmask, 1) + lag(college, 1) + lag(pshelter, 1) + lag(pgather50, 1) | lag(logdc, 2:99) + lag(school, 1:99) + lag(dlogtests, 1:99) + lag(pmask, 1:99) + lag(college, 1:99) + lag(pgather50, 1:99) + lag(pshelter, 1:99) 
fit.ab = pgmm(form.ab, dataset.ab, model = "twosteps", effect = "twoways")
theta.hat.ab = fit.ab$coefficients[[2]][1:(l+6)]  # 6 is the number of covariables (excluding lags of Y)
std.hat.ab = summary(fit.ab)$coefficients[1:(l+6),"Std. Error"]
stat.ab = theta.hat.ab/std.hat.ab

coefs.ab  = coef(fit.ab)
HCV.coefs = vcovHC(fit.ab, cluster = 'group')
cse.ab    = sqrt(diag(HCV.coefs)) # Clustered std errors

lr.ab     = coefs.ab[1]/(1-sum(coefs.ab[2:(l+1)]))
jac.lr    = c(1,rep(lr.ab,l))/(1-sum(coefs.ab[2:(l+1)]))
cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs[1:(l+1),1:(l+1)] %*% jac.lr)

results = list(fit.ab = fit.ab, theta.hat.ab = theta.hat.ab, std.hat.ab = std.hat.ab, stat.ab = stat.ab, lr.ab = lr.ab, cse.lr.ab = cse.lr.ab)
save(results, file = "ab_covid.dat")

load("ab_covid.dat")
# short-run effects
round(results$theta.hat.ab,2)
round(results$std.hat.ab,2)
round(results$stat.ab,2)

# long-run effects of D
round(results$lr.ab, 2)
round(results$cse.lr.ab, 2)
round(results$lr.ab/results$cse.lr.ab, 2) # T-stat

fit.ab = results$fit.ab
coefs.ab  = coef(fit.ab)
HCV.coefs = vcovHC(fit.ab, cluster = 'group')
cse.ab    = sqrt(diag(HCV.coefs)) # Clustered std errors

# long-run effects of others (C1)
lr.ab     = coefs.ab[7]/(1-sum(coefs.ab[2:(l+1)]))
jac.lr    = c(1,rep(lr.ab,l))/(1-sum(coefs.ab[2:(l+1)]))
cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs[c(7,2:(l+1)),c(7,2:(l+1))] %*% jac.lr)

round(lr.ab, 2)
round(cse.lr.ab, 2)
round(lr.ab/cse.lr.ab, 2) # T-stat

lr.ab     = coefs.ab[8]/(1-sum(coefs.ab[2:(l+1)]))
jac.lr    = c(1,rep(lr.ab,l))/(1-sum(coefs.ab[2:(l+1)]))
cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs[c(8,2:(l+1)),c(8,2:(l+1))] %*% jac.lr)

round(lr.ab, 2)
round(cse.lr.ab, 2)
round(lr.ab/cse.lr.ab, 2) # T-stat

lr.ab     = coefs.ab[9]/(1-sum(coefs.ab[2:(l+1)]))
jac.lr    = c(1,rep(lr.ab,l))/(1-sum(coefs.ab[2:(l+1)]))
cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs[c(9,2:(l+1)),c(9,2:(l+1))] %*% jac.lr)

round(lr.ab, 2)
round(cse.lr.ab, 2)
round(lr.ab/cse.lr.ab, 2) # T-stat

lr.ab     = coefs.ab[10]/(1-sum(coefs.ab[2:(l+1)]))
jac.lr    = c(1,rep(lr.ab,l))/(1-sum(coefs.ab[2:(l+1)]))
cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs[c(10,2:(l+1)),c(10,2:(l+1))] %*% jac.lr)

round(lr.ab, 2)
round(cse.lr.ab, 2)
round(lr.ab/cse.lr.ab, 2) # T-stat