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

############### DAB ################
l = 4
dataset.ab = cbind(as.factor(fips), as.factor(week), as.vector(dlogdc), as.vector(school), as.vector(logdc), as.vector(pmask), as.vector(pgather50), as.vector(college), as.vector(pshelter), as.vector(dlogtests))
colnames(dataset.ab) = c("fips", "week", "dlogdc", "school", "logdc", "pmask", "pgather50", "college", "pshelter", "dlogtests")

form.ab = logdc ~ lag(school, 1) + lag(logdc, 1:l) + dlogtests + lag(pmask, 1) + lag(college, 1) + lag(pshelter, 1) + lag(pgather50, 1) | lag(logdc, 2:99) + lag(school, 1:99) + lag(dlogtests, 1:99) + lag(pmask, 1:99) + lag(college, 1:99) + lag(pgather50, 1:99) + lag(pshelter, 1:99) 
sample1 = sample(unique(dataset.ab[,"fips"]), ceiling(N/2), replace = FALSE)
fit.dab1 = pgmm(form.ab, dataset.ab[which(dataset.ab[,"fips"] %in% sample1), ], model = "twosteps", effect = "twoways")
fit.dab2 = pgmm(form.ab, dataset.ab[-which(dataset.ab[,"fips"] %in% sample1), ], model = "twosteps", effect = "twoways")
fit.ab = pgmm(form.ab, dataset.ab, model = "twosteps", effect = "twoways")

theta.hat.dab = 2*fit.ab$coefficients[[2]][1:(l+6)] - (fit.dab1$coefficients[[2]][1:(l+6)] + fit.dab2$coefficients[[2]][1:(l+6)])/2  # 6 is the number of treatment variables
std.hat.dab = summary(fit.ab)$coefficients[1:(l+6),"Std. Error"]
stat.dab = theta.hat.dab/std.hat.dab

coefs.dab  = 2*coef(fit.ab) - (coef(fit.dab1) + coef(fit.dab2))/2
HCV.coefs = vcovHC(fit.ab, cluster = 'group')
cse.dab    = sqrt(diag(HCV.coefs)) # Clustered std errors
lr.dab     = coefs.dab[1]/(1-sum(coefs.dab[2:(l+1)]))
jac.lr    = c(1,rep(lr.dab,l))/(1-sum(coefs.dab[2:(l+1)]))
cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs[1:(l+1),1:(l+1)] %*% jac.lr)

results = list(fit.dab1 = fit.dab1, fit.dab2 = fit.dab2, fit.ab = fit.ab, theta.hat.dab = theta.hat.dab, std.hat.dab = std.hat.dab, stat.dab = stat.dab, lr.dab = lr.dab, cse.lr.dab = cse.lr.dab)
save(results, file = "dab_covid.dat")

load("dab_covid.dat")
# short-run effects
round(results$theta.hat.dab,2)
round(results$std.hat.dab,2)
round(results$stat.dab,2)

# long-run effects of D
round(results$lr.dab, 2)
round(results$cse.lr.dab, 2)
round(results$lr.dab/results$cse.lr.dab, 2) # T-stat

fit.ab = results$fit.ab
fit.dab1 = results$fit.dab1
fit.dab2 = results$fit.dab2

coefs.dab  = 2*coef(fit.ab) - (coef(fit.dab1) + coef(fit.dab2))/2
HCV.coefs = vcovHC(fit.ab, cluster = 'group')
cse.dab    = sqrt(diag(HCV.coefs)) # Clustered std errors

# long-run effects of others (C1)
lr.dab     = coefs.dab[7]/(1-sum(coefs.dab[2:(l+1)]))
jac.lr    = c(1,rep(lr.dab,l))/(1-sum(coefs.dab[2:(l+1)]))
cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs[c(7,2:(l+1)),c(7,2:(l+1))] %*% jac.lr)

round(lr.dab, 2)
round(cse.lr.dab, 2)
round(lr.dab/cse.lr.dab, 2) # T-stat

lr.dab     = coefs.dab[8]/(1-sum(coefs.dab[2:(l+1)]))
jac.lr    = c(1,rep(lr.dab,l))/(1-sum(coefs.dab[2:(l+1)]))
cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs[c(8,2:(l+1)),c(8,2:(l+1))] %*% jac.lr)

round(lr.dab, 2)
round(cse.lr.dab, 2)
round(lr.dab/cse.lr.dab, 2) # T-stat

lr.dab     = coefs.dab[9]/(1-sum(coefs.dab[2:(l+1)]))
jac.lr    = c(1,rep(lr.dab,l))/(1-sum(coefs.dab[2:(l+1)]))
cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs[c(9,2:(l+1)),c(9,2:(l+1))] %*% jac.lr)

round(lr.dab, 2)
round(cse.lr.dab, 2)
round(lr.dab/cse.lr.dab, 2) # T-stat

lr.dab     = coefs.dab[10]/(1-sum(coefs.dab[2:(l+1)]))
jac.lr    = c(1,rep(lr.dab,l))/(1-sum(coefs.dab[2:(l+1)]))
cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs[c(10,2:(l+1)),c(10,2:(l+1))] %*% jac.lr)

round(lr.dab, 2)
round(cse.lr.dab, 2)
round(lr.dab/cse.lr.dab, 2) # T-stat