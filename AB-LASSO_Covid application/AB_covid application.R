setwd("~/Documents/GitHub/Arellano-Bond-LASSO/AB-LASSO_Covid application")

library(plm)

data = load("data_weekly_balanced.Rdata")
attach(sdf_week)

############### AB ################
l = 4
dataset.ab = cbind(as.factor(fips), as.factor(week), as.vector(dlogdc), as.vector(school), as.vector(logdc), as.vector(pmask), as.vector(pgather50), as.vector(college), as.vector(pshelter), as.vector(dlogtests))
colnames(dataset.ab) = c("fips", "week", "dlogdc", "school", "logdc", "pmask", "pgather50", "college", "pshelter", "dlogtests")

form.ab = logdc ~ lag(school, 1) + lag(logdc, 1:l) + dlogtests + lag(pmask, 1) + lag(college, 1) + lag(pshelter, 1) + lag(pgather50, 1) | lag(logdc, 2:99) + lag(school, 1:99) + lag(dlogtests, 1:99) + lag(pmask, 1:99) + lag(college, 1:99) + lag(pgather50, 1:99) + lag(pshelter, 1:99) 
fit.ab = pgmm(form.ab, dataset.ab, model = "twosteps", effect = "twoways")
theta.hat.ab = fit.ab$coefficients[[2]][1:(l+6)]  # 6 is the number of treatment variables

HCV.coefs = vcovHC(fit.ab, cluster = 'group')[1:(l+6),1:(l+6)]
std.hat.ab = sqrt(diag(HCV.coefs))
stat.ab = theta.hat.ab/std.hat.ab

results = list(theta.hat.ab = theta.hat.ab, std.hat.ab = std.hat.ab, stat.ab = stat.ab, HCV.coefs = HCV.coefs)
save(results, file = "ab_covid_fat.dat")

load("ab_covid_fat.dat")
# short-run effects
round(results$theta.hat.ab,2)
round(results$std.hat.ab,2)
round(results$stat.ab,2)

# long-run effects of D
HCV.coefs = results$HCV.coefs
lr.ab     = results$theta.hat.ab[1]/(1-sum(results$theta.hat.ab[2:(l+1)]))
jac.lr    = c(1,rep(lr.ab,l))/(1-sum(results$theta.hat.ab[2:(l+1)]))
cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs[c(1,2:(l+1)),c(1,2:(l+1))] %*% jac.lr)

round(lr.ab, 2)
round(cse.lr.ab, 2)
round(lr.ab/cse.lr.ab, 2) # T-stat

# long-run effects of others (C1)
lr.ab     = results$theta.hat.ab[8]/(1-sum(results$theta.hat.ab[2:(l+1)]))
jac.lr    = c(1,rep(lr.ab,l))/(1-sum(results$theta.hat.ab[2:(l+1)]))
cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs[c(8,2:(l+1)),c(8,2:(l+1))] %*% jac.lr)

round(lr.ab, 2)
round(cse.lr.ab, 2)
round(lr.ab/cse.lr.ab, 2) # T-stat

lr.ab     = results$theta.hat.ab[7]/(1-sum(results$theta.hat.ab[2:(l+1)]))
jac.lr    = c(1,rep(lr.ab,l))/(1-sum(results$theta.hat.ab[2:(l+1)]))
cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs[c(7,2:(l+1)),c(7,2:(l+1))] %*% jac.lr)

round(lr.ab, 2)
round(cse.lr.ab, 2)
round(lr.ab/cse.lr.ab, 2) # T-stat

lr.ab     = results$theta.hat.ab[9]/(1-sum(results$theta.hat.ab[2:(l+1)]))
jac.lr    = c(1,rep(lr.ab,l))/(1-sum(results$theta.hat.ab[2:(l+1)]))
cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs[c(9,2:(l+1)),c(9,2:(l+1))] %*% jac.lr)

round(lr.ab, 2)
round(cse.lr.ab, 2)
round(lr.ab/cse.lr.ab, 2) # T-stat

lr.ab     = results$theta.hat.ab[10]/(1-sum(results$theta.hat.ab[2:(l+1)]))
jac.lr    = c(1,rep(lr.ab,l))/(1-sum(results$theta.hat.ab[2:(l+1)]))
cse.lr.ab = sqrt(t(jac.lr) %*% HCV.coefs[c(10,2:(l+1)),c(10,2:(l+1))] %*% jac.lr)

round(lr.ab, 2)
round(cse.lr.ab, 2)
round(lr.ab/cse.lr.ab, 2) # T-stat
