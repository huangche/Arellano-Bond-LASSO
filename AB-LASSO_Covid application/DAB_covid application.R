setwd("~/Documents/GitHub/Arellano-Bond-LASSO/AB-LASSO_Covid application")

library(plm)

data = load("data_weekly_balanced.Rdata")
attach(sdf_week)

############### DAB ################
l = 4
dataset.ab = cbind(as.factor(fips), as.factor(week), as.vector(dlogdc), as.vector(school), as.vector(logdc), as.vector(pmask), as.vector(pgather50), as.vector(college), as.vector(pshelter), as.vector(dlogtests))
colnames(dataset.ab) = c("fips", "week", "dlogdc", "school", "logdc", "pmask", "pgather50", "college", "pshelter", "dlogtests")
set.seed(202302)
form.ab = logdc ~ lag(school, 1) + lag(logdc, 1:l) + dlogtests + lag(pmask, 1) + lag(college, 1) + lag(pshelter, 1) + lag(pgather50, 1) | lag(logdc, 2:99) + lag(school, 1:99) + lag(dlogtests, 1:99) + lag(pmask, 1:99) + lag(college, 1:99) + lag(pgather50, 1:99) + lag(pshelter, 1:99) 
sample1 = sample(unique(dataset.ab[,"fips"]), ceiling(N/2), replace = FALSE)
fit.dab1 = pgmm(form.ab, dataset.ab[which(dataset.ab[,"fips"] %in% sample1), ], model = "twosteps", effect = "twoways")
fit.dab2 = pgmm(form.ab, dataset.ab[-which(dataset.ab[,"fips"] %in% sample1), ], model = "twosteps", effect = "twoways")
fit.ab = pgmm(form.ab, dataset.ab, model = "twosteps", effect = "twoways")

theta.hat.dab = 2*fit.ab$coefficients[[2]][1:(l+6)] - (fit.dab1$coefficients[[2]][1:(l+6)] + fit.dab2$coefficients[[2]][1:(l+6)])/2  # 6 is the number of treatment variables
HCV.coefs = vcovHC(fit.ab, cluster = 'group')[1:(l+6),1:(l+6)]
std.hat.dab = sqrt(diag(HCV.coefs))
stat.dab = theta.hat.dab/std.hat.dab

results = list(theta.hat.dab = theta.hat.dab, std.hat.dab = std.hat.dab, stat.dab = stat.dab, HCV.coefs = HCV.coefs)
save(results, file = "dab_covid_fat.dat")

load("dab_covid_fat.dat")
# short-run effects
round(results$theta.hat.dab,2)
round(results$std.hat.dab,2)
round(results$stat.dab,2)

# long-run effects of D
HCV.coefs = vcovHC(fit.ab, cluster = 'group')
lr.dab     = results$theta.hat.dab[1]/(1-sum(results$theta.hat.dab[2:(l+1)]))
jac.lr    = c(1,rep(lr.dab,l))/(1-sum(results$theta.hat.dab[2:(l+1)]))
cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs[c(1,2:(l+1)),c(1,2:(l+1))] %*% jac.lr)

round(lr.dab, 2)
round(cse.lr.dab, 2)
round(lr.dab/cse.lr.dab, 2) # T-stat

# long-run effects of others (C1)
lr.dab     = results$theta.hat.dab[8]/(1-sum(results$theta.hat.dab[2:(l+1)]))
jac.lr    = c(1,rep(lr.dab,l))/(1-sum(results$theta.hat.dab[2:(l+1)]))
cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs[c(8,2:(l+1)),c(8,2:(l+1))] %*% jac.lr)

round(lr.dab, 2)
round(cse.lr.dab, 2)
round(lr.dab/cse.lr.dab, 2) # T-stat

lr.dab     = results$theta.hat.dab[7]/(1-sum(results$theta.hat.dab[2:(l+1)]))
jac.lr    = c(1,rep(lr.dab,l))/(1-sum(results$theta.hat.dab[2:(l+1)]))
cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs[c(7,2:(l+1)),c(7,2:(l+1))] %*% jac.lr)

round(lr.dab, 2)
round(cse.lr.dab, 2)
round(lr.dab/cse.lr.dab, 2) # T-stat

lr.dab     = results$theta.hat.dab[9]/(1-sum(results$theta.hat.dab[2:(l+1)]))
jac.lr    = c(1,rep(lr.dab,l))/(1-sum(results$theta.hat.dab[2:(l+1)]))
cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs[c(9,2:(l+1)),c(9,2:(l+1))] %*% jac.lr)

round(lr.dab, 2)
round(cse.lr.dab, 2)
round(lr.dab/cse.lr.dab, 2) # T-stat

lr.dab     = results$theta.hat.dab[10]/(1-sum(results$theta.hat.dab[2:(l+1)]))
jac.lr    = c(1,rep(lr.dab,l))/(1-sum(results$theta.hat.dab[2:(l+1)]))
cse.lr.dab = sqrt(t(jac.lr) %*% HCV.coefs[c(10,2:(l+1)),c(10,2:(l+1))] %*% jac.lr)

round(lr.dab, 2)
round(cse.lr.dab, 2)
round(lr.dab/cse.lr.dab, 2) # T-stat
