
# =========================================================== 
# check glmIAL and glmNB.IAL by simulated data
# ===========================================================

library(MASS)
library(glmnet)
library(asSeq)

n = 500
p = 500

X  = matrix(rnorm(n*p), nrow=n, ncol=p)
xT = c(11,22,88,99)
mu = rowSums(X[,xT])

y = rnegbin(exp(mu), theta=1.0)

date()

l1 = glm.nb(y ~ X)

date()

g1 = glmIAL("poisson", y, X, trace=1, pMax=20)

date()

g2 = glmnet(X, y, family="poisson")

date()

g3 = glmNB.IAL(y, X, trace=1, pMax=20)

date()

par(mfrow=c(2,1))
plot(g1$wss)
plot(g3$likelihood)


par(mfrow=c(2,2), mar=c(5,4,1,1), cex=0.8)

plot(l1$coef, ylim=c(-0.2,1.1), ylab="Coef", xlab="glm.nb")
abline(v=xT+1)
abline(h=1.0, lty=2)

plot(g1$w2use, g1$b2use, ylim=c(-0.2,1.1), ylab="Coef", xlab="glmIAL")
abline(v=xT)
abline(h=1.0, lty=2)

widx = which.min(abs(g2$df - length(g1$w2use)))
w2kp = which(g2$beta[,widx]!=0)
plot(w2kp, g2$beta[w2kp,widx], ylim=c(-0.2,1.1), ylab="Coef", xlab="glmnet")
abline(v=xT)
abline(h=1.0, lty=2)

plot(g3$w2use, g3$b2use, ylim=c(-0.2,1.1), ylab="Coef", xlab="glmNB.IAL")
abline(v=xT)
abline(h=1.0, lty=2)


