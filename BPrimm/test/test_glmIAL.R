
# =========================================================== 
# check glmIAL by simulated data
# ===========================================================

library(glmnet)
library(asSeq)

n = 500
p = 1000

X  = matrix(rnorm(n*p), nrow=n, ncol=p)
xT = c(11,22,88,99)
mu = rowSums(X[,xT])

y = numeric(length(mu))
for(i in 1:400){
  y[i] = rpois(1, exp(mu[i]))
}

date()

g1 = glmIAL("poisson", y, X, trace=1, pMax=50)

date()

g2 = glmnet(X, y, family="poisson")

date()


quartz(width=9, height=3)
par(mfrow=c(1,3), mar=c(5,4,1,1), cex=0.8)

plot(g1$wss)

plot(g1$w2use, g1$b2use, ylim=c(-0.2,1.1), ylab="Coef", xlab="glmIAL")
abline(v=xT)
abline(h=1.0, lty=2)

widx = which.min(abs(g2$df - length(g1$w2use)))
w2kp = which(g2$beta[,widx]!=0)
plot(w2kp, g2$beta[w2kp,widx], ylim=c(-0.2,1.1), ylab="Coef", xlab="glmnet")
abline(v=xT)
abline(h=1.0, lty=2)

