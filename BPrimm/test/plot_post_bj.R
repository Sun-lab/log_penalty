sigma.j2 = 1.0
kappa.j  = 0.1
bj.bar   = 0.8

library(mlm)

n = 1000
bj = numeric(n)
# in R, X is a matrix of p*n, in c, X is a matrix of n*p
Z = .C("bj_sampler", as.integer(n), as.double(bj.bar),
  as.double(sigma.j2), as.double(kappa.j),
  bj = bj, PACKAGE="mlm")

hist(Z$bj, freq=FALSE, breaks=seq(-1,1,by=0.05))

bj = seq(-1, 1, by=0.01)
dy = exp(-(bj - bj.bar)^2/(2*sigma.j2) - abs(bj)/kappa.j)
plot(bj, 4*dy, type="l")

bj = seq(-1, 1, by=0.01)
dz = -(bj - bj.bar)^2/(2*sigma.j2) - abs(bj)/kappa.j
plot(bj, dz, type="l")

