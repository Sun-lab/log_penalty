
library(BPrimm)
n = 3
maxit = 500
trace = 0
dim = c(n, maxit, trace)
xin = rep(0, n)
x   = rep(0, n)
Fmin = -9999.99
fail = 0

abstol = -1e8
reltol = sqrt(.Machine$double.eps)
alpha  = 1
beta   = 0.5
gamma  = 2

control = c(abstol, reltol, alpha, beta, gamma)

f <- function(b, bj.bar, wj, lambda){
  sum((b - bj.bar)^2/wj) + lambda*sqrt(sum(b*b))
}

fncount = 0

bj.bar = c(9.41e-03, -1.83e-03, -7.29e-03 )
wj = c(3.51e-03, 4.88e-03, 2.54e-03 )
lambda =1.50e+01

beta2 = beta3 = matrix(nrow=100, ncol=3)

for(i in 1:1e8){
  if(i %% 1e6 == 0){
    cat(i, date(), "\n")
  }
  bj.bar = runif(n, -1.0, 1.0)
  wj     = runif(n,  0.0, 1.0)
  lambda = runif(1,  0.0, 5.0)
  ex = c(bj.bar, wj, lambda)

  o1 = .C("solve_bj", as.integer(dim), as.double(xin),
    par = as.double(x), value=as.double(Fmin),
    convergence=as.integer(fail), as.double(control),
    counts=as.integer(fncount), as.double(ex), PACKAGE="BPrimm")
}

  o2 = optim(rep(0,3), fn=f, bj.bar=bj.bar, wj=wj, lambda=lambda)

  o3 = list()
  nms = c("par", "value", "convergence", "counts")
  for(nm in nms){
    o3[[nm]] = o1[[nm]]
  }
  
  beta2[i,] = o2$par
  beta3[i,] = o3$par
}
