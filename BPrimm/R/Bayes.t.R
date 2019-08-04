`Bayes.t` <-
function(y, X, b=NULL, n.burn=2000, n.thin=20, n.iter=100,
  delta=NULL, tau=NULL, b.update.order=1)
{
  n = length(y)
  if(nrow(X) !=n){
    stop("dimensions of y and X do not match!\n")
  }
  p = ncol(X)

  dims = c(n,p,n.burn,n.thin,n.iter,is.null(delta),is.null(tau),!is.null(b))
  bSample = numeric(p*n.iter)
  deltaSam = tauSam = numeric(n.iter)

  if(is.null(delta)){
    delta1 = 1.0
  }else{
    delta1 = delta
  }

  if(is.null(tau)){
    tau1 = 0.01
  }else{
    tau1 = tau
  }

  if(is.null(b)){
    b = 0
  }else{
    if(length(b) != p){
      stop("length of b is different from number of columns of X\n")
    } 
  }

  Z = .C("bayes_t", as.double(y), as.double(t(X)), as.double(b),
     bSample = as.double(bSample), deltaSam = as.double(deltaSam),
     tauSam = as.double(tauSam), as.integer(dims), as.double(delta1), 
     as.double(tau1), as.integer(b.update.order), PACKAGE="BPrimm")
     
  bSample = matrix(Z[["bSample"]], nrow=p, ncol=n.iter, byrow=TRUE)
  
  if(is.null(delta)){ delta1 = Z$deltaSam }
  if(is.null(tau)){ tau1 = Z$tauSam }

  list(b=bSample, delta=delta1, tau=tau1)
}
