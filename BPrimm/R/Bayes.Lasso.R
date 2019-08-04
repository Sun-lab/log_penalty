`Bayes.Lasso` <-
function(y, X, b=NULL, n.burn=2000, n.thin=20, n.iter=100, s=0.01, r=0.01,
  b.update.order=1, method=1)
{
  n = length(y)
  if(nrow(X) !=n){
    stop("dimensions of y and X do not match!\n")
  }
  p = ncol(X)

  if(method !=1 && method != 2){
    stop("method is neither 1 or 2!\n")
  }
  
  dims = c(n, p, n.burn, n.thin, n.iter, method, !is.null(b))
  bSample  = numeric(p*n.iter)
  kappaSam = numeric(n.iter)
  
  if(is.null(b)){
    b = 0
  }else{
    if(length(b) != p){
      stop("length of b is different from number of columns of X\n")
    } 
  }

  # in R, X is a matrix of p*n, in c, X is a matrix of n*p
  Z = .C("bayes_Lasso", as.double(y), as.double(t(X)), as.double(b),
     bSample = as.double(bSample), kappaSam = as.double(kappaSam),
     as.integer(dims), as.double(s), as.double(r),
     as.integer(b.update.order), PACKAGE="BPrimm")
  bSample = matrix(Z[["bSample"]], nrow=p, ncol=n.iter, byrow=TRUE)

  list(b=bSample, kappa=Z$kappaSam)
}
