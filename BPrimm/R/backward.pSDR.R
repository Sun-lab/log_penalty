`backward.pSDR` <-
function(X, IA, pcut)
{
  if( class(IA) != "pSDR"){
    stop("IA must be of class 'pSDR'\n")
  }
  
  w = IA$w
  h = ncol(IA$yt)
  n = nrow(X)
  
  if(length(w)==0){
    return(IA)
  }
	
  bnorm = IA$b

  od = order(bnorm, decreasing=TRUE)
  ww = w[od]
  
  wuse = 0
  nw   = length(w)
    
  X2use = cbind(rep(1,n), X[,ww])
  
  for(i in nw:1){
    X1 = X2use[,1:(i+1)]
    X0 = X2use[,1:i]
    
    resid1 = resid0 = matrix(NA, nrow=n, ncol=h)
    
    for(s in 1:h){
      y = IA$yt[,s]
      
      l1 = lm(y ~ -1 + X1)
      l0 = lm(y ~ -1 + X0)
      
      wy = which(!is.na(y))
      resid1[wy,s] = l1$resid
      resid0[wy,s] = l0$resid      
    }
    
    det1 = det(cov(resid1, use="pair"))
    det0 = det(cov(resid0, use="pair"))
    lr   = n*(log(det0) - log(det1))
    if(lr < 0) stop("negative log likelihood ratio\n")
    
    pval = 1 - pchisq(lr, df=h)
    
    if(pval < pcut){
      wuse = i
      break
    }
  }
  
  IA1  = IA

  if(wuse > 0){
    wmin  = od[1:wuse]
    IA1$b = IA$b[wmin]
    IA1$w = IA$w[wmin]
  }else{
    IA1$b = numeric(0)
    IA1$w = numeric(0)
  }
    
  IA1
}
