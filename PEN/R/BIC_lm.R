BIC_lm <-
function(Betas,y,xx){
  n = nrow(xx)
  p = ncol(xx)
  
  if(length(Betas)!=(p + 1)){
      stop("The coefficient vector length not match!")
  }

  E = matrix(y,n,1) - Betas[1] - as.matrix(xx)%*%matrix(Betas[2:(p+1)],p,1) 
  df = length(which(Betas[2:(p+1)]!=0))
  sigma.e = var(E)*(n-1)/n
  bic = n*log(sigma.e) + df*log(n);
  return(bic)
}
