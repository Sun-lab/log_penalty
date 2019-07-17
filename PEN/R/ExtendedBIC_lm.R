ExtendedBIC_lm <-
function(Betas,y,xx){
  n = nrow(xx)
  p = ncol(xx)
  
  if(length(Betas)!=(p + 1)){
      stop("The coefficient vector length not match!")
  }

  BIC.gamma = 1 - 1/(2*log(p)/log(n))
  E = matrix(y,n,1) - Betas[1] - as.matrix(xx)%*%matrix(Betas[2:(p+1)],p,1) 
  df = length(which(Betas[2:(p+1)]!=0))
  sigma.e = var(E)*(n-1)/n
  tmp = 2*BIC.gamma*lchoose(p, df);
  bic = n*log(sigma.e) + df*log(n) + tmp;
  return(bic)
}
