ExtendedBIC_glm <-
function(Betas,y,xx){
  n = nrow(xx)
  p = ncol(xx)
  int = Betas[1]
  linear = int + as.matrix(xx)%*%matrix(Betas[2:(p+1)],p,1) 
  pi = exp(linear)/(1 + exp(linear))
  
  loglik = 0
  for(k in 1:n){
    loglik = loglik + y[k]*log(pi[k]) + (1 - y[k])*log(1 - pi[k])
  }
  
  
  df = length(which(Betas[2:(p+1)]!=0))
  bic = -2*loglik + df*log(n) + 2*df*0.5*log(p)
  return(bic)
}
