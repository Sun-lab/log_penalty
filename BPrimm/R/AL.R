`AL` <-
function(y, X, b=NULL, lambda=NULL, gamma=seq(0.5,2,by=0.5), pMax=20, alpha=0.95){
  
  n = length(y)
  if(nrow(X) !=n){
    stop("dimensions of y and X do not match!\n")
  }
  p = ncol(X)  
  
  one   = rep(1, n)
  meanx = drop(one %*% X)/n
  X     = scale(X, meanx, FALSE)
  mu    = mean(y)
  y     = drop(y - mu)

  if(is.null(b)){
    for(i in 1:p){
      l = lm(y~X[,i])
      b[i] = l$coef[2]
    }
  }else{
    if(length(b) != p){
      stop("length of b is different from number of columns of X\n")
    }
  }
  
  if(length(alpha) > 1){ stop("only scalar alpha are allowed\n") }
  # ------------------------------------------------------- 
  # If there is only one lambda and one gamma, 
  # the parameter selection is not needed
  # -------------------------------------------------------

  if(length(lambda) == 1 && length(gamma) == 1){
    lambda1 = lambda[1]
    gamma1  = gamma[1]
  }else{
    # ------------------------------------------------------- 
    # Use BIC to select variables
    # -------------------------------------------------------
    lam2Use = BIC = rep(NA, length(gamma))
    
    for(k2 in 1:length(gamma)){
      
      gamma1 = gamma[k2]
      ws     = abs(b)^gamma1
      Xn     = t(t(X) * ws)
      
      if(is.null(lambda)){
        gn = glmnet(Xn, y, family="gaussian", alpha=alpha)
      }else{
        gn = glmnet(Xn, y, family="gaussian", alpha=alpha, lambda=lambda)
      }
      
      wk1  = which(gn$df > 0 & gn$df <= pMax)
      BIC1 = rep(NA, length(gn$lambda))

      for(k1 in wk1){
        lambda1 = gn$lambda[k1]
        bhat    = gn$beta[,k1]*ws
      
        w2kp = which(abs(bhat) > 1e-10)
        b2kp = bhat[w2kp]
        
        if(length(w2kp)==1){
          resd = y - X[,w2kp] * b2kp
        }else{
          resd = y - X[,w2kp] %*% b2kp
        }
        BIC1[k1] = log(sum(resd*resd)/n) + length(w2kp)*log(n)/(n)        
      }
      
      if(! all(is.na(BIC1))){
        wMin = which.min(BIC1)
        BIC[k2] = BIC1[wMin]
        lam2Use[k2] = gn$lambda[wMin]
      }
      
    }
    
    if(all(is.na(as.vector(BIC)))){
      warning("no acceptable lambda/gamma are found\n")
      return(NULL)
    }
    
    wMin = which.min(BIC)
    lambda1 = lam2Use[wMin]
    gamma1  = gamma[wMin]
    
  }
  
  ws = abs(b)^gamma1
  Xn = t(t(X) * ws)
  
  gn = glmnet(Xn, y, family="gaussian", alpha=alpha, lambda=lambda1)
  bhat = gn$beta[,1]*ws

  w2kp = which(abs(bhat) > 1e-10)
  b2kp = bhat[w2kp]

  ll = list(BIC=BIC, w=w2kp, b=b2kp, lambda=lambda1, gamma=gamma1)
  class(ll) = "AL"
  return(ll)  
}

