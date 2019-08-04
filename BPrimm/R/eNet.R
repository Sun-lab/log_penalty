`eNet` <-
function(y, X, lambda=NULL, pMax=20, alpha=0.95, backward.pval=NULL){
  
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

  # ------------------------------------------------------- 
  # If there is only one lambda and one gamma, 
  # the parameter selection is not needed
  # -------------------------------------------------------
  
  BICmin = 1e16
  alpha2use = NA
  
  BIC2use = NULL
  eN3 = NULL
  
  for(alpha1 in alpha){
    # ------------------------------------------------------- 
    # Use BIC to select variables
    # -------------------------------------------------------
    
    if(is.null(lambda)){
      gn = glmnet(X, y, family="gaussian", alpha=alpha1)
    }else{
      gn = glmnet(X, y, family="gaussian", alpha=alpha1, lambda=lambda)
    }
    
    wk1 = which(gn$df > 0 & gn$df <= pMax)
    BIC = rep(NA, length(gn$lambda))

    for(k1 in wk1){
      lambda1 = gn$lambda[k1]
      bhat    = gn$beta[,k1]
    
      w2kp = which(abs(bhat) > 1e-10)
      b2kp = bhat[w2kp]
      
      if(length(w2kp)==1){
        resd = y - X[,w2kp] * b2kp
      }else{
        resd = y - X[,w2kp] %*% b2kp
      }
      BIC[k1] = log(sum(resd*resd)/n) + length(w2kp)*log(n)/(n)        
    }
    
    if(all(is.na(BIC))){
      warning("no acceptable lambda is found\n")
      return(NULL)
    }
    
    wMin = which.min(BIC)
    
    bhat = gn$beta[,wMin]
    w2kp = which(abs(bhat) > 1e-10)
    b2kp = bhat[w2kp]

    eN1 = list(BIC=BIC, w=w2kp, b=b2kp, lambda=gn$lambda[wMin], alpha=alpha1)
    class(eN1) = "eNet"

    BIC1 = BIC[wMin]
    if(!is.null(backward.pval)){
      eN2 = backward(y, X, eN1, pcut=backward.pval)
      BIC1 = BIC1 + (length(eN2$w) - length(eN1$w))*log(n)/n 
    }
    
    BIC2use = c(BIC2use, BIC1)

    if(BIC1 < BICmin){
      BICmin = BIC1
      eN3    = eN2
    }
  }

  eN3$BIC = BIC2use
  eN3
}

