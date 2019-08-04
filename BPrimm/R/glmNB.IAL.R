glmNB.IAL <-
function(y,  X, delta=seq(0, 5.0, by=0.5), 
        tau = c(0.01*c(10,8,6,4,2,1), 0.001*c(5,2,1)), 
        pMax=20, offset=NULL, naPercent=0.4, 
        nTotal=NULL, maxit=20, maxitIAL=20, nReEstimate=pMax, 
        conv=1e-5, scoreTestP=0.05, trace=1)
{
  if(!is.numeric(y)){
    stop("y must be a numeric vector\n")
  }
  
  if(!is.matrix(X)){
    stop("X must be a vector\n")
  }
  
  M = ncol(X)
  N = length(y)
  
  if(nrow(X) != N){
    stop("the dimension of X and y do not match\n")
  }
  
  if(!is.null(offset)){
    useOffset = 1
    if((! is.numeric(offset)) || length(offset) != N){
      stop("offset must be a numeric vector of the same length as y\n")
    }
  }else{
    useOffset = 0
    offset    = 0.0 
  }

  isNA = apply(X, 1, function(v){ any(is.na(v)) })
  isNA = is.na(y) | isNA
  
  if(length(which(isNA))/N > naPercent){
    stop("percent of missing data is too high\n")
  }
  
  w2kp  = which(!isNA)
  yR    = y[w2kp]
  XR    = X[w2kp,]
  
  if(useOffset){
    offset = offset[w2kp]
  }

  N  = length(yR)
  init = 0
  
  dims    = numeric(10)
  dims[1] = N
  dims[2] = M
  dims[3] = maxit
  dims[4] = maxitIAL
  dims[5] = init
  dims[6] = useOffset
  dims[7] = length(delta)
  dims[8] = length(tau)
  dims[9] = pMax
  dims[10] = nReEstimate

  nIter   = 0
  phi     = 0.0
  
  w2use   = rep(-9, pMax)
  b2use   = numeric(pMax)
  score   = matrix(0, nrow=length(delta), ncol=length(tau))
  b0      = matrix(0, nrow=length(delta), ncol=length(tau))
  b02use  = 0.0
  convg   = 0
  n2use   = 0
  score2use = 0.0
  delta2use = 0.0
  tau2use   = 0.0
  family    = 0
  likelihood = rep(0, maxit)
  
  Z = .C("glmNB_IAL", as.integer(dims), 
         as.double(yR), as.double(offset), 
         as.double(XR), as.double(conv), as.double(scoreTestP), 
         as.integer(trace), as.double(delta), as.double(tau), 
         nIter=as.integer(nIter), phi=as.double(phi),  
         n2use=as.integer(n2use), w2use=as.integer(w2use), 
         b2use=as.double(b2use), b02use=as.double(b02use),
         b0=as.double(b0), score=as.double(score), 
         score2use=as.double(score2use), delta2use=as.double(delta2use), 
         tau2use=as.double(tau2use), likelihood=as.double(likelihood),
         family=as.integer(family), PACKAGE="asSeq")
  
  n2use=Z$n2use
  
  if(n2use > 0){
    w2use = Z$w2use[1:n2use] + 1
    b2use = Z$b2use[1:n2use]
  }else{
    w2use = b2use = NULL
  }
  
  if (Z$nIter >= maxit){
    warning(sprintf("reach mamximum %d iterations\n", maxit))  
  }
  
  score  = matrix(Z$score, nrow=length(delta), ncol=length(tau), byrow=TRUE)
  b0     = matrix(Z$b0,    nrow=length(delta), ncol=length(tau), byrow=TRUE)
  likelihood = Z$likelihood[1:Z$nIter]
  
  result = list(n2use=n2use, w2use=w2use, b2use=b2use, b0=b0, 
                score2use=Z$score2use, delta2use=Z$delta2use, 
                tau2use=Z$tau2use, score=score, nIter=Z$nIter, 
                family=Z$family, phi=Z$phi, likelihood=likelihood)
  
  result  
}

