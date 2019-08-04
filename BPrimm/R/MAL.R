
MAL <- function(y, X, method=c("SIR", "cov", "spline", "power"),
  covMat=c("ridge", "pls"), dimY=3, pMax=20, bs.ord=3, 
  u.max=15, threshold=2, bs.inknots=max(1, dimY - bs.ord+1), 
  taus=c(seq(0.5, 0.9, by=0.1), 1:15), 
  lambda0s=c(0.1, 0.2, 0.5, 1:10, 20, 50)*log(length(y)))
{
  
  # -------------------------------------------------------
  # utility functions
  # -------------------------------------------------------
  
  # power of a matrix
  mat.power <- function(A, a)
  {
    ei = eigen(A)
    d  = ei$values
    d  = (d+abs(d))/2
    d2  = d^a
    d2[d == 0] = 0
    ans = ei$vectors %*% diag(d2) %*% t(ei$vectors)
    return(ans)
  }
  
  # -------------------------------------------------------
  # check data
  # -------------------------------------------------------
  
  if(!is.vector(y)){
    stop("y is not a vector\n")
  }
    
  if(!is.matrix(X)){
    stop("X must be a data matrix.\n")
  }
  
  n = length(y)
  p = ncol(X)
  
  if(nrow(X) !=n ){
    stop("dimensions of y and X do not match!\n")
  }


  # First set the mean value of y as 0
  y = y - mean(y)
  
  if(dimY==1){
    yt = matrix(y, ncol=1)
  }else if(method=="SIR"){
    # --------------------------------------------------
    # method: SIR
    # --------------------------------------------------
    sy = dr.slices(y, dimY+1)
    yt = NULL
    
    for(s in 1:(dimY+1)) {
      Js = rep(0, n)
      ws = which(sy$slice.indicator == s)
      Js[ws] = n/length(ws)
      yt = cbind(yt, Js)
    }
    
    w2drop = ceiling((dimY + 1)/2)
    yt = yt[,-w2drop]
  }else if(method=="cov"){
    # --------------------------------------------------
    # method: cov
    # --------------------------------------------------
    sy = dr.slices(y, dimY)
    yt = NULL
    
    for(s in 1:dimY) {
      Js = rep(0, n)
      Js[sy$slice.indicator == s] = 1
      y1 = Js*y
      yt = cbind(yt, y1)
    }
  }else if(method=="spline"){
    # --------------------------------------------------
    # method: spline
    # (length(knots) - ord) columns
    # that is (bs.inknots + bs.ord) columns
    # --------------------------------------------------
    knots = quantile(y, probs=seq(0, 1, length=(bs.inknots+2)))
    knots = c(rep(knots[1],bs.ord-1), knots, rep(knots[bs.inknots+2], bs.ord-1))
    yt    = splineDesign(knots, y, ord=bs.ord)
    yt    = yt[,-ncol(yt)]
  }else if(method=="power"){
    # --------------------------------------------------
    # method: power
    # --------------------------------------------------
    yt = NULL
    for(s in 1:dimY) {
      yt = cbind(yt, y^s)
    }
  }else{
    stop("invalid method\n")
  }
  
  h  = ncol(yt)
  yt = scale(yt, center = TRUE, scale = TRUE)
  
  # -------------------------------------------------------
  # standardize X
  # -------------------------------------------------------
  
  X       = scale(X)
  Sigma.x = cov(X)
  
  # since we use Sigma.x, which is t(X)%*%X/n
  B   = t(X) %*% yt / n
  tau = u = NULL
  
  # --------------------------------------------------
  # if use ridge to handel singularity of cov matrix
  # --------------------------------------------------
  
  if(covMat == "ridge"){
    gcv = NULL

    for(k in 1:length(taus)) {
      tau.k = taus[k]
      
      St1   = solve(Sigma.x + diag(tau.k, p))      
      At1   = X %*% St1 %*% t(X)/n

      Bt1   = diag(n) - At1
      comp1 = Bt1 %*% y
      comp1 = sum(comp1*comp1)/n
      
      comp2 = sum(diag(Bt1))/n
      comp2 = comp2*comp2
      
      gcv   = c(gcv, comp1/comp2)
    }
      
    tau = taus[which.min(gcv)]

    theta = solve(Sigma.x + diag(tau, p)) %*% B
  }else{
    # power transformation components
    ei = eigen(Sigma.x)
    dx = ei$values
    dx = (dx+abs(dx))/2
    R  = nu = t(X) %*% y / nrow(X)
    
    for(k in 1:(u.max-1)){
      R = cbind(R, ei$vectors %*% diag(dx^k) %*% t(ei$vectors) %*% nu)
    }
    
    ev   = abs(eigen(R%*%t(R))$values)
    ev.r = ev/c(ev[-1], 1)
    ev.r[is.nan(ev.r)] = 1
    u = min(which(ev.r < threshold))
    R = R[, seq(1,u)]
    theta = R %*% mat.power(t(R) %*% Sigma.x %*% R, -1) %*% t(R) %*% B
  }
  
  # theta are the coefficients estiamte
  weight = sqrt(rowSums(theta*theta))
  
  # final set of coefficients
  coefsFinal = matrix(0, nrow=p, ncol=h)

  # after scaling, \|X\|^2 = n - 1
  xnorm  = n-1
  BICs   = rep(NA, length(lambda0s))
  minBIC = 1e16
  lambda2use = NA
  
  for(k in 1:length(lambda0s)){
    
    lambda0 = lambda0s[k]
    resids  = yt
    lambdas = lambda0/weight
    coefs   = matrix(0, nrow=p, ncol=h)
    coefOld = coefs
    
    for(it in 1:20){
      coefNorms = sqrt(rowSums(coefs*coefs))

      for(j in 1:p){
        
        if(coefNorms[j] > 0){
          resids    = resids + X[,j] %*% t(coefs[j,])
        }
        
        Sj = drop(t(resids) %*% X[,j])
        Sjnorm = sqrt(sum(Sj*Sj))
        
        diffj = Sjnorm - lambdas[j]/2
        
        if(diffj > 0){
          coefs[j,] = diffj*Sj/xnorm/Sjnorm
          resids    = resids - X[,j] %*% t(coefs[j,])
        }else{
          coefs[j,] = rep(0, h)
        }
        
      }
      
      maxDiff = max(as.vector(abs(coefOld - coefs)))
      # cat(it, maxDiff, "\n")
      
      if(maxDiff < 0.001){ next } 
      coefOld = coefs
    }
    
    coefNorms = sqrt(rowSums(coefs*coefs))
    ww4b = which(coefNorms > 1e-10)

    if(length(ww4b) > pMax || length(ww4b) == 0){ next }
    
    df = length(ww4b) + sum((h-1)*coefNorms[ww4b]/weight[ww4b])
    df = df/h
    BICs[k] = n*log(det(cov(resids))) + df*log(n)
    
    if(BICs[k] < minBIC){
      coefsFinal = coefs
      lambda2use = lambda0
      minBIC = BICs[k]
    }
    
  }
  
  coefNorms = sqrt(rowSums(coefsFinal*coefsFinal))
  ww4b = which(coefNorms > 1e-10)
  b2kp = coefNorms[ww4b]

  # return
  ans=list(w=ww4b, b=b2kp, method=method, covMat=covMat,
           dimY=h, bs.ord=bs.ord, bs.inknots=bs.inknots, 
           BIC=BICs, tau=tau, u=u, lambda0=lambda2use, 
           yt=yt)
  
  class(ans) = "MAL"
  
  return(ans)
}

