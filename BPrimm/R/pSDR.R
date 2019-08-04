
pSDR <- function(y, X, method=c("SIR", "cov", "spline", "power"),
  covMat=c("ridge", "pls"), dimY=3, d=1, alpha=0.95,  pMax=20,
  bs.ord=3, bs.inknots=max(1, dimY - bs.ord+1), 
  taus=c(seq(0.1, 1, by=0.1), 1.2, seq(1.6, 3, by=0.4)), 
  u.max=15, threshold=2)
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
  
  # normalize a vector
  norm<-function(v)  
  { 
    sumv2 = sum(v^2)
    if(sumv2 == 0) sumv2 = 1
    v/sqrt(sumv2)
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
  
  B   = t(X) %*% yt / nrow(X)
  tau = u = NULL
  
  # --------------------------------------------------
  # if use ridge to handel singularity of cov matrix
  # --------------------------------------------------
  
  if(covMat == "ridge"){
    gcv = NULL
    St2 = mat.power(Sigma.x, 0.5)

    for(k in 1:length(taus)) {
      tau.k = taus[k]
      
      St1   = solve(Sigma.x + diag(tau.k, p))      
      theta = St1 %*% B 
      etah  = as.matrix(eigen(theta%*%t(theta))$vectors[,1:d])
      gamh  = t(etah) %*% theta
      
      comp1 = t(gamh) %*% solve(gamh%*%t(gamh)) %*% gamh
      comp2 = St2 %*% St1 %*% St2
      S     = comp1 %x% comp2
      
      comp3 = sum(((diag(1, p*h) - S) %*% as.vector(B))^2)
      comp4 = p*h*(1 - sum(diag(S))/(p*h))^2
      gcv   = c(gcv, comp3/comp4)
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
  
  # initial estimate, b.est = eta
  b.est = eigen(theta%*%t(theta))$vectors[,1:d, drop=FALSE]

  # form data for Lasso problem
  U = as.vector(theta)
  V = NULL
  
  for(s in 1:h) {
    # b.est is eta, and t(b.est) %*% theta[,s] is gamma
    V = rbind(V, diag(as.vector(b.est %*% (t(b.est) %*% theta[,s]))))
  }

  # solution path
  outLasso = glmnet(V, U, family="gaussian", alpha=alpha)
  
  w2kp  = which(outLasso$df < pMax & outLasso$df > 0)
  
  # BIC  
  df2=outLasso$df[w2kp]
  
  n.e = length(U)
  # for elastic net, outLasso$dev is R2
  RSS = (1 - outLasso$dev[w2kp])*var(U)*(n.e - 1)
  bic = n.e * log(RSS) + log(n.e) * df2
  w2use = w2kp[which.min(bic)]
  omega = outLasso$beta[,w2use]
  b.est = drop(apply(b.est * omega, 2, norm))

  if(is.matrix(b.est)){
    b.est = sqrt(rowSums(b.est*b.est))
  }
  
  ww4b = which(abs(b.est) > 1e-10)
  b2kp = b.est[ww4b]
  
  # return
  ans=list(w=ww4b, b=b2kp, method=method, covMat=covMat,
           dimY=h, bs.ord=bs.ord, bs.inknots=bs.inknots, 
           BIC=bic[which.min(bic)], tau=tau, u=u, d=d, yt = yt)
  
  class(ans) = "pSDR"
  
  return(ans)
}
