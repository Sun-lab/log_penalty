
`mGAP` <- function(y, X, method=c("SIR", "cov", "spline", "power"),
 dimY=3, bs.ord=3, bs.inknots=max(1, dimY - bs.ord+1), 
 lambda=seq(1, 15, by=1), 
 tau=c(seq(.01,.08,by=.01), 0.1, 0.2), recursive=TRUE,
 extbicgamma=1, bic.vec=TRUE, 
 pMax=20, L=10, epsilon=1e-4, n.max=500, b.update.order=1,
 p=NULL, offsetrow=0, offsetcol=0, transposeX=TRUE)
{

  if(!is.vector(y)){
    stop("y is not a vector\n")
  }

  n = length(y)
  if(is.matrix(X)){
    if(nrow(X) !=n ){
      stop("dimensions of y and X do not match!\n")
    }
    p = ncol(X)
  }else if(is.character(X)){
    if(file.access(X, mode = 0)!=0){
      stop(X, "is not a valid file name.\n")
    }
    if(is.null(p)){
      stop("when X is a file name, p must be specified\n")
    }
  }else{
    stop("X must be either a data matrix or a valid input file name.\n")
  }
  
  # First set the mean value of y as 0
  y = y - mean(y)
    
  if(dimY==1){
    yt = matrix(y, ncol=1)
  }else if(method=="SIR"){
    # ------------------------------------- #
    # method: SIR
    # ------------------------------------- #
    sy = dr.slices(y, dimY+1)
    yt = NULL

    for(s in 1:(dimY+1)) {
      Js = rep(0, n)
      ws = which(sy$slice.indicator == s)
      Js[ws] = 1/length(ws)
      yt = cbind(yt, Js)
    }
    
    w2drop = ceiling((dimY + 1)/2)
    yt = yt[,-w2drop]
  }else if(method=="cov"){
    # ------------------------------------- #
    # method: cov
    # ------------------------------------- #
    sy = dr.slices(y, dimY)
    yt = NULL

    for(s in 1:dimY) {
      Js = rep(0, n)
      Js[sy$slice.indicator == s] = 1
      y1 = Js*y
      yt = cbind(yt, y1)
    }
  }else if(method=="spline"){
    # ------------------------------------- #
    # method: spline
    # length(knots) - ord columns
    # that is bs.inknots + bs.ord columns
    # ------------------------------------- #
    knots = quantile(y, probs=seq(0, 1, length=(bs.inknots+2)))
    knots = c(rep(knots[1],bs.ord-1), knots, rep(knots[bs.inknots+2], bs.ord-1))
    yt    = splineDesign(knots, y, ord=bs.ord)
    yt    = yt[,-ncol(yt)]
  }else if(method=="power"){
    # ------------------------------------- #
    # method: power
    # ------------------------------------- #
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
  # preparation for C code
  # -------------------------------------------------------
  n.iter = 0

  score = matrix(0, nrow=length(lambda), ncol=length(tau))
  scoNA = matrix(1, nrow=length(lambda), ncol=length(tau))
  
  score2 = matrix(0, nrow=length(lambda), ncol=length(tau))
  scoNA2 = matrix(1, nrow=length(lambda), ncol=length(tau))
  
  score3 = matrix(0, nrow=length(lambda), ncol=length(tau))
  scoNA3 = matrix(1, nrow=length(lambda), ncol=length(tau))
  
  dims = c(n, p, h, L, n.max, offsetrow, offsetcol, transposeX)
  dims = c(dims, length(lambda), length(tau), pMax, recursive, bic.vec)
  
  score2use = lambda2use = tau2use = -9999.9999
  score2use2 = lambda2use2 = tau2use2 = -9999.9999
  score2use3 = lambda2use3 = tau2use3 = -9999.9999
  
  B = matrix(0, nrow=p, ncol=h)
  B2 = matrix(0, nrow=p, ncol=h)
  B3 = matrix(0, nrow=p, ncol=h)
  
  # ------------------------------------------------------- 
  # Call the C function
  # -------------------------------------------------------

  if(is.character(X)){
    Z = .C("mGAPc", yt=as.double(yt), as.character(X), B = as.double(t(B)), 
      as.double(lambda), as.double(tau), as.integer(dims),
      as.double(epsilon), n.iter=as.integer(n.iter),
      as.integer(b.update.order), score=as.double(t(score)),
      scoNA=as.integer(t(scoNA)), score2use=as.double(score2use), 
      lambda2use=as.double(lambda2use), tau2use=as.double(tau2use),
      B2 = as.double(t(B2)), score2=as.double(t(score2)),
      scoNA2=as.integer(t(scoNA2)), score2use2=as.double(score2use2), 
      lambda2use2=as.double(lambda2use2), tau2use2=as.double(tau2use2),
      B3 = as.double(t(B3)), score3=as.double(t(score3)),
      scoNA3=as.integer(t(scoNA3)), score2use3=as.double(score2use3), 
      lambda2use3=as.double(lambda2use3), tau2use3=as.double(tau2use3),
      extbicgamma=as.double(extbicgamma), 
      PACKAGE="BPrimm")
  }else{
    Z = .C("mGAPr", yt=as.double(yt), as.double(X), B = as.double(t(B)), 
      as.double(lambda), as.double(tau), as.integer(dims),
      as.double(epsilon), n.iter=as.integer(n.iter),
      as.integer(b.update.order), score=as.double(t(score)),
      scoNA=as.integer(t(scoNA)), score2use=as.double(score2use), 
      lambda2use=as.double(lambda2use), tau2use=as.double(tau2use),
      B2 = as.double(t(B2)), score2=as.double(t(score2)),
      scoNA2=as.integer(t(scoNA2)), score2use2=as.double(score2use2), 
      lambda2use2=as.double(lambda2use2), tau2use2=as.double(tau2use2),
      B3 = as.double(t(B3)), score3=as.double(t(score3)),
      scoNA3=as.integer(t(scoNA3)), score2use3=as.double(score2use3), 
      lambda2use3=as.double(lambda2use3), tau2use3=as.double(tau2use3),
      extbicgamma=as.double(extbicgamma), 
      PACKAGE="BPrimm")
  }
    
  score = matrix(Z$score, nrow=length(lambda), ncol=length(tau), byrow=TRUE)
  scoNA = matrix(Z$scoNA, nrow=length(lambda), ncol=length(tau), byrow=TRUE)
  score[scoNA==1] = NA
  
  score2 = matrix(Z$score2, nrow=length(lambda), ncol=length(tau), byrow=TRUE)
  scoNA2 = matrix(Z$scoNA2, nrow=length(lambda), ncol=length(tau), byrow=TRUE)
  score2[scoNA2==1] = NA
  
  score3 = matrix(Z$score3, nrow=length(lambda), ncol=length(tau), byrow=TRUE)
  scoNA3 = matrix(Z$scoNA3, nrow=length(lambda), ncol=length(tau), byrow=TRUE)
  score3[scoNA3==1] = NA
  
  if(all(is.na(as.vector(score)))){
    warning("no acceptable lambda/tau are found using BIC\n")
    return(NULL)
  }
  
  if(all(is.na(as.vector(score2)))){
    warning("no acceptable lambda/tau are found using extBIC\n")
    return(NULL)
  }
  
  if(all(is.na(as.vector(score3)))){
    warning("no acceptable lambda/tau are found using extBICgg\n")
    return(NULL)
  }
  
  
  bS      = matrix(Z$B, nrow=p, ncol=h, byrow=TRUE)
  bNom    = apply(bS, 1, function(v){sqrt(sum(v*v))})
  w2kp    = which(bNom > 1e-10)
  b2kp    = bS[w2kp,]
  lambda1 = Z$lambda2use
  tau1    = Z$tau2use

  bS2      = matrix(Z$B2, nrow=p, ncol=h, byrow=TRUE)
  bNom2    = apply(bS2, 1, function(v){sqrt(sum(v*v))})
  w2kp2    = which(bNom2 > 1e-10)
  b2kp2    = bS2[w2kp2,]
  lambda2  = Z$lambda2use2
  tau2     = Z$tau2use2
  
  bS3      = matrix(Z$B3, nrow=p, ncol=h, byrow=TRUE)
  bNom3    = apply(bS3, 1, function(v){sqrt(sum(v*v))})
  w2kp3    = which(bNom3 > 1e-10)
  b2kp3    = bS3[w2kp3,]
  lambda3  = Z$lambda2use3
  tau3     = Z$tau2use3
  
  rownames(score) = lambda
  colnames(score) = tau
  
  rownames(score2) = lambda
  colnames(score2) = tau
  
  rownames(score3) = lambda
  colnames(score3) = tau
  
  if(length(w2kp) > 1){ rownames(b2kp) = w2kp }
  
  if(length(w2kp2) > 1){ rownames(b2kp2) = w2kp2 }
  
  if(length(w2kp3) > 1){ rownames(b2kp3) = w2kp3 }
  
  
  ll = list(score_bic=score, w_bic=w2kp, b_bic=b2kp, 
            lambda_bic=lambda1, tau_bic=tau1, 
            score2use_bic = Z$score2use,
            score_extbic=score2, w_extbic=w2kp2, b_extbic=b2kp2, 
            lambda_extbic=lambda2, tau_extbic=tau2,
            score2use_extbic = Z$score2use2,
            score_extbicgg=score3, w_extbicgg=w2kp3, b_extbicgg=b2kp3, 
            lambda_extbicgg = lambda3, tau_extbicgg=tau3,
            score2use_extbicgg = Z$score2use3,
            method = method, dimY=h, bs.ord=bs.ord, bs.inknots=bs.inknots, 
            yt = yt)
  
  class(ll) = "mGAP"
  ll
}

