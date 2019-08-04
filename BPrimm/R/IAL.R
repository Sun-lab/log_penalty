`IAL` <-
function(y, X, b=NULL, delta=seq(0, 5.0, by=0.5), 
 tau = c(1e-3, 5e-3, 0.01, 0.02, 0.05, 0.1, 0.2), 
 criterion="BIC", cv=10, times.cv=10, BIC.gamma=0, pMax=20, 
 L=10, epsilon=1e-4, n.max=1000, b.update.order=1, 
 p=NULL, offsetrow=0, offsetcol=0, transposeX=FALSE, Rratio=1)
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
  
  if(is.null(BIC.gamma)){
    BIC.gamma = 1 - 1/(2*log(p)/log(n))
  }
  if(BIC.gamma > 1 | BIC.gamma < 0){
    stop("BIC.gamma must in the range of [0,1]\n")
  }
  
  if(is.null(b)){
    b = numeric(p)
  }else{
    if(length(b) != p){
      stop("length of b is different from number of columns of X\n")
    } 
  }

  if(criterion[1] == "cv"){
    if(cv%%1 != 0 || cv <= 1){
      stop("cv must be a positive integer larger than 1\n")
    }
  }
    
  # ------------------------------------------------------- 
  # preparation for C code
  # -------------------------------------------------------
  n.iter  = 0
  bSample = numeric(p)
  b0      = 0

  score = matrix(0, nrow=length(delta), ncol=length(tau))
  scoNA = matrix(1, nrow=length(delta), ncol=length(tau))
  
  # ------------------------------------------------------- 
  # If use cross-validation to select variables
  # idc is indicator of the cv groups
  # -------------------------------------------------------
  if(cv==n){
    idc = matrix(c(0:(n-1)),n,1)
    times.cv = 1
  }else{
    sam = sample(n)
    nb  = n %% cv
    na  = cv - nb
    nk  = floor(n/cv)
    nob = c(rep(nk,na), rep(nk+1,nb))
    
    indc = rep(0:(cv-1), times=nob)
    idc  = matrix(NA, n, times.cv)

    for (s in 1:times.cv){
      idc[,s] = sample(indc)
    }
  }
  
  dims = c(n, p, L, n.max, !is.null(b), offsetrow, offsetcol, transposeX)
  dims = c(dims, length(delta), length(tau), pMax, cv, times.cv)
  
  score2use = delta2use = tau2use = -9999.9999
    
  # ------------------------------------------------------- 
  # Call the C function
  # -------------------------------------------------------

  if(is.character(X)){
    Z = .C("IALc", as.double(y), as.character(X), as.double(b),
      bSample = as.double(bSample), b0 = as.double(b0),
      as.double(delta), as.double(tau), as.integer(dims),
      as.double(epsilon), n.iter=as.integer(n.iter),
      as.integer(b.update.order), score=as.double(t(score)),
      scoNA=as.integer(t(scoNA)), as.double(BIC.gamma), 
      as.integer(t(idc)), as.character(criterion), 
      score2use=as.double(score2use), delta2use=as.double(delta2use),
      tau2use=as.double(tau2use),Rratio=as.double(Rratio), PACKAGE="BPrimm")
  }else{
    Z = .C("IALr", as.double(y), as.double(t(X)), as.double(b), 
      bSample = as.double(bSample), b0 = as.double(b0),
      as.double(delta), as.double(tau), as.integer(dims),
      as.double(epsilon), n.iter=as.integer(n.iter),
      as.integer(b.update.order), score=as.double(t(score)),
      scoNA=as.integer(t(scoNA)), as.double(BIC.gamma), 
      as.integer(t(idc)), as.character(criterion), 
      score2use=as.double(score2use), delta2use=as.double(delta2use),
      tau2use=as.double(tau2use),Rratio=as.double(Rratio), PACKAGE="BPrimm")
  }
    
  score = matrix(Z$score, nrow=length(delta), ncol=length(tau), byrow=TRUE)
  scoNA = matrix(Z$scoNA, nrow=length(delta), ncol=length(tau), byrow=TRUE)
  score[scoNA==1] = NA
  
  if(all(is.na(as.vector(score)))){
    warning("no acceptable delta/tau are found\n")
    return(NULL)
  }
    
  bS = Z$bSample
  w2kp = which(abs(bS) > 1e-10)
  b2kp = bS[w2kp]
  delta1 = Z$delta2use
  tau1   = Z$tau2use
  if(Z$n.iter>=n.max){warning("The algorithm doesn't converge\n")}
  
  ll = list(score=score, b0=Z$b0, w=w2kp, b=b2kp, delta=delta1, tau=tau1)
  ll[["score2use"]] = Z$score2use
  class(ll) = "IAL"
  ll
}
