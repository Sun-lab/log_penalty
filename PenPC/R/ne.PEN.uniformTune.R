ne.PEN.uniformTune <-
  function(dat,nlambda,ntau,V,order=FALSE,verbose=FALSE, 
           Model.selection="ExtendedBIC")
  {
    ## Purpose: Fitting Gaussian Graphical Model using log penalty
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## - dat : n x p data matrix with n (# of sample) and p (# of variables)
    ## - nLambda : # of candidate lambda as a tuning parameter of Log penalty
    ## - nTau : # of candidate tau as a tuning parameter of Log penalty
    ## - V : the set of vertices which neighbors are estimated
    ## - order : if TURE covariate order are in the order of marginal correlation
    ## ----------------------------------------------------------------------
    ## Author: Min Jin Ha, Date: 18 March 2013
    ## Modified by Wei Sun, Date: 4 April 2017
    
    stopifnot(is.matrix(dat),nrow(dat)>1,ncol(dat)>1)
    p = ncol(dat)
    n = nrow(dat)
    
    # --------------------------------------------------------------
    # standardize data
    # --------------------------------------------------------------
    
    meanx   = apply(dat, 2, mean)
    normx   = sqrt(rowSums((t(dat)-meanx)^2)/n)
    nDat    = scale(dat, meanx, normx)
    
    # --------------------------------------------------------------
    # Determine the range of tuning parameters
    # --------------------------------------------------------------
    
    corr = abs(crossprod(nDat,nDat))
    diag(corr) = 0
    
    lamMax  = max(corr)
    thresh  = 2*exp(seq(log(lamMax), log(1/n), len=nlambda))
    lambda  = thresh/10
    tau     = 10^(seq(-6, -1, length.out=ntau))
    
    lambda  = rep(lambda, each=ntau)
    tau     = rep(tau, times=nlambda)
    
    key0 = sprintf("%.5e_%.5e", lambda, tau)
    
    if(Model.selection=="None"){
      coefNI = array(dim=c(p, length(V), length(lambda)))
    }else{
      coefNI = matrix(0, nrow=p, ncol=length(V))
    }

    for (i in 1:length(V)) {
      v = V[i]  
      if (verbose) cat("variable=",v," ")
      X = nDat[,-v]
      y = nDat[,v]
      corXy = abs(drop(cor(X,y)))
      
      # --------------------------------------------------------------
      # estimate skeleton by penalized regression
      # --------------------------------------------------------------
      
      if (order) {
        o = order(corXy,decreasing=T)
        o.back = order(o)
        wp     = PEN(X[,o], y, family="gaussian", penalty="LOG", 
                      lambda=lambda, tau=tau, Model.selection=Model.selection)
        
        key1 = sprintf("%.5e_%.5e", wp$lambda, wp$tau)
        
        if(Model.selection=="None"){
          wpBetas = wp$beta[-1,]
          coefNI[-v,i,match(key1, key0)] = wpBetas[o.back,]
        }else{
          wpBetas = wp$beta[-1]
          coefNI[-v,i] = wpBetas[o.back]
        }
        
        if (verbose) cat(date(), "\n")
      }else {
        wp  = PEN(X, y, family="gaussian", penalty="LOG", 
                  lambda=lambda, tau = tau, Model.selection=Model.selection)
        
        key1 = sprintf("%.5e_%.5e", wp$lambda, wp$tau)
        
        if(Model.selection=="None"){
          coefNI[-v,i,match(key1, key0)] = wp$beta[-1,]
        }else{
          coefNI[-v,i] = wp$beta[-1]
        }
        
        if (verbose) cat( date(), "\n")
      }
    }
    return(coefNI)
  }
