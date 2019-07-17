ne.PEN <-
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
    
    stopifnot(is.matrix(dat),nrow(dat)>1,ncol(dat)>1)
    p = ncol(dat)
    n = nrow(dat)
    coefNI = matrix(0, nrow=p, ncol=length(V))
    
    for (i in 1:length(V)) {
        v = V[i]  
        if (verbose) cat("variable=",v," ")
        X = dat[,-v]
        y = dat[,v]
        corXy = abs(drop(cor(X,y)))
        
        # --------------------------------------------------------------
        # estimate skeleton by penalized regression
        # --------------------------------------------------------------
        if (order) {
            o = order(corXy,decreasing=T)
            o.back = order(o)
            wp     = PEN(X[,o], y, family="gaussian", penalty="LOG", 
                         n.lambdas=nlambda, n.tau = ntau, 
                         Model.selection=Model.selection)
            
            wpBetas = wp$beta[-1]
            coefNI[-v,i] = wpBetas[o.back]
            if (verbose) cat(date(), "\n")
        }else {
            wp = PEN(X, y, family="gaussian", penalty="LOG", 
                     n.lambdas=nlambda, n.tau = ntau, 
                     Model.selection=Model.selection)
          
          coefNI[-v,i] = wp$beta[-1]
          if (verbose) cat( date(), "\n")
        }
    }
    return(coefNI)
}
