PEN <-
function(X, y, family=c("gaussian","binomial"),
  penalty=c("SCAD","MCP","LOG","SICA"), gamma=3, lambda.min=NULL,
  n.lambdas=100, n.tau=6, lambda=NULL, tau=NULL, eps=.001,
  max.iter=1000, convex=TRUE, dfmax=NULL, ChooseXindex=NULL,
  pfactor=0.1, Model.selection="ExtendedBIC", restriction="none")
{
  
  if(! restriction %in% c("none", "non-negative")){
    stop("unexpected restriction\n")
  }
  
  restrictionR = as.integer(restriction == "non-negative")
  
  ## Index for non-penalized covariates
  if(is.null(ChooseXindex)){
      ChooseXindex = rep(1,ncol(X))      
  }else{
    if(length(ChooseXindex)!=ncol(X)){
        stop("error in length of ChooseXindex")
    }
  }
  
  ## Error checking
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (n.lambdas < 2) stop("n.lambdas must be at least 2")

  ## Set up XX, yy, lambda
  n <- length(y)
  p <- ncol(X)
  
  if(is.null(lambda.min)){
    lambda.min = ifelse(n>p,.001,.05) 
  }
  if(is.null(dfmax)){
    dfmax = p + 1
  }
  
  meanx <- apply(X,2,mean)
  normx <- sqrt(apply((t(X)-meanx)^2,1,sum)/n)
  if (any(normx < 0.0001)) stop("X contains columns which are numerically 
                                constant, please remove them; an intercept 
                                is included automatically")
  XX <- scale(X,meanx,normx)
  if (family=="gaussian") yy <- y - mean(y)
  else yy <- y
  if (is.null(lambda)) {
    lambda <- TuningGenerator(XX, y,family ,penalty, n.tau = n.tau, 
                              n.lambdas = n.lambdas, lambda.min = 
                              lambda.min, pfactor = pfactor)[["lambda"]]
    user.lambda <- FALSE
    
  } else {
    n.lambdas <- length(lambda)
    user.lambda <- TRUE
    
  }

  if (penalty == "LOG" | penalty == "SICA"){
    
    if (is.null(tau)){
      tau <- TuningGenerator(XX, y,family ,penalty, n.tau=n.tau,
                             n.lambdas = n.lambdas, lambda.min = lambda.min,
                             pfactor = pfactor)[["tau"]]
      user.tau <- FALSE
      
    }else{
      n.tau <- length(tau)
      user.tau <- TRUE
      
    }
  }
  
  ## Fit
  if (family=="gaussian"){
  
    fit <- .C("gaussian_model", double(p*length(lambda)),
              integer(length(lambda)), as.double(XX), as.double(yy),
              as.integer(n),  as.integer(p), penalty, as.double(lambda),
              as.double(tau), as.integer(length(lambda)),
              as.double(eps), as.integer(max.iter), as.double(gamma),
              as.integer(dfmax), as.integer(user.lambda),
              as.integer(ChooseXindex), as.integer(restrictionR))
              
    beta <- rbind(0,matrix(fit[[1]],nrow=p))
    iter <- fit[[2]]
  }
  
  if (family=="binomial"){
    fit <- .C("binomial_model", double(length(lambda)), double(p*length(lambda)),
              integer(length(lambda)), as.double(XX), as.double(yy),
              as.integer(n), as.integer(p), penalty, as.double(lambda),
              as.double(tau), as.integer(length(lambda)), as.double(eps),
              as.integer(max.iter), as.double(gamma), as.integer(dfmax),
              as.integer(user.lambda), as.integer(ChooseXindex),
              as.integer(restrictionR))
              
    beta <- rbind(fit[[1]],matrix(fit[[2]],nrow=p))
    iter <- fit[[3]]
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(beta[p,])
  
  if(length(ind)>0){
    beta <- beta[,ind,drop=FALSE]
    iter <- iter[ind]
    lambda <- lambda[ind]
    tau <- tau[ind]}
  
    
  #if (any(iter==max.iter)){
  # warning("Algorithm failed to converge for all values of lambda")
  #}
  
  ## Eliminate models which are not converged
  
  keep.index = which(iter!=max.iter)
  
  if(length(keep.index)==0){
    val= list(beta=NULL)
  }else{
    beta = beta[,keep.index,drop=FALSE]
    iter = iter[keep.index]
    lambda = lambda[keep.index]
    tau = tau[keep.index]
  
    convex.min <- NULL
      
    if(penalty=="MCP" || penalty=="SCAD"){ 
      if (convex) convex.min <- convexMin(beta, XX, penalty, gamma, 0, family)
      else convex.min <- NULL
    }else{
      convex.min=NULL
    }

    ## Unstandardize
    beta[-1,] <- beta[-1,]/normx
    
    if (family=="gaussian") {
      beta[1,] <- mean(y) - crossprod(meanx,beta[-1,,drop=FALSE])
    }
    
    if (family=="binomial") {
      beta[1,] <- beta[1,] - crossprod(meanx,beta[-1,,drop=FALSE])
    }

    ## Names
    if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
    else varnames <- colnames(X)
    varnames <- c("(Intercept)",varnames)
    dimnames(beta) <- list(varnames,round(lambda,digits=4))
      
      
      
    ## Model selection   
    if(Model.selection=="ExtendedBIC"){
      if(family=="gaussian"){
        extbics = apply(beta,2,ExtendedBIC_lm, y=y, xx=X)
      }
      
      if(family=="binomial"){
        extbics = apply(beta,2,ExtendedBIC_glm, y=y, xx=X)
      }
      
      if(penalty=="MCP" || penalty=="SCAD"){
        tau.output = NULL
      }else{
          tau.output = tau[which.min(extbics)]
      }
      
      val <- list(beta=beta[,which.min(extbics)], iter=iter[which.min(extbics)],
                  lambda.output=lambda[which.min(extbics)],
                  tau.output = tau.output, penalty=penalty,
                  family=family, gamma=gamma, convex.min=convex.min)

    }else if(Model.selection=="BIC"){
      if(family=="gaussian"){
        bics = apply(beta,2,BIC_lm, y=y, xx=X)
      }
      
      if(family=="binomial"){
        bics = apply(beta,2,BIC_glm, y=y, xx=X)
      }
      
      if(penalty=="MCP" || penalty=="SCAD"){
        tau.output = NULL
      }else{
        tau.output = tau[which.min(bics)]
      }
      
      val <- list(beta=beta[,which.min(bics)], iter=iter[which.min(bics)],
                  lambda.output=lambda[which.min(bics)],
                  tau.output = tau.output, penalty=penalty,
                  family=family, gamma=gamma, convex.min=convex.min)
      
    }else if(Model.selection=="None"){
      val <- list(beta=beta, iter=iter, lambda=lambda, tau=tau, penalty=penalty,
                  family=family, gamma=gamma, convex.min=convex.min)
    }else{
      stop("wrong value for Model.selection!\n")
    }

    ## Output
    class(val) <- "PenalizedEstimation"
  }
  
  return(val)
}
