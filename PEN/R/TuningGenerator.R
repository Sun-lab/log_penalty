TuningGenerator <-
function(X, y,family ,Penalty, n.tau, n.lambdas, lambda.min, pfactor){
  
  output = list()
  output[["lambda"]] = NULL
  output[["tau"]] = NULL
  n = length(y)
  p = ncol(X)

  if (family=="gaussian"){

    r <- y - mean(y)
    l1.max <- max(abs(crossprod(X,r))) 
    
    if(Penalty == "SCAD" | Penalty == "MCP"){
      l1.max <- l1.max/n

      if (lambda.min==0) lambda <- c(exp(seq(log(l1.max),log(pfactor*l1.max),
                                             len=n.lambdas-1)),0)
      else lambda <- exp(seq(log(l1.max),log(pfactor*l1.max),len=n.lambdas))
      output[["lambda"]] = lambda
      
    }
  
    if(Penalty == "LOG"){
    
      maxTau = l1.max/n
      minLambda = l1.max*pfactor
      Thresholds <- c(exp(seq(log(l1.max),log(minLambda),len=n.lambdas)))
      tauset = c(exp(seq(log(1e-6),log(maxTau),len=n.tau)))

      tau = lambda = c()

      for(t in 1:length(Thresholds)){
        
        thres = Thresholds[t]  
        slambda = tauset*thres
  
        tau = c(tau, tauset)
        lambda = c(lambda, slambda)
      }

      if(length(tau)!=length(lambda)){
        stop("error")}
      
      output[["lambda"]] = lambda
      output[["tau"]] = tau
    }
    
    if(Penalty == "SICA"){
      maxTau = l1.max/n
      minLambda = l1.max*pfactor
      Thresholds <- c(exp(seq(log(l1.max),log(minLambda),len=n.lambdas)))
      
      tauset = c(exp(seq(log(1e-6),log(maxTau),len=n.tau)))
      tau = c()
      lambda = c()
      for(t in 1:length(Thresholds)){
        thres = Thresholds[t]
        
        for(tt in 1:length(tauset)){
          slambda = tauset[tt]^2*thres/(tauset[tt]*(tauset[tt]+1))
          tau = c(tau, tauset[tt])
          lambda = c(lambda, slambda)
        }
      }
      if(length(tau)!=length(lambda)){print("error")}
      
      output[["lambda"]] = lambda
      output[["tau"]] = tau

    }
  } 
  
  if (family=="binomial"){
    
    fit <- glm(y~1,family="binomial")
    pi. <- fit$fitted.values
    w <- pi.*(1-pi.)
    r = (y - pi.)/w
    l1.max <- max(abs(crossprod(X,w*r)))
    
    if(Penalty == "SCAD" | Penalty == "MCP"){
      
      l1.max <- l1.max/n
      
      if (lambda.min==0) lambda <- c(exp(seq(log(l1.max),log(pfactor*l1.max),len=n.lambdas-1)),0)
      else lambda <- exp(seq(log(l1.max),log(pfactor*l1.max),len=n.lambdas))
      output[["lambda"]] = lambda
      
    }
  
    if(Penalty == "LOG"){
      
      maxTau = l1.max/n
            
      Thresholds <- c(exp(seq(log(l1.max),log(pfactor*l1.max),len=n.lambdas)))
      tauset = c(exp(seq(log(1e-6),log(maxTau),len=n.tau)))

      tau = lambda = c()
      
      for(t in 1:length(Thresholds)){
        thres = Thresholds[t]
        slambda = thres*tauset
        tau = c(tau, tauset)
        lambda = c(lambda, slambda)
      }
      if(length(tau)!=length(lambda)){print("error")}

      output[["lambda"]] = lambda
      output[["tau"]] = tau
      
    }
    
    if(Penalty == "SICA"){
      maxTau = l1.max/n

      Thresholds <- c(exp(seq(log(l1.max),log(l1.max*pfactor),len=n.lambdas)))
      
      tauset = c(exp(seq(log(1e-6),log(maxTau),len=n.tau)))

      tau = c()
      lambda = c()
      
      for(t in 1:length(Thresholds)){
        thres = Thresholds[t]
        
        for(tt in 1:length(tauset)){
          slambda = tauset[tt]^2*thres/(tauset[tt]*(tauset[tt]+1))
          tau = c(tau, tauset[tt])
          lambda = c(lambda, slambda)
        }
        
      }
      if(length(tau)!=length(lambda)){print("error")}
      
      output[["lambda"]] = lambda
      output[["tau"]] = tau
    }
  }
 
  return(output)
}
