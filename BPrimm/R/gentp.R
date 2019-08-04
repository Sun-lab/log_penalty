`gentp` <- function(y, X, method=c("SIR", "cov", "spline", "power"),
                    dimY=3, bs.ord=3, bs.inknots=max(1, dimY - bs.ord+1), 
                    ntau = 5, nlambda = 50, lambda_min = 1e-6){
  
  if(!is.vector(y)){
    stop("y is not a vector\n")
  }
  
  n <- length(y)
  if(is.matrix(X)){
    if(nrow(X) !=n ){
      stop("dimensions of y and X do not match!\n")
    }
    p <- ncol(X)
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
  y <- y - mean(y)
  
  if(dimY == 1){
    yt <- matrix(y, ncol = 1)
  }else if(method == "SIR"){
    # ------------------------------------- #
    # method: SIR
    # ------------------------------------- #
    sy <- dr.slices(y, dimY+1)
    yt <- NULL
    
    for(s in 1:(dimY+1)) {
      Js <- rep(0, n)
      ws <- which(sy$slice.indicator == s)
      Js[ws] <- 1/length(ws)
      yt <- cbind(yt, Js)
    }
    
    w2drop <- ceiling((dimY + 1)/2)
    yt <- yt[,-w2drop]
  }else if(method == "cov"){
    # ------------------------------------- #
    # method: cov
    # ------------------------------------- #
    sy <- dr.slices(y, dimY)
    yt <- NULL
    
    for(s in 1:dimY) {
      Js <- rep(0, n)
      Js[sy$slice.indicator == s] <- 1
      y1 <- Js*y
      yt <- cbind(yt, y1)
    }
  }else if(method == "spline"){
    # ------------------------------------- #
    # method: spline
    # length(knots) - ord columns
    # that is bs.inknots + bs.ord columns
    # ------------------------------------- #
    knots <- quantile(y, probs=seq(0, 1, length=(bs.inknots+2)))
    knots <- c(rep(knots[1],bs.ord-1), knots, rep(knots[bs.inknots+2], bs.ord-1))
    yt <- splineDesign(knots, y, ord=bs.ord)
    yt <- yt[,-ncol(yt)]
  }else if(method == "power"){
    # ------------------------------------- #
    # method: power
    # ------------------------------------- #
    yt <- NULL
    for(s in 1:dimY) {
      yt <- cbind(yt, y^s)
    }
  }else{
    stop("invalid method\n")
  }
  
  h <- ncol(yt)
  yt <- scale(yt, center = TRUE, scale = TRUE)
  
  X <- scale(X)
  
  tau_grid <- 10^(seq(-6, -1, length.out = ntau)) 
  lambda_max_grid <- vector()
  
  b_bar_store = list()
  # for each tau, calculate the max lambda 
  # calculate max lambda for given tau
  for(tau_iter in 1:length(tau_grid)){
    kappa <- tau_grid[tau_iter]
    eta <- matrix(nrow = ncol(X), ncol = dimY)
    b_bar_store[[tau_iter]] = matrix(nrow = ncol(X), ncol = dimY)
    
    for(j in 1:ncol(X)){
      Xj <- X[, j, drop = FALSE]
      for(s in 1:dimY){
        Ts <- yt[, s, drop=FALSE]
        b_bar <- (t(Xj)%*%Xj)^-1 * t(Xj)%*%Ts
        b_bar_store[[tau_iter]][j,s] <- b_bar
        omega <- (norm(Ts)^2/(n*(norm(Xj)^2)))
        eta[j,s] <- (2*b_bar*kappa)/omega
      }
    }
    
    eta_norm <- apply(eta, 1, norm)
    #for the first iteration, we assume b = 0
    
    lambda_max <- max(eta_norm)
    lambda_max_grid[tau_iter] <- lambda_max
  }
  
  lambda_grid <- exp(seq(log(lambda_min), 
                         log(max(lambda_max_grid)), 
                         length.out = nlambda))
  
  res <- list(lambda = lambda_grid, tau = tau_grid)
  res
}

