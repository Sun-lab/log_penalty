MajorityVote_PEN = function(X, y, Xinfo, UpdatingNum = 20,family, penalty, gamma=3, lambda.min=NULL, 
  n.lambdas=100, n.tau = 6, lambda = NULL, tau = NULL, eps=.001, max.iter=1000, convex=TRUE, 
  dfmax=NULL, ChooseXindex= NULL, pfactor = 0.1,Model.selection="Extended BIC"){
            
  if(length(Xinfo)!= ncol(X)){
    stop("Xinfo dim incorrect")
  }
            
  ## Set up XX, yy, lambda
  n <- length(y)
  p <- ncol(X)
            
  if(is.null(lambda.min)){
    lambda.min = ifelse(n>p,.001,.05) 
  }
  if(is.null(dfmax)){
    dfmax = p + 1
  }
            
  Chrs = sort(unique(Xinfo))
  if(any(Chrs<0)){
    stop("Chromosome index should not be negative")
  }
            
  if(any(Chrs==0)){
    Chrs = Chrs[-1] # Chrs==0 corresponds to non-genetic covariates
  }
  
  # The function "Create_Updating_Order" creates the number of "UpdatingNum" of random shuffling orders of chromosomes index of "Chrs".
  
  U.orders = Create_Updating_Order(ChromosomesOrder = Chrs, updateNum = UpdatingNum)
  
  
  testlist = c()
  KeepX = NULL
  
  # Reorder X matrix following the random shuffling orders of chromosomes index of "Chrs".
  
  for(uu in 1:UpdatingNum){
    if(length(which(Xinfo==0))>0){
      RX = X[,which(Xinfo==0),drop=F]
    }else{
      RX = NULL
    }              
              
    for(chr in U.orders[[uu]]){
      RX = cbind(RX, X[,which(Xinfo==chr),drop=F])   
    }
  
    output =  PEN(X = RX, y = y, family = family, penalty = penalty,gamma = gamma, 
      lambda.min = lambda.min, n.lambdas = n.lambdas, n.tau = n.tau, 
      lambda = lambda, tau = tau, eps = eps, max.iter = max.iter, 
      convex = convex, dfmax = dfmax, ChooseXindex = ChooseXindex, pfactor = pfactor,
      Model.selection = Model.selection)

    wkp = which(output$beta[2:(ncol(RX)+1)]!=0)
    
    # KeepX saves the chosen X names for each updating orders. 
    KeepX = c(KeepX, colnames(RX)[wkp])
  }
  
            
  # The X names that are chosen for more than 50% times will be outputted.                              
  if(length(KeepX)>0){
    KeepXTable = as.data.frame(table(KeepX))
    VoteX = KeepXTable[which(KeepXTable[,2]>UpdatingNum/2),1]
              
    val <- list(VoteX = VoteX)
  }else{
    val <- list(VoteX = NULL)}
            
  return(val)
  
}

  
  