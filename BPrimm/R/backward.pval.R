#
# instead of backward based on the order the coefficient size
# backward based on the p-value
#

`backward.pval` <-
function(y, X, IA, pcut)
{
  if(! class(IA) %in% c("Lasso", "AL", "IAL") ){
    stop("IA must be of class Lasso, AL or IAL\n")
  }
  
  n = length(y)
  b = IA$b
  w = IA$w
  
  if(length(w)==0){
    return(IA)
  }
  
  Xn = X[,w,drop=FALSE]
  colnames(Xn) = paste(1:ncol(Xn), sep="")
	l0 = lm(y ~ Xn)
  todrop = character(0)
  
  for(i in 1:length(w)){
    pv   = summary(l0)$coef[-1,4]
    pmax = max(pv)
    wmax = which.max(pv)
    
    if(pmax < pcut){ break }
    
    todrop = c(todrop, colnames(Xn)[wmax])
    if(ncol(Xn) == 1) { break }

    Xn = Xn[,-wmax, drop=FALSE]
    l0 = lm(y ~ Xn)
  }
  
  if(length(todrop) == 0){
    IA1  = IA
  }else if (length(todrop) == length(w)){
    IA1  = IA
    IA1$b = numeric(0)
    IA1$w = numeric(0)
  }else{
    IA1 = IA
    todrop = as.numeric(todrop)
    IA1$b  = b[-todrop]
    IA1$w  = w[-todrop]
  }

  IA1
}
