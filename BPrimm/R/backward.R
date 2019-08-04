`backward` <-
function(y, X, IA, pcut)
{
  if(! class(IA) %in% c("eNet", "AL", "IAL") ){
    stop("IA must be of class Lasso, AL or IAL\n")
  }
  
  n = length(y)
  b = IA$b
  w = IA$w
  
  if(length(w)==0){
    return(IA)
  }
	
  od = order(abs(b), decreasing=TRUE)
  
  wuse = 0
  l0 = lm(y~X[,w])
  
  for(i in length(od):1){
  
    l1 = l0
    
    if(i == 1){
      l0 = lm(y~1)
    }else{
      ww = w[od[1:(i-1)]]
      l0 = lm(y~X[,ww])
    }
    
    pval = anova(l0, l1)$Pr[2]
    
    if(is.na(pval)){
      next
    }
    
    if(pval < pcut){
      wuse = i
      break
    }
  }
  
  IA1  = IA

  if(wuse > 0){
    wmin  = od[1:wuse]
    IA1$b = b[wmin]
    IA1$w = w[wmin]
  }else{
    IA1$b = numeric(0)
    IA1$w = numeric(0)
  }
    
  IA1
}
