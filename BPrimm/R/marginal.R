`marginal` <-
function(y, X, chrs, nper=1000, pcut=0.05)
{
  if(!is.null(chrs) && length(chrs) != ncol(X)){
    stop("length of chrs does not match the columns of X\n")
  }
  
  cy = as.numeric(cor(y, X))
  cy = cy^2
  
  # use permutations to find the cutoff
  pcy = numeric(nper)
  for(per in 1:nper){
    yp  = sample(y, length(y))
    cyp = as.numeric(cor(yp, X))
    cyp = cyp^2
    pcy[per] = max(cyp)
  }

  if(is.null(chrs)){
    r2Cut = quantile(pcy, probs=1-pcut)
    wkps = which(cy > r2Cut)
    nomps = pkps = numeric(length(wkps))
    
    for(j in 1:length(wkps)){
      lj = lm(y~X[,wkps[j]])
      fj = summary(lj)$fstat
      nomps[j] = pf(fj[1], fj[2], fj[3], lower.tail = FALSE)
      pkps[j]  = length(which(pcy > cy[wkps[j]]))/nper
    }
    
  }else{
    wkps = pkps = nomps = numeric(0)
    for(unchr in unique(chrs)){
      wchr = which(chrs==unchr)
      wkp1 = wchr[which.max(cy[wchr])]
      pkp1 = length(which(pcy > cy[wkp1]))/nper
      wkps = c(wkps, wkp1)
      pkps = c(pkps, pkp1)
      lj   = lm(y~X[,wkp1])
      fj   = summary(lj)$fstat
      nomps = c(nomps, pf(fj[1], fj[2], fj[3], lower.tail = FALSE))
    }
  }
  
  list(ww=wkps, pp=pkps, nomp=nomps)
}
