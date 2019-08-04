`forward` <-
function(y, X, nper=1000, p.cut=0.05)
{
  ynew = y
  Xnew = X
  ww = pp = numeric(0)
  remains = 1:ncol(Xnew)
  cy = as.numeric(cor(ynew, Xnew))
  cy = cy^2
  r2max = max(cy)
  wy = remains[which.max(cy)]

  # use permutations to find the cutoff
  pcy = numeric(nper)
  for(per in 1:nper){
    yp  = sample(ynew, length(ynew))
    cyp = as.numeric(cor(yp, Xnew))
    cyp = cyp^2
    pcy[per] = max(cyp)
  }
  perP = length(which(pcy > r2max))/nper

  while(perP <= p.cut){
    ww = c(ww, wy)
    pp = c(pp, perP)
    
    re2drop = which(remains==wy)
    xchosen = Xnew[,re2drop]
    remains = remains[-re2drop]
    Xnew    = Xnew[,-re2drop]

    l1 = lm(ynew ~ xchosen)
    ynew = l1$resid

    cy = as.numeric(cor(ynew, Xnew))
    cy = cy^2
    r2max = max(cy)
    wy = remains[which.max(cy)]

    # use permutations to find the cutoff
    pcy = numeric(nper)
    for(per in 1:nper){
      yp  = sample(ynew, length(ynew))
      cyp = as.numeric(cor(yp, Xnew))
      cyp = cyp^2
      pcy[per] = max(cyp)
    }
    perP = length(which(pcy > r2max))/nper
  }
  
  list(ww=ww, pp=pp)
}

