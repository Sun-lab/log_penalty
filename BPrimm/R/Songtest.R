Songtest <- function(x,y,z,trans='exp',alpha=0.05, VERBOSE = FALSE){
  
  n = length(x);
  xtemp = x;
  x = z;
  z = xtemp;
  
  B=2000;
  n99=99;
  n95=95;
  n90=90;
  
  bet=1;
  sig=1;
  
  # A. Grid points for computing the Kolmogorov-Smirnov test statistic.
  
  K=10; 
  wd=t(seq(0, 1-1/K, length = K))
  
  # End of A.
  
  
  # B. Definitions of variables that are to be of repeated use.
  
  oK2=t(rep(1,K^3));
  owk1=kronecker(wd,t(rep(1,K^2)));
  owk2=kronecker(kronecker(t(rep(1,K)),wd),t(rep(1,K)));
  owk3=kronecker(t(rep(1,K^2)),wd);
  
  count = 0
  while(any(is.na(owk3)) & count < 1000){
    owk3=kronecker(t(rep(1,K^2)),wd);
    count = count+1
  }
  
  if(count == 1000){
    print(paste("owk3 has NA's, wd =", wd))  
    result = list()
    result[["crit"]] = NA
    result[["problem"]] = wd
    return(result)
  }
  
  a1=-(5^.5-1)/2;
  a2=(5^.5+1)/2;  # Two points for Rademacher variables for wild bootstrap.
  
  # End of B.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  result=vector();
  
  
  # C. Definitions of variables that are to be of repeated use.<<<<<<<<<<<<<<
  
  gr=seq(1/n, 1, by=(1/n));               # Grid for use in quantile transform of x.
  on=rep(1,n);
  on2=matrix(1, n, n);
  ond=on2-diag(n);
  ow1=on%*%owk1;
  ow2=on%*%owk2;
  ow3=on%*%owk3;
  
  # End of C.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                           
  
  bpb99=0;
  bpb95=0;
  bpb90=0;
  
  # D. Arrange the data
  
  idx = sort.int(y, index.return = T)$ix
  if(ncol(as.matrix(x)) > 1){
    x = x[idx,]
  } else {
    x = x[idx]
  }
  y = y[idx]
  z = z[idx]
  
  gy=(y%*%t(on)<=on%*%t(y));       # n by n : Indicator functions of y's.
  gz=(z%*%t(on)<=on%*%t(z));       # n by n : Indicator functions of y's.
  
  # E. Quantile transform of x <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  tx=t(rep(0,1));
  idx = sort.int(as.matrix(x)[,1], index.return = T)$ix
  tx[idx] = gr;
  
  # End of E.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  xdiff=tx%*%t(on)-on%*%t(tx);     # n by n : differences in tx 
  # (used for kernel estimation)
  
  # H. Construction of the Test Statistics
  wb99=vector();
  wb95=vector();
  wb90=vector();
  
  h1 = 0.25;
  
  h=h1*(n^(-1/5));   # Bandwidth computed from the CV
  
  xh=xdiff/h;
  nx=((15/16)*((on2-xh^2)^2)*(abs(xh)<=1)/h)*ond;
  snx=colSums(nx);
  
  tyx=colSums(gy*nx)/snx; # Rosenblatt transform
  tzx=colSums(gz*nx)/snx; # Rosenblatt transform
  
  ryx=exp((tyx%*%oK2)*ow1)*ow1-(exp(ow1)-1);     
  rzx=exp((tzx%*%oK2)*ow2)*ow2-(exp(ow2)-1);     
  
  gxe=exp((tx%*%oK2)*ow3)/n^.5;  # Exponential function of tx
  gxi=+((tx%*%oK2)<=ow3)/n^.5;     # Indicator function of tx
  
  tser=mas(gxe,rzx*ryx);      # Test statistic based on the exp. function.
  tsir=mas(gxi,rzx*ryx);      # Test statistic based on the ind. function.
  
  # End of H.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  # Wild Bootstrap Procedure. <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  wber=rep(0, B); 
  wbir=rep(0, B);
  
  ber=gxe*ryx*rzx;
  bir=gxi*ryx*rzx;
  
  for(j in 1:B){
    lot=(runif(n)<=(5^.5+1)/(2*5^.5));
    om=(a1*lot+a2*(1-lot))%*%oK2;
    wber[j]=mas(om,ber);
    wbir[j]=mas(om,bir);
  }
  
  if(VERBOSE == TRUE & sum(is.na(wber))/length(wber) > 0){
    print(paste("Proportion of Bootstram samples returning NA:", sum(is.na(wber))/length(wber)))  
    #browser()
  }
  
  if(sum(is.na(wber))/length(wber) > 0.5){
    result = list()
    result[["crit"]] = NA
    result[["problem"]] = wber
    return(result)
  }
  
  try(pwer95<-quantile(wber,1-alpha,na.rm=T));  # Bootstrap percentile
  
  try(pwir95<-quantile(wbir,1-alpha,na.rm=T));  # Bootstrap percentile
  
  # End of H.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  result = list()
  
  if(trans == 'exp'){
    try(result[["crit"]] <- pwer95)
    try(result[["Sta"]] <- tser)
  } else if(trans == 'id'){
    try(result[["crit"]] <- pwir95)
    try(result[["Sta"]] <- tsir)
  }
  return(result)
}
