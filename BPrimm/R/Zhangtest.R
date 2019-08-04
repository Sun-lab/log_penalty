Zhangtest  <- function(x, y, z, alpha=0.05, VERBOSE = FALSE, trans = NULL){
  # To test if x and y are independent.
  # INPUT:
  #   The number of rows of x and y is the sample size.
  #   alpha is the significance level (we suggest 1#).
  #   width contains the kernel width.
  # Output:
  #   Cri: the critical point at the p-value equal to alpha obtained by bootstrapping.
  #   Sta: the statistic Tr(K_{\ddot{X}|Z} * K_{Y|Z}).
  #   p_val: the p value obtained by bootstrapping.
  #   Cri_appr: the critical value obtained by Gamma approximation.
  #   p_apppr: the p-value obtained by Gamma approximation.
  # If Sta > Cri, the null hypothesis (x is independent from y) is rejected.
  # Copyright (c) 2010-2011  ...
  # All rights reserved.  See the file COPYING for license terms.
  
  if(any(is.null(z))){
    result = list()
    result[["crit"]] = NA 
    return(result)
  }
  
  width = 0;
  
  # Controlling parameters
  IF_unbiased = 0;
  IF_GP = 1;
  Approximate = 1;
  Bootstrap = 1; # Note: set to 0 to save time if you do not use simulation to generate the null !!!
  
  T = length(y); # the sample size
  # Num_eig = floor(T/4); # how many eigenvalues are to be calculated?
  Num_eig = T;
  T_BS = 5000; # 5000
  lambda = 1E-3; # the regularization paramter  ####Problem
  Thresh = 1E-5;
  # normalize the data
  x = x - mean(x); x = x/sd(x);
  y = y - mean(y); y = y/sd(y);
  z = apply(z, 2, scale, center = TRUE, scale = TRUE)
  
  D = ncol(z);
  logtheta_x = vector(); logtheta_y = vector();  df_x = vector(); df_y = vector();
  Cri = vector(); Sta = vector(); p_val = vector(); Cri_appr = vector(); p_appr = vector();
  
  if(width ==0){
    if(T <= 200){
      width = 1.2;
    } else if(T < 1200){
      width = 0.7;
    }
    else{
      width = 0.4;
    }
  }
  ##################################################################################
  # Developing the statistic   
  theta = 1/(width^2 * D);
  
  H =  diag(1,T) - matrix(1,T,T)/T; # for centering of the data in feature space
  
  xtemp = x;
  thetatemp = theta;
  Kx = kernel(cbind(x, z/2), cbind(x, z/2), c(theta,1))[["kx"]]
  Kx = H%*%Kx%*%H
  #browser()
  Ky = kernel(y, y, c(theta,1))[["kx"]]
  Ky = H%*%Ky%*%H 
  
  if(IF_GP){
    
    # learning the hyperparameters
    eig_Kx = eigen((Kx+Conj(t(Kx)))/2)$values[1:min(400, floor(T/4))]
    eix = eigen((Kx+Conj(t(Kx)))/2)$vectors[,1:min(400, floor(T/4))]
    
    eig_Ky = eigen((Ky+Conj(t(Ky)))/2)$values[1:min(200, floor(T/5))]
    eiy = eigen((Ky+Conj(t(Ky)))/2)$vectors[,1:min(200, floor(T/5))]
    
    #covfunc = "covSum,covSEiso,covNoise"
    covfunc = "covSum,covSEard,covNoise"
    logtheta0 = c(log(width * sqrt(D))*rep(1,D), 0, log(sqrt(0.1)))
    
    #old gpml-toolbox
    IIx = which(eig_Kx > max(eig_Kx) * Thresh)
    eig_Kx = eig_Kx[IIx]
    eix = eix[,IIx]
    IIy = which(eig_Ky > max(eig_Ky) * Thresh)
    eig_Ky = eig_Ky[IIy]
    eiy = eiy[,IIy]
    
    t = minimize(logtheta0, 'gpr_multi', -350, covfunc, z, 
                 2*sqrt(T)*eix%*%diag(sqrt(eig_Kx))/(sqrt(eig_Kx[1])))
    logtheta_x = t[[1]]
    fvals_x = t[[2]]
    iter_x = t[[3]]
    
    t2 = minimize(logtheta0, 'gpr_multi', -350, covfunc, z, 
                  2*sqrt(T)*eiy%*%diag(sqrt(eig_Ky))/(sqrt(eig_Ky[1])))
    logtheta_y = t2[[1]]
    fvals_y = t2[[2]]
    iter_y = t2[[3]]
    covfunc_z = {'covSEard'};
    
    Kz_x = eval(call(covfunc_z, logtheta_x, z))[[1]]
    Kz_y = eval(call(covfunc_z, logtheta_y, z))[[1]]
    
    # Note: in the conditional case, no need to do centering, as the regression
    # will automatically enforce that.
    
    # Kernel matrices of the errors
    temp = Kz_x + exp(2*logtheta_x[length(logtheta_x)])*diag(T)
    tt <- try(P1_x <- (diag(T) - Kz_x%*%solve(temp)), silent = TRUE)
    if(inherits(tt,"try-error")){
      P1_x <- (diag(T) - Kz_x%*%ginv(temp))
    }
    Kxz = P1_x%*%Kx%*%t(P1_x) 
    temp = Kz_y + exp(2*logtheta_y[length(logtheta_y)])*diag(T)
    tt <- try(P1_y <- (diag(T) - Kz_y%*%solve(temp)) , silent = TRUE)
    if(inherits(tt,"try-error")){
      P1_y <- (diag(T) - Kz_y%*%ginv(temp))
    }
    Kyz = P1_y%*%Ky%*%t(P1_y);
    # calculate the statistic
    Sta = sum(diag(Kxz%*%Kyz));
    
    # degrees of freedom
    df_x = sum(diag(diag(T)-P1_x));
    df_y = sum(diag(diag(T)-P1_y));
    
  } else {
    Kz = kernel(z, z, c(theta,1))
    Kz = H%*%Kz%*%H #*4 # as we will calculate Kz
    # Kernel matrices of the errors
    P1 = (diag(T) - Kz%*%solve(Kz + lambda*diag(T)));
    Kxz = P1%*%Kx%*%t(P1);
    Kyz = P1%*%Ky%*%t(P1);
    # calculate the statistic
    Sta = sum(diag(Kxz%*%Kyz))
    
    # degrees of freedom
    df = sum(diag(diag(T)-P1));
    
  }
  ##################################################################################### 
  # Approximating the null distribution
  # calculate the eigenvalues
  # Due to numerical issues, Kxz and Kyz may not be symmetric:
  eig_Kxz = eigen((Kxz+t(Kxz))/2)$values[1:Num_eig]
  eivx = eigen((Kxz+t(Kxz))/2)$vectors[,1:Num_eig]
  
  eig_Kyz = eigen((Kyz+t(Kyz))/2)$values[1:Num_eig]
  eivy = eigen((Kyz+t(Kyz))/2)$vectors[,1:Num_eig]
  
  # calculate the product of the square root of the eigvector and the eigen
  # vector
  IIx = which(eig_Kxz > max(eig_Kxz) * Thresh);
  IIy = which(eig_Kyz > max(eig_Kyz) * Thresh);
  eig_Kxz = eig_Kxz[IIx];
  eivx = eivx[,IIx];
  eig_Kyz = eig_Kyz[IIy];
  eivy = eivy[,IIy];
  
  eiv_prodx = eivx%*%diag(sqrt(eig_Kxz),nrow=length(eig_Kxz),ncol=length(eig_Kxz));
  eiv_prody = eivy%*%diag(sqrt(eig_Kyz),nrow=length(eig_Kyz),ncol=length(eig_Kyz));
  
  rm(eivx,eig_Kxz,eivy,eig_Kyz)
  # calculate their product
  Num_eigx = ncol(eiv_prodx)
  Num_eigy = ncol(eiv_prody)
  Size_u = Num_eigx * Num_eigy
  uu = matrix(0,nrow = T, ncol = Size_u)
  
  for(i in 1:Num_eigx){
    for(j in 1:Num_eigy){
      uu[,(i-1)*Num_eigy + j] = eiv_prodx[,i] * eiv_prody[,j]
    }
  }
  
  if(Size_u > T){
    uu_prod = uu %*% t(uu);
  } else {
    uu_prod = t(uu) %*% uu;
  }
  if(Bootstrap){
    eig_uu = eigen(uu_prod)$values[1:min(T,Size_u)]
    II_f = which(eig_uu > (max(eig_uu) * Thresh))
    eig_uu = eig_uu[II_f]
  }
  
  Cri=-1;
  p_val=-1;
  
  if(Bootstrap){
    # use mixture of F distributions to generate the Null dstr
    if(length(eig_uu) * T < 1e6){
      #     f_rand1 = frnd(1,T-2-df, length(eig_prod),T_BS);
      #     Null_dstr = eig_prod'/(T-1) * f_rand1;
      f_rand1 = matrix(rchisq(n = length(eig_uu)*T_BS, df = 1), nrow = length(eig_uu), ncol = T_BS)
      
      if(IF_unbiased){
        Null_dstr = T^2/(T-1-df_x)/(T-1-df_y) * eig_uu%*%f_rand1;
      } else {
        Null_dstr = eig_uu%*%f_rand1;
      }
    } else {
      # iteratively calcuate the null dstr to save memory
      Null_dstr = zeros(1,T_BS);
      Length = max(floor(1e6/T),100);
      Itmax = floor(length(eig_uu)/Length);
      for(iter in 1:Itmax){
        f_rand1 = matrix(rchisq(n = Length*T_BS, df = 1), nrow = Length, ncol = T_BS)
        if(IF_unbiased){
          Null_dstr = Null_dstr + T^2/(T-1-df_x)/(T-1-df_y) * eig_uu[(iter-1)*Length+(1:(iter*Length))] %*% f_rand1
        } else {
          Null_dstr = Null_dstr + (eig_uu[(iter-1)*Length+(1:(iter*Length))]) * f_rand1;
        }                
      }
      
    }
    
    sort_Null_dstr = sort(Null_dstr);
    Cri = sort_Null_dstr[ceiling((1-alpha)*T_BS)];
    p_val = sum(Null_dstr>Sta)/T_BS;
  }
  Cri_appr=-1;
  p_appr=-1;
  if(Approximate){
    mean_appr = sum(diag(uu_prod));
    var_appr = 2*sum(diag(uu_prod%*%uu_prod));
    k_appr = mean_appr^2/var_appr;
    theta_appr = var_appr/mean_appr;
    Cri_appr = qgamma(1-alpha, shape=k_appr, scale=theta_appr)
    p_appr = 1-pgamma(Sta, shape=k_appr, scale=theta_appr)
  }
  
  result = list()
  result[["Sta"]] = Sta
  result[["crit"]] = Cri
  result[["p_val"]] = p_val
  result[["Cri_appr"]] = Cri_appr
  result[["p_appr"]] = p_appr
  result[["eig_uu"]] = eig_uu
  return(result)
}



