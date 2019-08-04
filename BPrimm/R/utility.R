linkVStructures <- function(graph){
  result <- graph
  falseind <- c(NA, NA)
  
  for(i in 1:(ncol(graph) - 1)){
    atrib.to <- which(graph[, i] == 1)
    for(j in (i + 1):ncol(graph)){
      if(!(j %in% atrib.to)){
        atrib.to.j <- which(graph[, j] == 1)
        if(any(atrib.to %in% atrib.to.j)){
          result[i, j] <- 1
        }
      } else {
        atrib.to.j <- which(graph[, j] == 1)
        if(any(atrib.to %in% atrib.to.j)){
          falseind <- rbind(falseind, c(i, j))
        }
      }
    }
  }
  
  ind <- which((graph - result) != 0, arr.ind = T)
  
  output <- list()
  output[["linkedGraph"]] <- result
  output[["vstruc"]] <- ind
  output[["falsevstruc"]] <- falseind
  
  return(output)
}

mas <- function(x, y){
  z=max(abs(colSums(x*y)));
  return(z)
}

norm <- function(x){
  sqrt(sum(x ^ 2))
}

cubic <- function(x){
  return(x * x * x)
}

iden <- function(x){
  return(x)
}

quad <- function(x){
  return(x * x)
}

convert <- function(res, G){
  for(i in 1:nrow(res)){
    res.row = unlist(res[i,])
    
    #print(paste("Res.row Check:",res.row,"(",res.row[1],";",res.row[2],";",res.row[3],")"))
    
    if(!is.na(res.row[1]) & !is.na(res.row[2])){
      G[res.row[1], res.row[2]] = res.row[3]
      G[res.row[2], res.row[1]] = res.row[3]
    }
  }
  
  return(G)
}

f.trace <- function(vec){
  return(which(vec > 0))
}


traceback <- function(G.tilde, index){
  result = vector()
  
  temp = which(G.tilde[index,] > 0) #initialize the current generation
  
  if(length(temp) == 0){
    return(result)
  }
  
  result = unique(c(temp, result)) #add to traceback
  
  diff = 4
  
  while(abs(diff) > 0){
    #pull the current generation of children
    temp = unlist(apply(G.tilde[temp,,drop=FALSE], 1, f.trace))
    
    #add them to the traceback list
    oldres = result
    result = unique(c(result, temp))
    
    #stop rule is when no new children are added.
    if(length(temp) > 0){
      diff = length(oldres) - length(result)
    } else {
      diff = 0
    }
  }
  
  return(result)
  
}

calcSepSet <- function(G.tilde, i, j){
  Aij = sort(unique(c(which(G.tilde[i,] == 1), which(G.tilde[j,] == 1))))
  if(length(which(Aij %in% c(i, j)))>0){
    Aij = Aij[-which(Aij %in% c(i, j))]
  }
  
  Bij = unique(c(which(G.tilde[i,] == 1 & G.tilde[j,] == 1)))
  if(length(which(Bij %in% c(i, j))) > 0){
    Bij = Bij[-which(Bij %in% c(i, j))]
  }
  
  Cij = Bij
  
  tracebacks = list()
  
  if(length(Aij) > 0){
    for(Aij.iter in 1:length(Aij)){
      print(Aij.iter)
      index = Aij[Aij.iter]
      temp = traceback(G.tilde, index)
      tracebacks[[Aij.iter]] = temp
      
      if(Aij.iter > 1 & length(temp) > 0){
        for(Aij.iter.iter in 1:(Aij.iter-1)){
          if(any(temp %in% tracebacks[[Aij.iter.iter]])){
            Cij = unique(c(Cij,Aij[Aij.iter], Aij[Aij.iter.iter]))
          }
        }
      }
    }
  }
  
  return(Cij)
}

calcABC <- function(G.tilde, i, j){
  Aij = sort(unique(c(which(G.tilde[,i] == 1), which(G.tilde[,j] == 1))))
  if(length(which(Aij %in% c(i, j)))>0){
    Aij = Aij[-which(Aij %in% c(i, j))]
  }
  
  Bij = unique(c(which(G.tilde[,i] == 1 & G.tilde[,j] == 1)))
  if(length(which(Bij %in% c(i, j))) > 0){
    Bij = Bij[-which(Bij %in% c(i, j))]
  }
  
  Cij = Bij
  
  if(length(which(Aij %in% Bij)) > 0){
    ConCg.Bij = Aij[-which(Aij %in% Bij)]
  } else {
    ConCg.Bij = Aij
  }
  
  if(length(ConCg.Bij) > 0){
    ### Find Cij exhautively see if there's any connection between these and items in Bij
    linked = unique(c(Bij, which(rowSums(cbind(G.tilde[,Bij],0)) > 0)))
    if(length(which(linked %in% c(i,j))) > 0){
      linked = linked[-which(linked %in% c(i,j))]
    }
    
    for(con.i in 1:length(ConCg.Bij)){
      checked = linked
      wave = vector()
      if(ConCg.Bij[con.i] %in% linked){
        Cij = c(Cij, ConCg.Bij[con.i])
        checked = c(checked, ConCg.Bij[con.i])
      } else {
        sw = 0
        wave = which(G.tilde[,ConCg.Bij[con.i]] == 1)
        if(length(which(wave %in% c(i,j))) > 0){
          wave = wave[-which(wave %in% c(i,j))]
        }
        
        lcheckv = 0
        while(sw == 0 & (length(checked) + 2) < p){
          wave.connection = which(rowSums(cbind(G.tilde[,wave],0)) > 0)
          if(length(which(wave.connection %in% c(i,j))) > 0){
            wave.connection = wave.connection[-which(wave.connection %in% c(i,j))]
          }
          if(any(wave.connection %in% linked)){
            sw = 1
            Cij = c(Cij, ConCg.Bij[con.i])
            linked = c(linked, ConCg.Bij[con.i])
          }
          
          lcheck = length(checked)
          checked = unique(c(checked, ConCg.Bij[con.i], wave))
          
          if((lcheck - length(checked)) == 0){
            lcheckv = lcheckv + 1
            
            if(lcheckv > 5){
              sw = 1
            }
          }
          
          wave = wave.connection   
        }
      }
    } #End Find for Cij
  }
  
  output = list()
  output[["Aij"]] = Aij
  output[["Bij"]] = Bij
  output[["Cij"]] = Cij
  
  return(output)
}

connectedComp <- function(G,x) {
  # purpose : search for connected components of the vertices set x
  # Input   :
  #     - G : undirected adjacency matrix (edges =TRUE, gaps = FALSE)
  #     - x : set of vertices
  # Output  : the union of connected components of all elements of x (x included)
  if (!identical(dim(G)[1],dim(G)[2])) stop('G must be square matrix')
  if (!all(G==t(G))) stop('G must be symmetric')
  if (sum(diag(G))!=0) stop('diag(G) must be all zero')
  if (!is.vector(x)) stop('x must be a vector')
  p = nrow(G)
  id = seq_len(p)
  
  num_addcomp = 1L
  Ccomp = upx = x # connected components to be updated
  while (num_addcomp != 0) {
    adjmat = as.matrix(G[,upx])
    candx = id[rowSums(adjmat)!=0]
    upx = candx[!candx%in%Ccomp]
    #print(upx)
    num_addcomp = length(upx)
    Ccomp = unique(c(Ccomp,upx))
  }
  return((id %in% Ccomp))
}

repmat <- function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  return(matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T))
}

tril <- function(a,i=0){a[row(a)+i<col(a)] <- 0;a}

dist2 <- function(x, c){
  ndata = dim(as.matrix(x))[1]
  dimx = dim(as.matrix(x))[2]
  
  ncentres = dim(as.matrix(c))[1]
  dimc = dim(as.matrix(c))[2]
  
  n2 = t(matrix(rep(1,ncentres), ncol=1) %*%t(colSums(t((x^2))))) + matrix(rep(1, ndata), ncol=1) %*% t(colSums(t((c^2)))) - 2*(x%*%t(c))
  
  #Rounding errors occasionally cause negative entries in n2
  if(any(any(n2<0))){
    n2[which(n2<0)] = 0;
  }
  
  return(n2)
}

sq_dist <- function(x, c = NULL, Q = NULL){
  if(is.null(c)){
    c = x
  }
  
  D = dim(as.matrix(x))[1]
  n = dim(as.matrix(x))[2]
  
  d = dim(as.matrix(c))[1]
  m = dim(as.matrix(c))[2]
  
  if(d!=D){
    print("Error: Column lengths must agree")
    return()
  }
  
  C = matrix(0, nrow = n, ncol = m)
  
  if(is.null(Q)){
    temp = matrix(nrow = n, ncol = 100)
    for(D.iter in 1:D){
      #C = C + (repmat(b(d,:), n, 1) - repmat(a(d,:)', 1, m)).^2;
      
      C = C + (repmat(c[D.iter,,drop=FALSE], n, 1) - repmat(t(x[D.iter,,drop=FALSE]), 1, m))^2
    }        
  } else {
    if(dim(as.matrix(x)) == dim(as.matrix(Q))){
      #         C = zeros(D,1);
      # for d = 1:D
      #   C(d) = sum(sum((repmat(b(d,:), n, 1) - repmat(a(d,:)', 1, m)).^2.*Q));
      # end
      C = matrix(0,nrow=D,ncol=1)
      for(D.iter in 1:D){
        C[D.iter,] = sum(sum(((repmat(c[D.iter,,drop=FALSE], n, 1) - repmat(t(x[D.iter,,drop=FALSE]), 1, m))^2)*Q))
      }
    }
  }
  
  return(C)    
}

kernel <- function(x, xKern, theta){
  n2 = dist2(x, xKern);
  
  if(theta[1]==0){
    theta[1]=2/median(n2[which(tril(n2)>0)]);
    theta_new=theta[1];    
  }
  
  wi2 = theta[1]/2;
  
  kx = theta[2]*exp(-n2*wi2);
  bw_new=1/theta[1];
  
  result = list()
  
  result[["kx"]] = kx
  result[["bw_new"]] = bw_new
  
  return(result)   
}

gpr_multi <- function(logtheta, covfunc.gpr, x, y, xstar = NULL,  partial.derivatives = FALSE){
  n= dim(as.matrix(x))[1]
  D = dim(as.matrix(x))[2]
  m = dim(as.matrix(y))[2]
  toprint = FALSE     
  tempF1=strsplit(covfunc.gpr, ",")
  covfunc1 = tempF1[[1]][1]
  covfunc23 = paste(tempF1[[1]][2] , tempF1[[1]][3], sep=",") 
  if (length(tempF1[[1]])==3)
  {
    mo= eval(parse(text = eval(call( covfunc1, covfunc23))))
    # if (mo[1]  != dim(as.matrix(logtheta))[1])  
    # {
    #     print("Error: Number of parameters do not agree with covariance function")
    #     return (-1) 
    # }
    K = eval(call(covfunc1, covfunc23, logtheta ,x))[[1]]
  }else{
    K = eval(call(covfunc.gpr, logtheta, x)) [[1]]    
  }
  t<-try(L <- t(chol(K)), silent = TRUE)
  if(inherits(t,"try-error")){
    try(nearK <- nearPD(K, keepDiag = TRUE)$mat,silent=TRUE)
    if(nearPD(K,keepDiag = TRUE)$converged == FALSE){
      try(nearK <- nearPD(K)$mat, silent=TRUE)
      
      if(nearPD(K)$converged == FALSE){
        stop("Non-convergence of PD.")
      }
    } 
    
    t2 <-try(L <- t(chol(nearK)))   
  }
  alphat = solve(t(L),solve(L,y))    
  if (is.null(xstar))  #nargin ==4
  {
    out1 = 0.5*sum(diag(t(y)%*%alphat)) + m*sum(log(diag(L))) + 0.5*m*n*log(2*pi)
    out2=0  
    if( partial.derivatives == TRUE)
    {         
      out2 = rep(0,length(logtheta))       # set the size of the derivative vector             
      W = m*solve( t(L) ,  (solve(L, diag(n)) )  )  -  ( alphat %*% t(alphat) )  # precompute for convenience
      for (i in 1:length(out2))
      {
        if (length(tempF1[[1]])==1)
        {
          out2[i] = sum(W*eval(call(covfunc.gpr, logtheta, x, i ,  FALSE)) [[1]] )/2      
        }
        if (length(tempF1[[1]])==3)
        {
          out2[i] = sum(W*eval(call(covfunc1, covfunc23, logtheta ,  x, i  ,  FALSE))[[1]])/2  
          
        }
      }
    }# end partial derivative   
    
  }else if(!is.null(xstar))
  {             
    K = eval(call(covfunc1, covfunc23, logtheta ,  x,  xstar ,TRUE ))
    Kss= K[[1]] 
    Kstar = K[[2]]                  
    out1 = t(Kstar) %*% alphat                                      # predicted means
    out2 = 0
    if(partial.derivatives == TRUE)
    {
      if(!is.null(dim(as.matrix(x)))){  #can not say if it is null, becasue if it is an array it gives an error
        v = solve(L,Kstar)#dim v: 27,51
        out2 = Kss - t(colSums(v * v))                      
      }
    }# end partial derivative   
  }
  mout1 =  out1
  mout2 = t(out2)
  result = list(mout1, mout2) 
  return(result)
}




covSEard <- function(loghyper=NULL, x=NULL, z=NULL, test=FALSE){
  
  if(is.null(loghyper) & is.null(x) & is.null(z)){
    A = 'D+1'; 
    return(A) 
  }        
  
  n = dim(as.matrix(x))[1]
  D = dim(as.matrix(x))[2]
  ell = exp(loghyper[1:D])
  sf2 = exp(2*loghyper[D+1])
  B = NULL
  
  if(is.null(z)){
    K <<- sf2*exp(-sq_dist(diag(1/ell,nrow = length(ell), ncol=length(ell))%*%t(x))/2)  
    A = K;  
  } else if(test == TRUE){
    A = sf2*rep(1, dim(as.matrix(z))[1]);
    B = sf2*exp(-sq_dist(diag(1/ell,nrow=length(ell),ncol=length(ell))%*%t(x),diag(1/ell)%*%t(z))/2)
  } else {
    if(any(dim(K)!=n)){
      #browser()
      K <<- sf2*exp(-sq_dist(diag(1/ell,nrow=length(ell),ncol=length(ells))%*%t(x))/2)
    }
    if(z <= D){
      A = K*(sq_dist(t(x[,z])/ell[z]));  
    } else {
      A = 2*K;
      suppressWarnings(rm(K));
    }
  }
  
  result = list()
  result[["A"]] = A
  result[["B"]] = B
  return(result)
}

covSum <- function (covfuncsum, logtheta = NULL, x = NULL, z = NULL, testset.covariances = FALSE) 
{
  
  
  if(is.null(logtheta) & is.null(x) & is.null(z)){
    
    A = B = 0
    temp1 = strsplit(covfuncsum, ",")
    covfuncs1 = temp1[[1]][1]
    covfuncs2 = temp1[[1]][2]
    covarray = as.vector(c(covfuncs1, covfuncs2))
    j = as.array(c(eval(call(covfuncs1)), eval(call(covfuncs2))))
    
    A = j[1]
    
    for(i in 2:length(temp1[[1]])){
      A = paste(A, '+', j[i])
    }
    
    return(A)
  }
  
  
  A = B = 0
  D = ncol(x)
  toprint = FALSE
  temp1 = strsplit(covfuncsum, ",")
  covfuncs1 = temp1[[1]][1]
  covfuncs2 = temp1[[1]][2]
  covarray = as.vector(c(covfuncs1, covfuncs2))
  j = as.array(c(eval(parse(text=eval(call(covfuncs1)))), eval(call(covfuncs2))))
  
  v = NULL
  for (i in 1:length(j)) {
    v = cbind(v, array(rep(i, j[i]), dim = c(1, j[i])))
  }
  
  if (is.null(logtheta)) {
    A = eval(call(covfuncs1)) + eval(call(covfuncs2))
    B = 0
  } else {
    n = dim(x)[1]
    D = dim(x)[2]
    if (is.null(z)) {
      A1 = eval(call(covfuncs1, logtheta[v == 1], x))[[1]]
      A2 = eval(call(covfuncs2, logtheta[v == 2], x))[[1]]
      A = A1 + A2
      B = 0
    } else {
      if (testset.covariances == TRUE) {
        ans1 = eval(call(covfuncs1, logtheta[v == 1], 
                         x, z, testset.covariances))
        ans2 = eval(call(covfuncs2, logtheta[v == 2], 
                         x, z, testset.covariances))
        A = ans1[[1]][1] + ans2[[1]]
        B = ans1[[2]] + ans2[[2]]
      } else if (testset.covariances == FALSE) {
        i = v[z]
        j = sum(v[1:z] == i)
        f = covarray[i]
        A = eval(call(f, logtheta[v == i], x, j, testset.covariances))[[1]]
        B = 0
      }
    }
  }
  result = list(A, B)
  return(result)
}

minimize <- function (X, f, .length, covfunc, x, y) {
  #Explanation from gpml Matlab package:
  # SIG and RHO are the constants controlling the Wolfe-
  # Powell conditions. SIG is the maximum allowed absolute ratio between
  # previous and new slopes (derivatives in the search direction), thus setting
  # SIG to low (positive) values forces higher precision in the line-searches.
  # RHO is the minimum allowed fraction of the expected (from the slope at the
  # initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
  # Tuning of SIG (depending on the nature of the function to be optimized) may
  # speed up the minimization; it is probably not worth playing much with RHO.
  
  toprint = FALSE
  INT = 0.1
  EXT = 3
  MAX = 20
  RATIO = 10
  SIG = 0.1
  RHO = SIG/2
  if (is.array(.length)) {
    if (max(dim(.length)) == 2) {
      red = .length[2]
      .length = .length[1]
    }
  } else {
    red = 1
  }
  if (.length > 0) {
    S = "Linesearch"
  } else {
    S = "Function evaluation"
  }
  i.m = 0
  ls_failed = 0
  f_out = eval(call(f, X, covfunc, x, y, NULL, TRUE))
  f0 = f_out[1][[1]]
  df0 = t(f_out[2][[1]])
  fX = f0
  i.m = i.m + (.length < 0)
  s = -df0
  s = round(s * 10000)/10000
  d0 = -t(s) %*% s
  x3 = red/(1 - d0)
  mainloop = TRUE
  while (i.m < abs(.length) && mainloop) {
    i.m = i.m + (.length > 0)
    X0 = X
    F0 = f0
    dF0 = df0
    if (.length > 0) {
      M = MAX
    } else {
      M = min(MAX, -.length - i.m)
    }
    whilerun = TRUE
    while (whilerun == TRUE) {
      x2 = 0
      f2 = f0
      d2 = d0
      f3 = f0
      df3 = df0
      success = FALSE
      while (success == FALSE && M > 0) {
        M = M - 1
        i.m = i.m + (.length < 0)
        #options(show.error.messages = FALSE)
        
        f_out2 = eval(call(f, X + x3[1] * s, covfunc, 
                           x, y, NULL, TRUE))
        
        f3 = f_out2[1][[1]][[1]]
        df3 = t(f_out2[2][[1]])
        f3 = round(f3 * 10000)/10000
        df3 = round(df3 * 10000)/10000
        if (is.na(f3) || is.infinite(f3) || is.nan(f3) || 
            any(is.nan(df3) || is.na(df3) || is.infinite(df3))) {
          cat(" ")
          x3 = (x2 + x3)/2
        }
        else {
          success = TRUE
        }
        #options(show.error.messages = TRUE)
      }
      if (f3 < F0) {
        X0 = X + x3[1] * s
        F0 = f3
        dF0 = df3
      }
      d3 = t(df3) %*% s
      if (d3 > SIG * d0 || f3 > f0 + x3 * RHO * d0 || M == 
          0) {
        whilerun = FALSE
        break
      }
      x1 = x2
      f1 = f2
      d1 = d2
      x2 = x3
      f2 = f3
      d2 = d3
      A = 6 * (f1 - f2) + 3 * (d2 + d1) * (x2 - x1)
      B = 3 * (f2 - f1) - (2 * d1 + d2) * (x2 - x1)
      x3 = x1 - d1 * (x2 - x1)^2/(B + sqrt(abs(B * B - 
                                                 A * d1 * (x2 - x1))))
      if ((B * B - A * d1 * (x2 - x1) < 0)[1] || is.nan(x3) || 
          is.infinite(x3) || x3 < 0) {
        x3 = x2 * EXT
      } else if (x3 > x2 * EXT) {
        x3 = x2 * EXT
      } else if (x3 < x2 + INT * (x2 - x1)) {
        x3 = x2 + INT * (x2 - x1)
      }
      #x3 = round(x3 * 10000)/10000
    }
    while ((abs(d3) > -SIG * d0 || f3 > f0 + x3 * RHO * d0) && 
           M > 0) {
      if (d3 > 0 || f3 > f0 + x3 * RHO * d0) {
        x4 = x3
        f4 = f3
        d4 = d3
      }
      else {
        x2 = x3
        f2 = f3
        d2 = d3
      }
      if (f4 > f0) {
        x3 = x2 - (0.5 * d2 * (x4 - x2)^2)/(f4 - f2 - 
                                              d2 * (x4 - x2))
      }
      else {
        A = 6 * (f2 - f4)/(x4 - x2) + 3 * (d4 + d2)
        B = 3 * (f4 - f2) - (2 * d2 + d4) * (x4 - x2)
        x3 = x2 + (sqrt(B * B - A * d2 * (x4 - x2)^2) - 
                     B)/A
      }
      if (is.nan(x3) || is.infinite(x3)) {
        x3 = (x2 + x4)/2
      }
      x3 = max(min(x3, x4 - INT * (x4 - x2)), x2 + INT * 
                 (x4 - x2))
      f_out3 = eval(call(f, X + x3 * s, covfunc, x, y, 
                         NULL, TRUE))
      f3 = f_out3[1][[1]][[1]]
      df3 = t(f_out3[2][[1]])
      if (f3 < F0) {
        x3 = x3[[1]]
        X0 = X + x3 * s
        F0 = f3
        dF0 = df3
      }
      M = M - 1
      i.m = i.m + (.length < 0)
      d3 = t(df3) %*% s
    }
    if (abs(d3) < -SIG * d0 && f3 < f0 + x3 * RHO * d0) {
      x3 = x3[[1]]
      X = X + x3 * s
      f0 = f3
      fX = t(cbind(t(fX), f0))
      s = (((t(df3) %*% df3 - t(df0) %*% (df3))[1])/((t(df0) %*% 
                                                        (df0))[[1]]) * s) - df3
      df0 = df3
      d3 = d0
      d0 = t(df0) %*% s
      if (d0 > 0) {
        s = -df0
        d0 = -t(s) %*% s
      }
      x3 = x3 * min(RATIO, d3/(d0 - (2^(-1022))))
      ls_failed = 0
    } else {
      X = X0
      f0 = F0
      df0 = dF0
      if (ls_failed || i.m > abs(.length)) {
        mainloop = 0
        break
      }
      s = -df0
      d0 = -t(s) %*% s
      x3 = 1/(1 - d0)
      ls_failed = 1
    }
  }
  return(list(X, fX, i.m))
}









