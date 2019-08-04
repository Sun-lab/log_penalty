`mofredsPC` <-
function(vGraphM, variables, VERBOSE = FALSE, 
         condtest = Zhangtest, trans = 'exp', dimRed = NULL, 
         drmethod = 'sir', alpha = 0.05){
  # Implement Stable PC Algorithm
  
  l = -1
  G = vGraphM
  if(sum(G) == 0){
    warning("Nothing was selected in Step 1: Your data may not satisfy assumptions of this method.")
    
    res = list()
    output2 = list()
    output2[["first"]] = NA
    output2[["second"]] = NA
    output2[["third"]] = "Nothing was selected in Step 1."
    
    res[[1]] = output2
    
    result = list()
    result[["res"]] = res
    result[["G"]] = G
    return(result)
  }
  S = list()
  storage = matrix(ncol = 6)
  
  connect <- as.matrix(which(G == 1, arr.ind=T))
  connect = connect[which(connect[,1] < connect[,2]),,drop=FALSE]
  
  if(nrow(connect) > 0){
    for(iterd in 1:nrow(connect)){
      i = connect[iterd,1]
      j = connect[iterd,2]
      
      if(vGraphM[i,j] == 1 & (j > i)){
        # Test Marginal Dependence:
        x = scale(variables[,i], center = FALSE)
        y = scale(variables[,j], center = FALSE)
        
        result.temp = hoeffd(x,y)$P[1,2]
        if(result.temp > 0.05){
          G[i,j] = 0
          G[j,i] = 0
          storage = rbind(storage,c(i,j,0,result.temp,NA,NA))
        }
      } # end vGraphM
    }
  }
  
  disCG = matrix(10, nrow = nrow(G), ncol=ncol(G))
  
  while(max(disCG, na.rm= T) > l){
    disCG = matrix(NA, nrow = nrow(G), ncol=ncol(G))
    
    l = l+1
    # print(paste("Filtering for conditional sizes:", l))
    G.tilde = G
    
    # Build list of connections:
    connect <- as.matrix(which(G.tilde == 1, arr.ind=T))
    connect = connect[which(connect[,1] < connect[,2]),,drop=FALSE]
    
    if(nrow(connect) > 0){
      # cl<-makeCluster(8, outfile=paste0("debug",bb,".txt"))
      # registerDoParallel(cl)
      # res = foreach(iterd = 1:nrow(connect), .packages='dr') %dopar% { 
      # tryCatch({
      res = list()
      for(iterd in 1:nrow(connect)){
        # if(iterd %in% floor(quantile(c(1:nrow(connect)), c(0.25, 0.5, 0.75)))){
        #   print(paste("Removal:", iterd,"out of",nrow(connect)))
        # }
        output = NA
        
        Sij = vector()
        
        i = connect[iterd, 1]
        j = connect[iterd, 2]
        # print(paste("i, j are",i, "and", j))
        
        # temp = calcABC(G.tilde, i, j)
        # Aij = temp$Aij
        # Bij = temp$Bij
        # Cij = temp$Cij
        
        nbrsUnion = (G.tilde[, i] | G.tilde[, j])
        nbrsUnion[i] = FALSE
        nbrsUnion[j] = FALSE
        nbrsIsect = (G.tilde[, i] & G.tilde[, j])
        nbrsIsect[i] = FALSE
        nbrsIsect[j] = FALSE
        
        nbrsdec = nbrsIsect
        if (sum(nbrsIsect) != 0) {
          Gsub = G.tilde # stability
          Gsub[,c(i,j)] = Gsub[c(j,i),] = FALSE
          nbrsdec = (connectedComp(Gsub,which(nbrsIsect)) & nbrsUnion)
        }
        
        Aij = which(nbrsUnion)
        Cij = which(nbrsdec)
        
        # step = 1
        # save(list = ls(), file=paste0("debug",iterd,".RData"))
        
        # Begin 3.2
        if(length(Cij) >= l){
          # Iteratively select a subset in Cij
          sw = 0 # switch for whether we removed the edge
          
          if(length(Cij) == l){
            Gamma = matrix(Cij,ncol = 1)
          } else {
            Gamma = combn(Cij, l) # all possible sets
          }
          
          Gamma.i = 0
          while(sw == 0){ 
            Gamma.i = Gamma.i + 1
            Gamma.temp = Gamma[,Gamma.i]
            
            if(length(which(Aij %in% Gamma.temp)) > 0){
              Kappa = Aij[-which(Aij %in% Gamma.temp)]  
            } else {
              Kappa = Aij
            }
            
            # step = 2
            # save(list=ls(), file=paste0("debug",iterd,".RData"))
            
            x = scale(variables[,i], center = FALSE)
            y = scale(variables[,j], center = FALSE)
            
            if(length(Kappa) > 1){
              if(is.null(dimRed)){
                z = apply(variables[,Kappa], 2, scale)
                
                drobj = dr(cbind(x,y) ~ z, method = drmethod, numdir = ncol(z))
                # drobj.min = min(which(dr.test(drobj, numdir=NULL)[,3] > 0.05))
                drobj.min = min(which(dr.test(drobj, numdir=length(Kappa))[,3] > 0.05))
                
                if(drobj.min == Inf){
                  drobj.min = ncol(z)
                }
                
                z.loadings = dr(cbind(x,y) ~ z, method = drmethod)$evectors[,1:drobj.min,drop=FALSE]
                z = z %*% z.loadings
              } else if(dimRed == FALSE){
                z = apply(variables[,Kappa], 2, scale, center = FALSE)
              } else if(length(Kappa) > dimRed){
                z = apply(variables[,Kappa], 2, scale, center = FALSE)
                z.loadings = dr(cbind(x,y) ~ z, method = drmethod, 
                                numdir = dimRed)$evectors[,1:dimRed,drop=FALSE]
                z = z %*% z.loadings
              } else {
                z = apply(variables[,Kappa], 2, scale, center = FALSE)
              }  
            } else if(length(Kappa) == 0){
              sw = 1
              next
            } else {
              z = scale(variables[,Kappa], center = FALSE)
            }
            
            # step = 3
            # save(list=ls(), file=paste0("debug",iterd,".RData"))
            result.temp = condtest(x,y,z, alpha, VERBOSE = VERBOSE, 
                                   trans = trans)
            if(is.na(result.temp[["crit"]])){
              # print(paste0("Error at: iterd = ", iterd))
              result = list()
              result[[1]] = result.temp
              result[[2]] = iterd
              result[['res']] = 'PROBLEM'
              return(result)
            }
            
            # step = 4
            # save(list=ls(), file=paste0("debug",iterd,".RData"))
            d <- try(if(result.temp$Sta < result.temp$crit){
              G[i,j] = 0
              G[j,i] = 0
              Sij = Kappa
              sw =  1
              storage = rbind(storage,c(i,j,0,result.temp$Sta,result.temp$crit,l))
            })
            
            if(class(d) == "try-error"){
              # changed on 1/16/2018; add paste0 function
              print(paste0("Error at: iterd = ", iterd))
              result = list()
              result[[1]] = result.temp
              result[[2]] = iterd
              result[['res']] = 'PROBLEM'
              return(result)
            }
            
            if(Gamma.i == ncol(Gamma)){
              sw = 1
            }
          }
        } # End 3.2
        
        # step = 5 
        # save(list=ls(), file=paste0("debug",iterd,".RData"))
        
        output2 = list()
        output2[["first"]] = c(i, j, G[i,j])
        output2[["second"]] = Sij
        output2[["third"]] = storage
        
        res[[iterd]] = output2
      }
      # }, error = function(e) {
      #   output2 = list()
      #   output2[["first"]] = NA #MAYBE THIS NEEDS TO BE VECTOR OF 3?
      #   output2[["second"]] = NA
      #   output2[["third"]] = NA
      #   output2[["error"]] = paste("The", iterd,
      #                              "th Iteration (with i=",i," and j=",j,") Caused Error:",e,".")
      #   print("Error in step2Func, check output2.")
      #   return(output2)
      # }) #end tryCatch
      # } #end parallel
      # stopCluster(cl)
    }
    
    # CHECK RES HERE
    
    res2 = t(sapply(res, "[[", "first"))
    G = convert(res2, G)
    
    # REMEMBER STILL NEED TO FIGURE OUT A WAY TO STORE SIJ
    
    disCG = vector()
    for(i in 1:(ncol(G) - 1)){
      for(j in (i+1):ncol(G)){
        if(G[i,j] != 0){
          # temp = calcABC(G, i, j)
          # disCG = c(disCG, length(temp$Cij))
          nbrsUnion = (G.tilde[, i] | G.tilde[, j])
          nbrsUnion[i] = FALSE
          nbrsUnion[j] = FALSE
          nbrsIsect = (G.tilde[, i] & G.tilde[, j])
          nbrsIsect[i] = FALSE
          nbrsIsect[j] = FALSE
          
          nbrsdec = nbrsIsect
          if (sum(nbrsIsect) != 0) {
            Gsub = G.tilde # stability
            Gsub[,c(i,j)] = Gsub[c(j,i),] = FALSE
            nbrsdec = (connectedComp(Gsub,which(nbrsIsect)) & nbrsUnion)
          }
          
          Cij = which(nbrsdec)
          
          disCG = c(disCG, length(Cij))
        }
      }
    }
  }
  
  result = list()
  result[["res"]] = res
  result[["G"]] = G
  return(result)
}  
