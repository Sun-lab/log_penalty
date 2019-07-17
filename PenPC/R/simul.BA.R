simul.BA <-
function(p,e,n) {
    # ----------------------------------------
    # generate BA model 
    # p : number of vertices
    # e : number of edges 
    # n : sample size for data generation
    # ----------------------------------------
   A=as(igraph.to.graphNEL(barabasi.game(n=p,m=e)),"matrix")
   w = which(A!=0)
     A[w] =  runif(n=length(w),min=0.1,max=1)
     Sigh = solve(diag(p) - A)
     X = rmvnorm(n=n,sigma = Sigh %*% t(Sigh))
    
     # shuffling
     id.shuffle = sample(1:p,p,replace=FALSE)
     A = A[id.shuffle,id.shuffle]
     X = X[,id.shuffle]
     return(list(A=A,X=X,G=as(t(A),"graphNEL")))
}
