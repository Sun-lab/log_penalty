simul.ER <-
function(p,pE,n) {
    # ----------------------------------------
    # generate ER model following Colombo et al., 2012
    # p : number of vertices
    # pE : connection probability of ER model
    # n : sample size for data generation
    # ----------------------------------------
    stopifnot(pE>0,pE<1)
    A = matrix(0,p,p)
    w = which(lower.tri(A))
    A[w] = rbinom(n=length(w),size=1,prob=pE) * runif(n=length(w),min=0.1,max=1)
    Sigh = solve(diag(p) - A)
    X = rmvnorm(n=n,sigma = Sigh %*% t(Sigh))
    
    # shuffling
    id.shuffle = sample(1:p,p,replace=FALSE)
    A = A[id.shuffle,id.shuffle]
    X = X[,id.shuffle]
    return(list(A=A,X=X,G=as(t(A),"graphNEL")))
}
