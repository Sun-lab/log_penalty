detMatrix <- function(X) {
  if(!is.matrix(X)){
    stop("X must be a matrix\n")
  }
  
  if(nrow(X) != ncol(X)){
    stop("X must be a square matrix\n")
  }
  
  det = 0.0;
  Z = .C("detR", as.double(X), as.integer(nrow(X)), det=as.double(det),
         PACKAGE="BPrimm")
  Z[["det"]]
}
