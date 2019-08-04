dinvGauss <- function(x, mu, lambda){

  n = length(x)
  y = numeric(n)
  
  Z = .C("dinvGauss", y=as.double(y), as.double(x),
     as.integer(n), as.double(mu), as.double(lambda),
     PACKAGE="BPrimm")
  Z[["y"]]
}
