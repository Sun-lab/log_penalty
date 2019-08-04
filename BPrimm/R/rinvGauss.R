rinvGauss <- function (n, mu, lambda){

  x = numeric(n)

  Z = .C("rinvGauss", x=as.double(x), as.integer(n),
     as.double(mu), as.double(lambda), PACKAGE="BPrimm")

  Z[["x"]]
}
