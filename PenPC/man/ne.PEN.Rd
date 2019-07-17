\name{ne.PEN}
\alias{ne.PEN}
\title{Neighborhood Selection using penalized regressions with log penalty
}
\description{Fit penalized regressions with log penalty to estimate neighborhoods for a set of vertices
}
\usage{
ne.PEN(dat, nlambda, ntau, V, order = FALSE, verbose = FALSE, Model.selection="ExtendedBIC")
}
\arguments{
  \item{dat}{ nxp numeric data matrix with n (number of sample) and p (number of variables)
}
  \item{nlambda}{
  	 a numeric constant, the number of candidate lambda as a tuning parameter of Log penalty
}
  \item{ntau}{
     a numeric constant, the number of candidate tau as a tuning parameter of Log penalty}
  \item{V}{
     a numeric vector, the set of vertices which neighbors are estimated
     }
  \item{order}{
      If \code{TRUE}, covariate order are in the decreasing order of absolute marginal correlation}
  \item{verbose}{
  	  If \code{TRUE}, detailed output is provided
 }
  \item{Model.selection}{
  	  the model selection criterion, allowed values are "BIC" or "ExtendedBIC".
 }
}
\details{
	The tuning parameters, lambda and tau, for a penalized regressions are estimated using extended BIC.
}
\value{ a px\code{length(V)} matrix. Each column corresponds to the estimated coefficients from the penalized regression 
}
\references{
 Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.
 
  Chen, J. and Chen, Z. (2008) Extended Bayesian information criteria for model selection with large model spaces. Biometrika., 95: 759-771.
  }
 \author{
Min Jin Ha (\email{minjin.ha@gmail.com})}
\examples{


# simulate a DAG following ER model
p=100
simul = simul.ER(p=p,pE=0.02,n=30)
dat =  simul$X
# Estimate neighbors of all vertices
coeff = ne.PEN(dat=dat,nlambda=100,ntau=10,V=1:p,order=TRUE,verbose=TRUE)
# Fit Partial correlation graph with OR rule
Gmat = matrix(0,p,p)
Gmat[coeff!=0 | t(coeff)!=0] =1
PCgraph = as(Gmat,"graphNEL")

}
