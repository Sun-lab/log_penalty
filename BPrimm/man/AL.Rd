\name{AL}
\alias{AL}
\title{ Adaptive Lasso }
\description{
  Adaptive Lasso
}
\usage{
AL(y, X, b=NULL, lambda=NULL, gamma=seq(0.5,2,by=0.5), pMax=20, alpha=0.95)
}
\arguments{
  \item{y}{ A numerical vector of the quantitative trait data }

  \item{X}{ A numerical marix of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one covariate/marker }

  \item{b}{ Initial estimates of regression coefficients, if it is NULL, 
  all the initial values are set by marginal regression.}

  \item{lambda}{ One of the tuning praramter of the adaptive Lasso, the coefficient of the penalziation term }
  
  \item{gamma}{ One of the tuning parameter of the adaptive Lasso, the power of the penalization term. gamma=1 is a special case of adaptive Lasso that is closely related with nonnegative garotte}
  
  \item{pMax}{ maximum number of covariates to be included }
  
  \item{alpha}{ The elasticnet mixing parameter used in function glmnet, with 0 < alpha <= 1. The penalty is defined as (1-alpha)*||b|| + alpha*|b|. Alpha=1 is the lasso penalty. By default, alpha=0.95. }
}

\value{
 A list of 'AL' with 5 components:
 \item{BIC}{ best BIC for each gamma.}
 \item{w}{ which marker/covariate to keep.}
 \item{b}{ regression coefficients.}
 \item{lambda}{ tuning parameter lambda.}
 \item{gamma}{ tuning parameter gamma.}
}
\references{ Zou, H., 2006. The adaptive lasso and its oracle properties. Journal of the American Statistical Association, 101:1418-1429. 
}
\author{ Wei Sun <wsun@bios.unc.edu> }
\note{ the Lasso calculation within this function is carried out by function glmnet in library glmnet. 
}
\seealso{ \code{\link{IAL}} }
\keyword{ methods }
