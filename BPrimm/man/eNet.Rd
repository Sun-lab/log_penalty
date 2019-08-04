\name{eNet}
\alias{eNet}
\title{ eNet }
\description{
  eLastic Net. 
}
\usage{
eNet(y, X, lambda=NULL, pMax=20, alpha=0.95, backward.pval=NULL)
}
\arguments{
  \item{y}{ A numerical vector of the quantitative trait data }

  \item{X}{ A numerical marix of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one covariate/marker }

  \item{lambda}{ the tuning praramter of the eNet, the coefficient of the penalziation term }
    
  \item{pMax}{ maximum number of covariates to be included }
  
  \item{alpha}{ The elastic-net mixing parameter used in function glmnet, with 0 < alpha <= 1. The penalty is defined as (1-alpha)*||b|| + alpha*|b|. Alpha=1 is the eNet penalty. By default, we set alpha=0.95 for high-dimensino low sample size setting. }
  
  \item{backward.pval}{ if backward.pval is not null, the backward function is run within eNet. }

}

\value{
 A list of 'eNet' with 4 components:
 \item{BIC}{ best BIC for each gamma.}
 \item{w}{ which marker/covariate to keep.}
 \item{b}{ regression coefficients.}
 \item{lambda}{ tuning parameter lambda.}
 \item{alpha}{ tuning parameter alpha.}
}
\references{ J Friedman, T Hastie, R Tibshirani 2008. Regularization paths for generalized linear models via coordinate descent
}
\author{ Wei Sun <wsun@bios.unc.edu> }
\note{ the eNet calculation within this function is carried out by function glmnet in library glmnet. 
}
\seealso{ \code{\link{AL}} }
\keyword{ methods }
