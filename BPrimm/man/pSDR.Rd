\name{pSDR}
\alias{pSDR}
\title{
  penalized Sufficient Dimension Reduction
}
\description{
  penalized Sufficient Dimension Reduction, using ridge or paritial least squares to circumvent the problem of singular covariance matrix
}
\usage{
pSDR(y, X, method = c("SIR", "cov", "spline", "power"), 
  covMat = c("ridge", "pls"), dimY = 3, d = 1, alpha = 0.95, 
  pMax = 20, bs.ord = 3, bs.inknots = max(1, dimY - bs.ord + 1), 
  taus=c(seq(0.1, 1, by=0.1), 1.2, seq(1.6, 3, by=0.4)), u.max = 15, 
  threshold = 2)
}
\arguments{
  \item{y}{ A numerical vector of the quantitative trait data }

  \item{X}{ A numerical marix. It consists of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one covariate/marker. } 
  
  \item{method}{ The method to transform the response variable }

  \item{covMat}{ The approach to avoid singularity of the covaraince matrix: either "ridge" or "pls"}
  
  \item{dimY}{ dimension of the transformed Y. For method "SIR", it is the number of slices minus 1; for method "cov", it is the number of slices; }

  \item{d}{ The dimension of the central sub space }

  \item{alpha}{ The elastic-net mixing parameter used in function glmnet, with 0 < alpha <= 1. The penalty is defined as (1-alpha)*||b|| + alpha*|b|. Alpha=1 is the lasso penalty. By default, we set alpha=0.95 for high-dimensino low sample size setting. }

  \item{pMax}{ maximum number of covariates to be included }

  \item{bs.ord}{ A positive integer giving the order of the spline function. This is the number of coefficients in each piecewise polynomial segment, thus a quadratic spline has order 3. Defaults to 3. }

  \item{bs.inknots}{ Number of inner konts in b-spline }

  \item{taus}{ The taus are the penalization parameters to be used for ridge correction of the covariance matrix: cov(X) + tau*I  }
  
  \item{u.max}{ partial least squares regression coefficient can be write as 
     bPLS = R %*% (t(R) %*% cov(X) %*% R)^{-1} R %*% cov(X,y). 
     R = (cov(X,y), cov(X)%*%cov(X,y), ..., cov(X)^{q-1}%*%cov(X,y)), 
     u.max = q, the maximum order we consider in the expansion.
     The actual number of u used is decided by threhold, such that the ratio of 
     the adjacent eigen-values of the matrix R%*%t(R) is larger than threshold.
   }
   
  \item{threshold}{see the explainaiton for u.max}
  
}
\value{
 A list with 12 components:
  \item{w}{ which marker/covariate to keep.}
  \item{b}{ regression coefficients.}  
  \item{method}{ The method to transform the response variable }
  \item{covMat}{ The approach to avoid singularity of the covaraince matrix }
  \item{dimY}{ the number of transfomred response varaibles }
  \item{bs.ord}{ b-spline order for spline method, default bs.ord=3 correponds to quadratic spline }
  \item{bs.inknots}{ number of inner konts in b-spline }
  \item{BIC }{ BIC of the chosen model }
  \item{tau }{ tuning parameter selected for ridge penalty of the covariance matrix }  
  \item{u }{ tuning parameter selected for pls approach to avoid the singular  covariance matrix }
  \item{d}{ The dimension of the central sub space }
  \item{yt}{ the trasformd response variables }
}

\author{ Wei Sun <wsun@bios.unc.edu> and Lexin Li <li@stat.ncsu.edu>}
\keyword{ methods }
