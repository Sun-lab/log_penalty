\name{backward.pSDR}
\alias{backward.pSDR}
\title{ backward.pSDR }
\description{
  backward filtering of the output from penalized Sufficient Dimension Reduction.
}
\usage{
backward.pSDR(X, IA, pcut)
}
\arguments{

  \item{X}{ A numerical marix of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one covariate/marker }

  \item{IA}{ Output of function pSDR, of class "pSDR" }
  
  \item{pcut}{ p-value thereshold in backward filtering }
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
\author{ Wei Sun <wsun@bios.unc.edu> }
\seealso{ \code{\link{pSDR}} }
\keyword{ methods }
