\name{backward.MAL}
\alias{backward.MAL}
\title{ backward.MAL }
\description{
  backward filtering of the output from Multivariate Adaptive Lasso.
}
\usage{
backward.MAL(X, IA, pcut)
}
\arguments{

  \item{X}{ A numerical marix of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one covariate/marker }

  \item{IA}{ Output of function MAL, of class "MAL" }
  
  \item{pcut}{ p-value thereshold in backward filtering }
}

\value{
 A list of class "MAL" with 12 components:
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
  \item{lambda0}{ The tuning parameter corresponding to the samllest BIC }
  \item{yt}{ the trasformd response variables }
}
\author{ Wei Sun <wsun@bios.unc.edu> }
\seealso{ \code{\link{MAL}} }
\keyword{ methods }
