\name{backward.mGAP}
\alias{backward.mGAP}
\title{ backward.mGAP }
\description{
  backward filtering of the output from Iterative Multivariate Adaptive Lasso.
}
\usage{
backward.mGAP(X, IA, pcut)
}
\arguments{

  \item{X}{ A numerical marix of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one covariate/marker }

  \item{IA}{ Output of function mGAP, of class "mGAP" }
  
  \item{pcut}{ p-value thereshold in backward filtering }
}

\value{
 A list with 10 components:
  \item{score}{ BIC or cv scores at each combination of lambda and tau.}
  \item{w}{ which marker/covariate to keep.}
  \item{b}{ regression coefficients.}
  \item{lambda}{ hyperparameter lambda corresponding to the best score.}
  \item{tau}{ hyperparameter tau corresponding to the best score.}
  \item{method}{ The method to transform the response variable }
  \item{dimY}{ the number of transfomred response varaibles }
  \item{bs.ord}{ b-spline order for spline method, default bs.ord=3 correponds to quadratic spline }
  \item{bs.inknots}{ number of inner konts in b-spline }
  \item{yt}{ the trasformd response variables }
}
\author{ Wei Sun <wsun@bios.unc.edu> }
\seealso{ \code{\link{mGAP}} }
\keyword{ methods }
