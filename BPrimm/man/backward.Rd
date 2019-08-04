\name{backward}
\alias{backward}
\title{ backward }
\description{
  backward filtering of the output from Lasso, adpative Lasso, and Iterative Adaptive Lasso.
}
\usage{
backward(y, X, IA, pcut)
}
\arguments{
  \item{y}{ A numerical vector of the quantitative trait data }

  \item{X}{ A numerical marix of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one covariate/marker }

  \item{IA}{ Output of function IAL, of class "eNet", "AL", or "IAL"}
  
  \item{pcut}{ p-value thereshold in backward filtering }
}

\value{
 A list of 'IAL' with 7 components:
 \item{score}{ BIC or cv scores at each combination of delta and tau.}
 \item{b0}{ intercept.}
 \item{w}{ which marker/covariate to keep.}
 \item{b}{ regression coefficients.}
 \item{delta}{ hyperparameter delta corresponding to the best score.}
 \item{tau}{ hyperparameter tau corresponding to the best score.}
 \item{score2use}{ the best score among all combinitions of delta and tau.}
}

\references{
Wei Sun, Joseph Ibrahim, and Fei Zou, Variable Selection by Bayesian Adaptive Lasso and Iterative Adaptive Lasso, with Application for Genome-wide Multiple Loci Mapping.
}
\author{ Wei Sun <wsun@bios.unc.edu> }
\seealso{ \code{\link{IAL}}, code{\link{AL}}, code{\link{eNet}} }
\keyword{ methods }
