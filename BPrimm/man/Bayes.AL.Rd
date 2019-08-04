\name{Bayes.AL}
\alias{Bayes.AL}
\title{ Bayesian shrinkage method, with adpative Lasso (L1) penalty }
\description{
  This Bayesian shrinkage method can be viewed as a Bayesian version
  of the multiple step adpative Lasso.
}
\usage{
Bayes.AL(y, X, b=NULL, n.burn=2000, n.thin=20, n.iter=100,
  delta=NULL, tau=NULL, b.update.order=1)
}
\arguments{
  \item{y}{ A numerical vector of the quantitative trait data }
  \item{X}{ A numerical marix of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one marker }
  \item{b}{ Initial estimates of regression coefficients, if it is NULL, 
  all the initial values are set as 0.}
  \item{n.burn}{ Number of burn in iterations }
  \item{n.thin}{ Number of iterations between two sampled b }
  \item{n.iter}{ Number of coefficients b to sample }
  \item{delta}{ hyperparameter, see details }
  \item{tau}{ hyperparameter, see details }
  \item{b.update.order}{ The order to update each elements of b.
   Currently b.update.order can be 1, 2, or 3.
   1 means b is ordered by their order in matrix X.
   2 means for each update, b is randomly ordered.
   3 means b is ordred by its absolute values,
   the initial order is based on marginal regression coefficients. }
}
\details{
This is a herarchical model.
  First, each ceofficeint \eqn{b_j} follows a double exponential
  distirbution with parameter \eqn{\kappa_j}.
  Second, \eqn{\kappa_j^2} follows an invers-Chi-square distribution with
  parameter \eqn{\delta} and \eqn{\tau}.
  \eqn{\delta} and \eqn{\tau} can be specificed by the user. 
  Otherwise (by default, they are NULL) they
  are estimated by imposing a prior:
  \eqn{\pi(\delta, \tau) \sim 1/\tau}.
}
\value{
 A list with three elements:
 \item{b}{ a matrix of the samples of coefficients b.
 Each column is a vector of b from one sample.}
 \item{delta}{ this is the input delta}
 \item{tau}{ this is the input tau}
}
\references{
Wei Sun, Joseph Ibrahim, Fei Zou, Genome-wide Multiple Loci Mapping
in Experimental Crosses Using Dense Genetic Markers
}
\author{ Wei Sun <wsun@bios.unc.edu> }
\seealso{ \code{\link{Bayes.t}}, \code{\link{Bayes.Lasso}} }
\keyword{ methods }
