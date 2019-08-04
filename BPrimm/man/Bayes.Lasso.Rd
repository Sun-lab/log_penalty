\name{Bayes.Lasso}
\alias{Bayes.Lasso}

\title{ Bayesian shrinkage method, with Lasso penalty }
\description{
  This Bayesian shrinkage method can be viewed as a Bayesian version
  of the Lasso.
}
\usage{
Bayes.Lasso(y, X, b=NULL, n.burn = 2000, n.thin = 20, n.iter = 100, s = 0.01,
 r = 0.01, b.update.order = 1, method=1)
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
  \item{s}{ hyperparameter, see details}
  \item{r}{ hyperparameter, see details }
  \item{b.update.order}{ The order to update each elements of b.
   Currently b.update.order can be 1, 2, or 3.
   1 means b is ordered by their order in matrix X.
   2 means for each update, b is randomly ordered.
   3 means b is ordred by its absolute values,
   the initial order is based on marginal regression coefficients. }
  \item{method}{only two values are allowed, 1 or 2, indicating
  two methods, see details}
}
\details{
This is a herarchical model.
If method equals to 1, the model is as follows.
  First, each ceofficeint \eqn{b_j} follows a normal distirbution with
  mean 0 and variance \eqn{\sigma_j^2}.
  Second, \eqn{\sigma_j^2} follows an exponential distribution with
  parameter \eqn{\kappa^2/2}.
  Third, \eqn{\kappa^2/2} follows a Gamma distribution with
  hyperparameter parameter s and r.

If method equals to 2, the model is the same as model 1, except
each ceofficeint \eqn{b_j} follows a normal distirbution with
  mean 0 and variance \eqn{\sigma_0^2\sigma_j^2}, where \eqn{\sigma_0^2} is
  the variance of residual error.
  
Method 1 and 2 corresponds to method I and II in Yi and Xu (2008).
}

\value{
A list with two elements:
 \item{b}{ a matrix of the samples of coefficients b.
 Each column is a vector of b from one sample.}
 \item{kappa}{ a vector of the samples of kappa}
}
\references{
Yi, N. and Xu, S., 2008. Bayesian LASSO for Quantitative Trait Loci
Mapping. Genetics, 179:1045-1055.

Wei Sun, Joseph Ibrahim, Fei Zou, Genome-wide Multiple Loci Mapping
in Experimental Crosses Using Dense Genetic Markers
}

\author{ Wei Sun <wsun@bios.unc.edu> }
\seealso{ \code{\link{Bayes.t}}, \code{\link{Bayes.AL}} }
\keyword{ methods }
