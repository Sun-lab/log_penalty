\name{dinvGauss}
\alias{rinvGauss}
\alias{dinvGauss}

\title{ Inverse Gaussian Distribution}
\description{
  Density function and random generation for the inverse Gaussian
  distribution.
}
\usage{
dinvGauss(x, mu, lambda)
rinvGauss(n, mu, lambda)
}

\arguments{
  \item{x}{ a vector of the values of the random variable }
  \item{n}{ the number of random numbers to generate }
  \item{mu}{ mean parameter }
  \item{lambda}{ shape parameter }
}
\value{
  dinvGauss gives density, rinvGauss genrates random deviates.
}
\references{
Chhikara, R.S. and Folks, J.L., 1989.
The inverse gaussian distribution: theory, methodology, and
  applications, Marcel Dekker, Inc., New York, NY, USA.
 }
\author{ Wei Sun }
\keyword{ utilities }
