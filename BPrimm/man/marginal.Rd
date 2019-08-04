\name{marginal}
\alias{marginal}
\title{ marginal regression for muliple loci mapping }
\description{
  a marginal regression model for muliple loci mapping
}
\usage{
marginal(y, X, chrs, nper = 1000, pcut=0.05)
}
\arguments{
  \item{y}{ A numerical vector of the quantitative trait data }
  \item{X}{ A numerical marix of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one marker }
  \item{chrs}{ A vector of the same length as one column of X,
  indicating the chromosome of each marker }
  \item{nper}{ The number of permutations in each forward selection step }
  \item{pcut}{ p-value cutoff }
}
\details{
  With the assumption of no missing data,
  r-square is used as the statistics to evaluate permuation p-value.
  In each permutation, we record the highest r-square between any
  genotype profile and the trait. The permutation p-value for an
  observed r-square is calcualted as the proportion of permutaions
  where we observed at least one r-square bigger than the observed one.
  In each chromosome, only the most significantly linked marker is
  returned.
}
\value{
 a list with two components
  \item{ww }{ the selected marker }
  \item{pp }{ permutation p-value of each selected marker}
}
\author{ Wei Sun <wsun@bios.unc.edu> }
\keyword{ methods }
