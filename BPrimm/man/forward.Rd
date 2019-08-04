\name{forward}
\alias{forward}
\title{ forward regression }
\description{
  A forward regression model for muliple loci mapping
}
\usage{
forward(y, X, nper = 1000, p.cut = 0.05)
}

\arguments{
  \item{y}{ A numerical vector of the quantitative trait data }
  \item{X}{ A numerical marix of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one marker }
  \item{nper}{ The number of permutations in each forward selection step }
  \item{p.cut}{ Permutation p-value cutoff, the forward selection
  is stopped when the minimum p-value from all the remaining marakers
  are bigger than p.cut }
}
\details{
  Permutations p-value is used in each forward step to decide the
  significance level. In the first step, the trait y is permuted.
  In the following steps, the residuals of current model is permuted.
  This is actually the residual empirical threshold (RET) in
  Doerge and Churchill (1996). With the assumption of no missing data,
  r-square is used as the statistics to evaluate permuation p-value.
}
\value{
 a list with two components
  \item{ww }{ the selected marker }
  \item{pp }{ permutation p-value of each selected marker}
}
\references{
Doerge, R. and Churchill, G., 1996,
Permutation tests for multiple loci affecting a quantitative character,
Genetics, 142, 285-294.
}
\author{ Wei Sun <wsun@bios.unc.edu> }
\keyword{ methods }
