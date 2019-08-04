\name{gentp}
\alias{gentp}
\title{ Generate the grid of tuning parameters for mGAP }
\description{
  Generate the grid of tuning parameters for mGAP
}
\usage{
gentp(y, X, method=c("SIR", "cov", "spline", "power"),
      dimY=3, bs.ord=3, bs.inknots=max(1, dimY - bs.ord+1), 
      ntau = 5, nlambda = 50, lambda_min = 1e-6)
}
\arguments{
  \item{y}{ A numerical vector of the quantitative trait data }

  \item{X}{ A numerical marix or a character. 
  If it is a numerical matirx, it consists of the marker genotype data, each row corresponds to one sample/individual and each column corresponds to one covariate/marker. If it is a character vector, it is the name of the genotype data file, which must be a text file with different fileds separated by "whitespace". The data read from the text file may need to be transposed, see the parameter 'transposeX' for details.}

  \item{method}{ The method to transform the response variable. Currently four transformation methods are implemented: SIR (sliced inverse regression), cov (instraslice covariance), spline (b-spline transformations of order 'bs.ord' and 'bs.inknots' knots), and power (power transformation of y). }
  
  \item{dimY}{ dimension of the transformed Y. For method "SIR", it is the number of slices minus 1; for method "cov", it is the number of slices; for method "spline" with k inner knots and order m, it is k+m-1, for method "power", it is the number of power transformations. }
  
  \item{bs.ord}{ a positive integer giving the order of the spline function. This is the number of coefficients in each piecewise polynomial segment, thus a cubic spline has order 4. The default value of bs.ord is 3. }
  
  \item{bs.inknots}{ number of inner konts in b-spline }

  \item{ntau}{ number of tau values - default is 5}

  \item{nlambda}{ number of lambda values - default is 50}

  \item{lambda_min}{ minimum of lambda}


}

\value{
A list of the following components:
 \item{lambda}{ A sequence of lambda values.}
 \item{tau}{ A sequence of lambda values.}
}

\author{ Wei Sun <wsun@bios.unc.edu> }
\keyword{ methods }
