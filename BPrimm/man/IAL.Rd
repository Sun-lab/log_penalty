\name{IAL}
\alias{IAL}
\title{ Iterative Adaptive Lasso }
\description{
  An iterative addaptive Lasso
}
\usage{
IAL(y, X, b=NULL, delta = seq(0, 5, by=0.5), 
 tau = c(1e-3, 5e-3, 0.01, 0.02, 0.05, 0.1, 0.2), 
 criterion="BIC", cv=10, times.cv=10, BIC.gamma=0, pMax=20, 
 L=10, epsilon=1e-4, n.max=1000, b.update.order=1,
 p=NULL, offsetrow=0, offsetcol=0, transposeX=FALSE)
}
\arguments{
  \item{y}{ A numerical vector of the quantitative trait data }

  \item{X}{ A numerical marix or a character. 
  If it is a numerical matirx, it consists of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one covariate/marker. 
  If it is a character vector, it is the name of the genotype data file. 
  The genotype data cannot be read from a file if tau and delta need to 
  be selected by cross validation. }

  \item{b}{ Initial estimates of regression coefficients, if it is NULL, 
  all the initial values are set as 0.}
  
  \item{delta}{ See explanation of \code{tau} }

  \item{tau}{ delta and tau are two vectors of possible values for hyper-parameter \eqn{\delta} and \eqn{\tau}, respectively. One pair of \eqn{\delta} and \eqn{\tau} will be selected from these possible values based on \code{criterion}}

  \item{criterion}{ which criterion to be used to select delta and tau, the default is BIC. Another criterion is cross validation: cv. }
  
  \item{cv}{ folds of cross-validation }
  
  \item{times.cv}{ how many times of "cv" folds of cross-validation }

  \item{BIC.gamma}{ the coefficient of the term log(number of models with k variables) in the modified BIC, see Chen and Chen (2008). By default, we set it as 0, so that it is traditional BIC, and we further select variables by backward filtering. }

  \item{pMax}{ maximum number of covariates to be included }

  \item{L}{ See the explanation of \code{epsilon} }

  \item{epsilon}{ If the maximum change of the regression coefficients is smaller than \code{epsilon} for \code{L} consecutive iterations, we say the IAL is converged}
    
  \item{n.max}{ max number of iterations. }
  
  \item{b.update.order}{ The order to update each elements of b.
   Currently b.update.order can be 1, 2, or 3.
   1 means b are updated according to the order of the corresponding covaraites in matrix X. 
   2 means b are updated by a random order. 
   3 means b are updated based on their absolute values in the previous iteration (from the largest to the smallest). The initial order is based on the marginal regression coefficients. }
   
  \item{p}{ if X is character, i.e., genotype data need to be read from a file, p indicate the number of genetic markers }

  \item{offsetrow}{ skip the first 'offsetrow' rows in the genotype data file }

  \item{offsetcol}{ skip the first 'offsetcol' columns in the genotype data file }
  
  \item{transposeX}{ If the genotype data is read from file, whether it need to be transposed. }
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
\seealso{ \code{\link{Bayes.t}}, \code{\link{Bayes.Lasso}} }
\keyword{ methods }
