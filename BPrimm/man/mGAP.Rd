\name{mGAP}
\alias{mGAP}
\title{ multivariate group-wise adaptive penalization }
\description{
  Dimension reduction regression and variable selection by Iterative Group-wise Adaptive Lasso
}
\usage{
mGAP(y, X, method=c("SIR", "cov", "spline", "power"),
 dimY=3, bs.ord=3, bs.inknots=max(1, dimY - bs.ord+1), 
 lambda=seq(1, 15, by=1), 
 tau=c(seq(.01,.08,by=.01), 0.1, 0.2), recursive=TRUE, 
 extbicgamma=1, bic.vec=TRUE,
 pMax=20, L=10, epsilon=1e-4, n.max=500, b.update.order=1,
 p=NULL, offsetrow=0, offsetcol=0, transposeX=TRUE)
}
\arguments{
  \item{y}{ A numerical vector of the quantitative trait data }

  \item{X}{ A numerical marix or a character. 
  If it is a numerical matirx, it consists of the marker genotype data, each row corresponds to one sample/individual and each column corresponds to one covariate/marker. If it is a character vector, it is the name of the genotype data file, which must be a text file with different fileds separated by "whitespace". The data read from the text file may need to be transposed, see the parameter 'transposeX' for details.}

  \item{method}{ The method to transform the response variable. Currently four transformation methods are implemented: SIR (sliced inverse regression), cov (instraslice covariance), spline (b-spline transformations of order 'bs.ord' and 'bs.inknots' knots), and power (power transformation of y). }
  
  \item{dimY}{ dimension of the transformed Y. For method "SIR", it is the number of slices minus 1; for method "cov", it is the number of slices; for method "spline" with k inner knots and order m, it is k+m-1, for method "power", it is the number of power transformations. }
  
  \item{bs.ord}{ a positive integer giving the order of the spline function. This is the number of coefficients in each piecewise polynomial segment, thus a cubic spline has order 4. The default value of bs.ord is 3. }
  
  \item{bs.inknots}{ number of inner konts in b-spline }
  
  \item{lambda}{ See explanation of \code{tau} }

  \item{tau}{ lambda and tau are two vectors of possible values for hyper-parameter \eqn{\lambda} and \eqn{\tau}, respectively. One pair of \eqn{\lambda} and \eqn{\tau} will be selected from these possible values based on \code{criterion}}

  \item{recursive}{ If TRUE, use the estimated coefficients from the last combination of lambda and tau as the initial values for the new combination of lambda and tau; otherwise use all zero's as the initial values for all lambdas and taus.}

  \item{extbicgamma}{ tuning parameter of the extended BIC }

  \item{bic.vec}{ If TRUE, the df of BIC is based on the vector forms of coefficients; otherwise based on the matrix forms. }

  \item{pMax}{ maximum number of covariates to be included }

  \item{L}{ See the explanation of \code{epsilon} }

  \item{epsilon}{ If the maximum change of the regression coefficients is smaller than \code{epsilon} for \code{L} consecutive iterations, we say the IAL is converged. }
    
  \item{n.max}{ max number of iterations. }
  
  \item{b.update.order}{ The order to update each elements of b.
   Currently b.update.order can be 1 or 2.
   1 means b are updated according to the order of the corresponding covaraites in matrix X. 2 means b are updated by a random order.}
  
  \item{p}{ if X is character of file name, i.e., genotype data need to be read from a file, p indicates the number of genetic markers }

  \item{offsetrow}{ skip the first 'offsetrow' rows in the genotype data file }

  \item{offsetcol}{ skip the first 'offsetcol' columns in the genotype data file }
  
  \item{transposeX}{ If the genotype data is read from file, whether it need to 	be transposed. By default, we assume the data file include n rows (samples/individuals) and p columns (markers/covariates), and transpose is needed. }
}

\value{
 A list with 23 components:
  \item{score_bic}{ BIC scores at each combination of lambda and tau.}
  \item{w_bic}{ which marker/covariate to keep using BIC.}
  \item{b_bic}{ regression coefficients using BIC.}
  \item{lambda_bic}{ hyperparameter lambda corresponding to the best score using BIC.}
  \item{tau_bic}{ hyperparameter tau corresponding to the best score using BIC.}
  \item{score2use_bic}{the smallest BIC score.}
  \item{score_extbic}{ Extended BIC (extBic) scores at each combination of lambda and tau.}
  \item{w_extbic}{ which marker/covariate to keep using extBIC.}
  \item{b_extbic}{ regression coefficients using extBIC.}
  \item{lambda_extbic}{ hyperparameter lambda corresponding to the best score using extBIC.}
  \item{tau_extbic}{ hyperparameter tau corresponding to the best score using extBIC.}
  \item{score2use_extbic}{the smallest extBIC score.}
  \item{score_extbicgg}{ Gaussian Graphical extended BIC (extBICgg) scores at each combination of lambda and tau.}
  \item{w_extbicgg}{ which marker/covariate to keep using extBICgg.}
  \item{b_extbicgg}{ regression coefficients using extBICgg.}
  \item{lambda_extbicgg}{ hyperparameter lambda corresponding to the best score using extBICgg.}
  \item{tau_extbicgg}{ hyperparameter tau corresponding to the best score using extBICgg.}
  \item{score2use_extbicgg}{the smallest extBICgg score.}
  \item{method}{ The method to transform the response variable }
  \item{dimY}{ the number of transfomred response varaibles }
  \item{bs.ord}{ b-spline order for spline method, default bs.ord=3 correponds to quadratic spline }
  \item{bs.inknots}{ number of inner konts in b-spline }
  \item{yt}{ the trasformd response variables }
}

\author{ Wei Sun <wsun@bios.unc.edu> }
\keyword{ methods }
