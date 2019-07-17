\name{PEN}
\alias{PEN}
\title{ Penalized linear regression and logistic regression.}
\description{ Penalized linear regression and logistic regression with SCAD, MCP, Log or SICA 
penalties.}
\usage{
  PEN(X, y, family=c("gaussian","binomial"),
  penalty=c("SCAD","MCP","LOG","SICA"), gamma=3, lambda.min=NULL,
  n.lambdas=100, n.tau=6, lambda=NULL, tau=NULL, eps=.001,
  max.iter=1000, convex=TRUE, dfmax=NULL, ChooseXindex=NULL,
  pfactor=0.1, Model.selection="ExtendedBIC", restriction="none")
}
\arguments{
  \item{X}{ A numerical design matrix.}
  \item{y}{ Continuous response for linear regression or binary (0, 1) response for logistic
  regression.}
  \item{family}{ The specified distribution function is gaussian or binomial.}
  \item{penalty}{ The applied penalty functions including SCAD, MCP, LOG and SICA.}
  \item{gamma}{ The tuning parameter of MCP or SCAD penalty.}
  \item{lambda.min}{ The minimum value of lambda.}
  \item{n.lambdas}{ The number of values of lambdas.}
  \item{n.tau}{ The number of values of taus for Log and SICA.}
  \item{lambda}{ A user-specified sequence of lambda values.  By default, a sequence of values of 
  length \code{nlambda} is computed, equally spaced on the log scale.}
  \item{tau}{ A user-specified sequence of tau values.}
  \item{eps}{ Convergence threshhold. The algorithm iterates until the relative change in any 
  coefficient is less than \code{eps}.  Default is \code{.001}.}
  \item{max.iter}{ Maximum number of iterations.  Default is 1000.}
  \item{convex}{ Calculate index for which objective function ceases to be locally convex?  
  Default is TRUE, only for SCAD and MCP.}
  \item{dfmax}{ Upper bound for the number of nonzero coefficients. Default is no upper bound.}
  \item{ChooseXindex}{ The user-specified index vector with length equal to number of covariates. 
  Covariates with index=1 are to be penalized and 0 to be unpenalized. }
  \item{pfactor}{ The user specified proportion to specify the minimum of lambdas equal to 
  pfactor*maximum of lambda corresponding to the largest marginal effect size, which is data 
  dependent and specified by TuningGenerator R function.
  Default is \code{0.10}.}
  \item{Model.selection}{ Model selection criterion with extended BIC as default.}
  \item{restriction}{ This argument allows users to restrict that only regression coefficients with positive signs are kept.}
 }

\value{
 A list of 'PenalizedEstimation' with 8 components:
 
  \item{beta}{ estimates of regression coefficients.}
  \item{iter}{ number of iteration reaching convergence.}
  \item{lambda.output}{used lambda for the outputed beta.}
  \item{tau.output}{used tau for the outputed beta.}
  \item{penalty}{ used penalty.}
  \item{family}{ distribution family.}
  \item{gamma}{ tuning parameter for SCAD or MCP.}
  \item{convex.min}{for SCAD and MCP, directly quoted from ncvreg pacakge: the last index for which 
    the objective function is locally convex.  The smallest value of lambda for which the
    objective function is convex is therefore \code{lambda[convex.min]},
    with corresponding coefficients \code{beta[,convex.min]} }
}
\references{
  Breheny, P. and Huang, J. (2011) Coordinate descent
  algorithms for nonconvex penalized regression, with applications to
  biological feature selection. Ann. Appl. Statist., 5: 232-253.
}
\author{Ting-Huei Chen}

\examples{
  data(XData)
  dim(XData)
  mode(XData)
  data(XIndex)
  data(Y)
  out <- PEN( X=XData,  y=Y, family=c("binomial"), penalty =c("SICA"), ChooseXindex=XIndex)
  out$beta[which(out$beta!=0)]    
}
