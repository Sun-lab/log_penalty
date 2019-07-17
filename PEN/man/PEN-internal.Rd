\name{PEN-internal}
\title{Internal PEN functions}
\alias{BIC_glm}
\alias{BIC_lm}
\alias{ExtendedBIC_glm}
\alias{ExtendedBIC_lm}
\alias{Create_Updating_Order}
\description{Internal PEN functions}
\usage{
convexMin(beta,X,penalty,gamma,l2,family)
}
\author{Ting-Huei Chen <tingstat22@gmail.com>}
\details{These functions are internal functions for PEN function except for MajorityVote_PEN.
  \code{convexMin} calculates the lowest index for which the penalized objective function
  is locally convex.  \code{TuningGenerator} creates tuning parameter input values.  
  \code{BIC_glm} calculates the value of the BIC for logistic regression.
  \code{BIC_lm} calculates the value of the BIC for linear regression.
  \code{ExtendedBIC_glm} calculates the value of the extended BIC for logistic regression.
  \code{ExtendedBIC_lm} calculates the value of the extended BIC for linear regression.
  \code{Create_Updating_Order} creating random updating orders.	
  \code{MajorityVote_PEN} Runnning PEN multiple times using multiple random updating orders, and the X names that are chosen for more than 50% times will be outputted. 
  }
\keyword{internal}
