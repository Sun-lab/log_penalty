\name{glmNB.IAL}
\alias{glmNB.IAL}
\title{
Fit a negative binomial regression model using iterative adaptive Lasso penalty
}
\description{
Fit a generalized linear model, with varaible selection by
penalized likelihood method. 
}
\usage{
glmNB.IAL(y, X, delta=seq(0, 5.0, by=0.5),  
        tau = c(0.01*c(10,8,6,4,2,1), 0.001*c(5,2,1)), 
        pMax=20, offset=NULL, naPercent=0.4, 
        nTotal=NULL, maxit=20, maxitIAL=20, nReEstimate=pMax, 
        conv=1e-5, scoreTestP=0.05, trace=1)
}
\arguments{
	
  \item{y}{
	a vector of the response variable. }
	
  \item{X}{
	a matrix of covaraite data, each column is one covariate and 
	each row is one observation. Note, do not include intercept in 
	this covariate matrix. }
	
  \item{delta}{
	one of the two penalization parameters for iterative adpative Lasso }
	
  \item{tau}{
	one of the two penalization parameters for iterative adpative Lasso }

  \item{pMax}{
	the maiximum number of covariates to be kept in the model }
	
  \item{offset}{
	a vector specifying offset in the linear predictor scale, eta = Xbeta + offset.}
	
  \item{naPercent}{
	the maximum percent of missing that are allowed in the data. }
	
  \item{nTotal}{
	if faimly is "binomail", nTotal specify the total number of trails 
	in each subject}
	
  \item{maxit}{
	the maximum number of iteration of Iterated Reweighted Least Square }
	 
  \item{maxitIAL}{
	the maximum number of iteration of for IAL }

  \item{nReEstimate}{
	If there is less than nReEstimate non-zero coefficients, re-estimate these 
  coefficients by least squares in IAL function. }

  \item{conv}{
	convergence tolerence. Let wss be the weighted sum of perason residual square, 
	the algorithm converges if |wss-wss_last|/wss_last < conv }
		
  \item{scoreTestP}{ p-value cutoff for score test of the overdispersion in 
  Poission model. For each associaiton regresion, we first fit a Poisson model, 
  and move to negative binomial if the score test p-value smaller than scoreTestP}

  \item{trace}{ Print out tracking information if trace > 1.}
}

\value{
a list of the following components:
  \item{n2use }{number of covariates selected}

  \item{w2use }{the selected covariates}

  \item{b2use }{the coefficients of the selected covariates}

  \item{b0 }{the intercepts, one for each tuning parameter combination}

  \item{score2use }{the smallest (BIC) score}
  
  \item{delta2use }{delta value corresponding to the smallest score}
  
  \item{score }{(BIC) scores across all the tunning parameters}
  
  \item{nIter }{number of iterations}
  
  \item{family }{Poissoin or Negative Binomial}

  \item{phi }{overdispersion parameter}

  \item{likelihood }{likelihoods across iterations}
}
\author{
Wei Sun weisun@email.unc.edu
}

\keyword{ methods }
