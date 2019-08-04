\name{mofredsPC}
\alias{mofredsPC}
\title{ Implement Stable PC algorithm for model free estimation of DAG skeletons }
\description{
  Estimate the DAG skeleton based on the moral graph. 
}
\usage{
mofredsPC(vGraphM, variables, VERBOSE = FALSE, 
          condtest = Zhangtest, trans = 'exp', dimRed = NULL, 
          drmethod = 'sir', alpha = 0.05)
}
\arguments{
  \item{vGraphM}{ The moral graph of the DAG }

  \item{variables}{ An n-by-p matrix of observations with the DAG structure }

  \item{VERBOSE}{ If TRUE, report the proportion of Bootstram samples returning NA when using Songtest. }

  \item{condtest}{ the conditional independence test - currently there are two tests: Songtest, Zhangtest }
  
  \item{trans}{ transformation function in Songtests}
  
  \item{dimRed}{ If TRUE, reduce the dimension before the test. }
  
  \item{drmethod}{ The method to reduce the dimension. There are two methods: SIR (sliced inverse regression), SAVE (sliced average variance estimation). }

  \item{alpha}{ the level of significance level used in the test }
}

\value{
A list of the following components:
 \item{res}{ The removed edges at each step.}
 \item{G}{ The estimated DAG skeleton.}
}
\author{ Wei Sun <wsun@bios.unc.edu> }
\keyword{ methods }
