\name{simul.ER}
\alias{simul.ER}
\title{
Simulation for Erdos and Renyi,A. (ER) model
}
\description{
		Generate data from a DAG following ER model
}
\usage{
simul.ER(p, pE, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
  	a numeric constant, the number of vertices.
}
  \item{pE}{
  	a numeric constant, the probability for an edge inclusion.
}
  \item{n}{
  	a numeric constant, the sample size.
}
}

\value{
a list
	\item{A}{a pxp matrix representing edge weights of the true graph}
	\item{X}{a nxp data matrix}
	\item{G}{a graphNEL object represening the true graph}
	}
\references{
Erdos, P. and Renyi, A. (1960). On the evolution of random graphs. Publications of the Mathematical Institute of the Hungarian Academy of Sciences, 5:17-61.
}
\author{
Min Jin Ha (\email{minjin.ha@gmail.com})
}
\examples{
	p = 100 # number of vertices
pE = 0.02 # edge inclusion probability
n = 30 # sample size
simul=simul.ER(p,pE,n)
}