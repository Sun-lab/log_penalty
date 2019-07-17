\name{simul.BA}
\alias{simul.BA}
\title{
Simulation for Barabasi and Albert (BA) model
}
\description{
	Generate data from a DAG following BA model
}
\usage{
simul.BA(p, e, n)
}
\arguments{
  \item{p}{
  	a numeric constant, the number of vertices.
}
  \item{e}{
  	a numeric constant, the number of edges to add in each time step.
}
  \item{n}{
  	a numeric constant, the sample size
}
}
\value{
	a list
	\item{A}{a pxp matrix representing edge weights of the true graph}
	\item{X}{a nxp data matrix}
	\item{G}{a graphNEL object represening the true graph}
}
\references{
	Barabasi, A. and Albert, R. (1999). Emergence of scaling in random networks. Science, 286(5439):509-512.
}
\author{
Min Jin Ha (\email{minjin.ha@gmail.com})}
\examples{
	p = 100 # number of vertices
    e = 1 # one edge added in each time step
    n = 30 # sample size
    simul=simul.BA(p,e,n)
    }
