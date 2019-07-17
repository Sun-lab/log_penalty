\name{connectedComp}
\alias{connectedComp}
\title{
Calculation of connected components of a set of vertices in an undirected graph}
\description{
This function calculates the vertices connected to a set of vertices by any length of chains in an undirected graph}
\usage{
connectedComp(G, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
A p x p connection matrix for an undirected graph. It must be a symmetric square matrix which has the index numbers from 1 to p.}
  \item{x}{
a numeric vector representing a set of vertices among 1,...,p. It is a vector corresponding to the index numbers of G.}
}
\details{
The set of vertices x must be corresponding to the column/row index numbers}
\value{
The length p logical vector corresponding to the index order of G with TRUE if the vertex is a connected component and FALSE otherwise}
\author{
Min Jin Ha (\email{minjin.ha@gmail.com})}
\examples{
# simulate a DAG following ER model
G = simul.ER(p=100,pE=0.02,n=30)$G
Gmat = as(G,"matrix")
Gmat[Gmat!=0|t(Gmat)!=0] = 1 # skeleton of the DAG G
(1:nrow(Gmat))[connectedComp(G=Gmat,x=c(1,2))] # indices for connected components of c(1,2)
}
