\name{skeletonPENstable}
\alias{skeletonPENstable}

\title{Estimation of the skeleton of a DAG using a modified PC-stable algorithm starting from a Gaussian graphical model}

\description{Estimation of the skeleton of a directed acyclic graph (DAG) from observational data using a modified PC-stable algorithm starting from a Gaussian graphical model}

\usage{
skeletonPENstable(suffStat, indepTest, p, alpha, verbose = FALSE, edgeWeights, fixedEdges = NULL, NAdelete = FALSE)
}

\arguments{
 \item{suffStat}{Sufficient statistics: List containing all necessary
    elements for the conditional independence decisions in the
    function \code{indepTest}.}
  \item{indepTest}{Predefined function for testing conditional independence. The
    function is internally called as \code{indepTest(x,y,S,suffStat)}, and
    tests conditional independence of \code{x} and \code{y} given
    \code{S}. Here, \code{x} and \code{y} are variables, and \code{S} is
    a (possibly empty) vector of variables (all variables are denoted
    by their column numbers
    in the adjacency matrix). \code{suffStat} is a list containing
    all relevant elements for the conditional independence
    decisions. The return value of \code{indepTest} is the p-value of
    the test for conditional independence.}
  \item{p}{Number of variables.}
  \item{alpha}{Significance level for the individual conditional
    independence tests.}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
  \item{edgeWeights}{
    A symmetric matrix of dimension p*p representing edge weights for a Gaussian graphical model. The modified PC-algorithm starts from the undirected graph made by this matrix.}
  \item{fixedEdges}{A symmetric logical matrix of dimension p*p. If entry
    \code{[i,j]} and \code{[j,i]} are TRUE, the edge i-j is
    never considered for removal. Therefore, this edge is
    guaranteed to be present in the resulting graph.}
 \item{NAdelete}{If indepTest returns \code{NA} and this option is
    \code{TRUE}, the corresponding edge is deleted. If this option is
    \code{FALSE}, the edge is not deleted.}
}
\value{An object of \code{\link{class}} \code{"pcAlgo"} (see
  \code{\linkS4class{pcAlgo}}) containing an estimate of the skeleton of
  the underlying DAG,
  the conditioning sets that led to edge removals (sepset) and several
  other parameters.
}
\references{
    Colombo, D. and Maathuis, M. (2012). A modification of the pc algorithm yielding order-independent skeletons. arXiv preprint arXiv:1211.3295.
    
    Ha et al. PenPC: a two-step approach to estimate the skeletons of high dimensional directed acyclic graphs. 
}
\author{
Min Jin Ha (\email{minjin.ha@gmail.com})}
\examples{
alpha = 0.01 # significance level for a partial correlation testing
p = 100 # number of vertices
e = 1 # one edge added in each time step
n = 30 # sample size
simul=simul.BA(p,e,n)

dat = simul$X
indepTest = gaussCItest
suffStat  = list(C = cor(dat), n = n) 

coeff = ne.PEN(dat=dat,nlambda=100,ntau=10,V=1:p,order=TRUE,verbose=TRUE)
edgeWeights = matrix(0,p,p)
edgeWeights[coeff!=0|t(coeff)!=0] =1

fit.penpc  = skeletonPENstable(suffStat, indepTest, as.integer(p), alpha, 
edgeWeights=edgeWeights, verbose=TRUE)
}