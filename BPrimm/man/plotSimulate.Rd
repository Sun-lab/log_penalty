\name{plotSimulate}
\alias{plotSimulate}
\title{ plot simulation result }
\description{
  plot simulation result
}
\usage{
plotSimulate(chr, pos, b.true, main, b, ylab="Coefficient",
 b.cut = 0.05, lty=1, ylim=c(-1,1), mar=c(3,4,3,1), chr.lab=1)
}
\arguments{
  \item{chr}{ chromosome of all the markers }
  \item{pos}{ positions of all the markers }
  \item{b.true}{ the true marker effect, a vector of the same length as chr and pos }
  \item{main}{ title of the plot }
  \item{b}{ the estimated marker effect }
  \item{ylab}{ a title for the y axis }
  \item{b.cut}{ if b<b.cut, it is only ploted as vertical line, but do not have a point at the end of the line. }
  \item{lty}{ line type for the estimated marker effect }
  \item{ylim}{ range of y-axis }
  \item{mar}{ margin, used in par() }
  \item{chr.lab}{ Specifiying the way to draw chromosome labels }
}
\author{ Wei Sun }

\keyword{ utilities }
