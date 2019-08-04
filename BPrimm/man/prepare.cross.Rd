\name{prepare.cross}
\alias{prepare.cross}
\title{ prepare data for function read.cross in R/qtl  }
\description{
  write data into a file which can be readed by read.cross of
  R/qtl using ``csv'' format
}
\usage{
prepare.cross(y, X, chr, pos=NULL, pheno.name, marker.name, file, sep = "\t")
}
\arguments{
  \item{y}{ trait value, could be a vector or a matrix.
   If it is a matrix, each row of y corresponds to one trait }
  \item{X}{ A numerical marix of the marker genotype data,
  each row corresponds to one sample/individual and
  each column corresponds to one marker }
  \item{chr}{ chromosomes of each marker }
  \item{pos}{ positions of each marker }
  \item{pheno.name}{ names of phenotype varialbes (traits) }
  \item{marker.name}{ names of genetic markers }
  \item{file}{ name of the output file }
  \item{sep}{ seperation between two columns }
}

\author{ Wei Sun }
\keyword{ utilities }
