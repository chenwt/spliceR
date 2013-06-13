\name{CDSSet}
\alias{CDSSet}
\alias{CDSSet-class}
\docType{class}
\title{
Container for coding sequence (CDS) annotation information
}
\description{
A container for coding sequence annotation information.
}
\usage{
CDSSet(cds)
}
\arguments{
	\item{cds}{A \code{data.frame} object contaning CDS annotation. See details for required columns.}
}
\details{
This object can be generated automatically from \code{\link{getCDS}}, or can be be generated manually by creating a new CDSSet from a data.frame with the following columns: 

\code{chrom}, the chromosome name (NB: chromosome names must match when running \code{\link{annotatePTC}}).
\code{strand}, the strand,
\code{cdsStart}, the genomic start of the coding sequence (beware of 0/1-frame issuses), and
\code{cdsEnd}, the genomic end of the coding sequence (beware of 0/1-frame issuses).

The \code{CDSset} object is required by \code{\link{annotatePTC}} for translation of transcripts and premature termination codon annotation.

For an example, see \code{\link{annotatePTC}}.

}

\value{
A \code{CDSSet} object.
}

\author{Kristoffer Knudsen, Johannes Waage}

\references{
Knudsen, K, Waage, J (2013).
spliceR, an R-package for splice class and PTC annotation of RNA-seq data
\emph{Journal} xx,xx
}

\keyword{classes}
