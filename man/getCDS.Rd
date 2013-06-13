\name{getCDS}
\alias{getCDS}
\title{Retrieve CDS information from UCSC}
\description{
Retrieve CDS information from a selected repository from UCSC genome browser repositories.
}
\usage{
getCDS(selectedGenome, repoName)
}
\arguments{
  \item{selectedGenome}{A character, giving the genome. Currently supported are "hg19" and "mm9".} 
  \item{repoName}{A character, giving the gene model repository. Currently supported are "ensemble", "UCSC" (knownGene), and "refseq".} 
 }
\details{
For other genomes and/or gene model repositories, please construct a \code{getCDS} directly (see \code{\link{CDSSet}}).

For a full example of how to use \code{getCDS} in a workflow, please see \code{\link{annotatePTC}}.
}

\value{
A \code{CDSSet} containing the annotated CDSs. For a description of the dataframe, see \code{\link{CDSSet}}. 
}

\author{Kristoffer Knudsen, Johannes Waage}

\references{
Knudsen, K, Waage, J (2013).
spliceR, an R-package for splice class and PTC annotation of RNA-seq data
\emph{Journal} xx,xx
}

\examples{
\dontrun{
mm9UCSC <- getCDS("mm9", "UCSC")
}
}
